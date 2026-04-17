//! On-disk k-mer index file (`.kmi`) — writer + mmap-backed reader.
//!
//! Format spec: `ferritin-search/docs/KMER_INDEX_FORMAT_SPEC.md`.
//!
//! This module handles only serialization of the same CSR-style layout
//! `KmerIndex` uses in-memory (offsets + entries) to/from a file. The
//! reader mmaps the file and exposes `offsets`, `entries_seq_id`, and
//! `entries_pos` as byte-cast slices into the mmap — no decode, no
//! copy. Peak resident memory for a UniRef50-scale index drops from
//! ~100 GB (in-memory) to ~file-size bounded by the OS page cache.
//!
//! Phase 3 of the ferritin-search memmap refactor — see commit history
//! for phases 1 (memmap DBReader) and 2 (SearchEngine drops targets_full).

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;
use thiserror::Error;

use crate::kmer::{KmerEncoder, KmerHit, KmerIndex};

pub const KMI_MAGIC: [u8; 4] = *b"FKMI";
pub const KMI_VERSION: u32 = 1;
/// Header is fixed at 64 bytes; the three `*_pos` fields are absolute
/// byte offsets into the file and their consistency is validated on
/// open.
pub const KMI_HEADER_SIZE: u64 = 64;

#[derive(Debug, Error)]
pub enum KmiWriterError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Error)]
pub enum KmiReaderError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("file too small to contain a .kmi header (got {0} bytes, need ≥ {1})")]
    TruncatedHeader(u64, u64),
    #[error("bad magic: expected FKMI, got {0:?}")]
    BadMagic([u8; 4]),
    #[error("unsupported .kmi version {found} (reader handles {supported})")]
    BadVersion { found: u32, supported: u32 },
    #[error(
        "layout mismatch in .kmi header: {field} = {found}, expected {expected}"
    )]
    LayoutMismatch {
        field: &'static str,
        found: u64,
        expected: u64,
    },
    #[error(
        "offsets tail invariant failed: offsets[table_size] = {found}, \
         expected n_entries = {expected}"
    )]
    OffsetsTail { found: u64, expected: u64 },
    #[error("file size {found} disagrees with layout (expected {expected})")]
    FileSizeMismatch { found: u64, expected: u64 },
}

/// Write an in-memory [`KmerIndex`] to a `.kmi` file at `path`.
///
/// Pass-through of the CSR arrays: first the `offsets` (u64 LE), then
/// two parallel SoA arrays for the entries — `seq_id` (u32 LE) and
/// `pos` (u16 LE). See the format spec for byte layout.
pub fn write_kmi(index: &KmerIndex, path: impl AsRef<Path>) -> Result<(), KmiWriterError> {
    let path = path.as_ref();
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    let alphabet_size = index.encoder.alphabet_size();
    let kmer_size = index.encoder.kmer_size();
    let table_size = index.encoder.table_size();
    let n_entries = index.entries.len() as u64;

    let offsets_byte_pos = KMI_HEADER_SIZE;
    let entries_seq_id_pos = offsets_byte_pos + 8 * (table_size + 1);
    let entries_pos_pos = entries_seq_id_pos + 4 * n_entries;

    // --- header (64 bytes) ---
    w.write_all(&KMI_MAGIC)?;
    w.write_all(&KMI_VERSION.to_le_bytes())?;
    w.write_all(&alphabet_size.to_le_bytes())?;
    w.write_all(&(kmer_size as u32).to_le_bytes())?;
    w.write_all(&table_size.to_le_bytes())?;
    w.write_all(&n_entries.to_le_bytes())?;
    w.write_all(&offsets_byte_pos.to_le_bytes())?;
    w.write_all(&entries_seq_id_pos.to_le_bytes())?;
    w.write_all(&entries_pos_pos.to_le_bytes())?;
    w.write_all(&[0u8; 8])?; // reserved

    // --- offsets ---
    for off in &index.offsets {
        w.write_all(&off.to_le_bytes())?;
    }

    // --- entries_seq_id ---
    for hit in &index.entries {
        w.write_all(&hit.seq_id.to_le_bytes())?;
    }

    // --- entries_pos ---
    for hit in &index.entries {
        w.write_all(&hit.pos.to_le_bytes())?;
    }

    w.flush()?;
    Ok(())
}

/// Memory-mapped view of a `.kmi` file.
///
/// Slices (`offsets()`, `entries_seq_id()`, `entries_pos()`) are
/// byte-cast directly into the mmap. Lookup by k-mer hash: see
/// [`KmerIndexFile::lookup_hash`].
pub struct KmerIndexFile {
    mmap: Mmap,
    alphabet_size: u32,
    kmer_size: usize,
    table_size: u64,
    n_entries: u64,
    offsets_byte_pos: usize,
    entries_seq_id_pos: usize,
    entries_pos_pos: usize,
}

impl KmerIndexFile {
    /// Open and validate a `.kmi` file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, KmiReaderError> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let file_size = file.metadata()?.len();

        if file_size < KMI_HEADER_SIZE {
            return Err(KmiReaderError::TruncatedHeader(file_size, KMI_HEADER_SIZE));
        }

        // SAFETY: the data file is not mutated for the lifetime of the
        // returned mmap; callers are expected to treat the file as
        // read-only on disk.
        let mmap = unsafe { Mmap::map(&file)? };

        let hdr = &mmap[..KMI_HEADER_SIZE as usize];
        let magic = [hdr[0], hdr[1], hdr[2], hdr[3]];
        if magic != KMI_MAGIC {
            return Err(KmiReaderError::BadMagic(magic));
        }
        let version = u32_at(hdr, 4);
        if version != KMI_VERSION {
            return Err(KmiReaderError::BadVersion {
                found: version,
                supported: KMI_VERSION,
            });
        }
        let alphabet_size = u32_at(hdr, 8);
        let kmer_size = u32_at(hdr, 12) as usize;
        let table_size = u64_at(hdr, 16);
        let n_entries = u64_at(hdr, 24);
        let offsets_byte_pos = u64_at(hdr, 32);
        let entries_seq_id_pos = u64_at(hdr, 40);
        let entries_pos_pos = u64_at(hdr, 48);

        let expected_offsets = KMI_HEADER_SIZE;
        if offsets_byte_pos != expected_offsets {
            return Err(KmiReaderError::LayoutMismatch {
                field: "offsets_byte_pos",
                found: offsets_byte_pos,
                expected: expected_offsets,
            });
        }
        let expected_seq_id_pos = offsets_byte_pos + 8 * (table_size + 1);
        if entries_seq_id_pos != expected_seq_id_pos {
            return Err(KmiReaderError::LayoutMismatch {
                field: "entries_seq_id_pos",
                found: entries_seq_id_pos,
                expected: expected_seq_id_pos,
            });
        }
        let expected_pos_pos = entries_seq_id_pos + 4 * n_entries;
        if entries_pos_pos != expected_pos_pos {
            return Err(KmiReaderError::LayoutMismatch {
                field: "entries_pos_pos",
                found: entries_pos_pos,
                expected: expected_pos_pos,
            });
        }
        let expected_file_size = entries_pos_pos + 2 * n_entries;
        if file_size != expected_file_size {
            return Err(KmiReaderError::FileSizeMismatch {
                found: file_size,
                expected: expected_file_size,
            });
        }

        // Offsets tail invariant: the last offset must equal n_entries.
        let offsets_start = offsets_byte_pos as usize;
        let last_off_pos = offsets_start + 8 * table_size as usize;
        let last_off = u64_at(&mmap[..], last_off_pos);
        if last_off != n_entries {
            return Err(KmiReaderError::OffsetsTail {
                found: last_off,
                expected: n_entries,
            });
        }

        Ok(Self {
            mmap,
            alphabet_size,
            kmer_size,
            table_size,
            n_entries,
            offsets_byte_pos: offsets_byte_pos as usize,
            entries_seq_id_pos: entries_seq_id_pos as usize,
            entries_pos_pos: entries_pos_pos as usize,
        })
    }

    pub fn alphabet_size(&self) -> u32 { self.alphabet_size }
    pub fn kmer_size(&self) -> usize { self.kmer_size }
    pub fn table_size(&self) -> u64 { self.table_size }
    pub fn n_entries(&self) -> u64 { self.n_entries }

    /// Reconstruct a `KmerEncoder` matching the file's parameters.
    pub fn encoder(&self) -> KmerEncoder {
        KmerEncoder::new(self.alphabet_size, self.kmer_size)
    }

    /// Byte slice over the `(table_size + 1) × u64 LE` offsets array.
    ///
    /// Returns a `&[u8]` rather than `&[u64]` to sidestep alignment
    /// concerns on platforms where mmap pages may not align to 8-byte
    /// boundaries. Callers decode per lookup via `u64_at`.
    fn offsets_bytes(&self) -> &[u8] {
        let start = self.offsets_byte_pos;
        let len = 8 * (self.table_size as usize + 1);
        &self.mmap[start..start + len]
    }

    fn entries_seq_id_bytes(&self) -> &[u8] {
        let start = self.entries_seq_id_pos;
        let len = 4 * self.n_entries as usize;
        &self.mmap[start..start + len]
    }

    fn entries_pos_bytes(&self) -> &[u8] {
        let start = self.entries_pos_pos;
        let len = 2 * self.n_entries as usize;
        &self.mmap[start..start + len]
    }

    /// Look up the `(seq_id, pos)` postings for a single k-mer hash.
    ///
    /// Returns a `Vec<KmerHit>` — the SoA layout on disk forces one
    /// small owned allocation per lookup. For hot prefilter loops the
    /// cost is tiny relative to the OS page-fault of the posting
    /// bucket itself, and the API stays compatible with in-memory
    /// `KmerIndex::lookup_hash`.
    pub fn lookup_hash(&self, hash: u64) -> Vec<KmerHit> {
        debug_assert!(hash <= self.table_size);
        let offsets = self.offsets_bytes();
        let start = u64_at(offsets, 8 * hash as usize) as usize;
        let end = u64_at(offsets, 8 * (hash as usize + 1)) as usize;
        if end <= start {
            return Vec::new();
        }
        let seq_ids = self.entries_seq_id_bytes();
        let positions = self.entries_pos_bytes();
        (start..end)
            .map(|i| KmerHit {
                seq_id: u32_at(seq_ids, 4 * i),
                pos: u16_at(positions, 2 * i),
            })
            .collect()
    }
}

// ---- internal LE readers over &[u8] ----

fn u16_at(buf: &[u8], off: usize) -> u16 {
    u16::from_le_bytes([buf[off], buf[off + 1]])
}
fn u32_at(buf: &[u8], off: usize) -> u32 {
    u32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]])
}
fn u64_at(buf: &[u8], off: usize) -> u64 {
    u64::from_le_bytes([
        buf[off], buf[off + 1], buf[off + 2], buf[off + 3],
        buf[off + 4], buf[off + 5], buf[off + 6], buf[off + 7],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{KmerEncoder, KmerIndex};
    use tempfile::tempdir;

    fn tiny_index() -> KmerIndex {
        // 13-letter reduced alphabet, k=3, three tiny targets.
        let encoder = KmerEncoder::new(13, 3);
        let targets: Vec<(u32, &[u8])> = vec![
            (1u32, &[0, 1, 2, 3, 4][..]),
            (2u32, &[2, 3, 4, 5, 6][..]),
            (3u32, &[0, 1, 2][..]),
        ];
        KmerIndex::build(encoder, targets.into_iter(), 12).unwrap()
    }

    #[test]
    fn round_trip_recovers_offsets_and_entries() {
        let idx = tiny_index();
        let dir = tempdir().unwrap();
        let path = dir.path().join("tiny.kmi");
        write_kmi(&idx, &path).unwrap();

        let file = KmerIndexFile::open(&path).unwrap();
        assert_eq!(file.alphabet_size(), idx.encoder.alphabet_size());
        assert_eq!(file.kmer_size(), idx.encoder.kmer_size());
        assert_eq!(file.table_size(), idx.encoder.table_size());
        assert_eq!(file.n_entries() as usize, idx.entries.len());

        // Every kmer hash yields byte-identical hits between the
        // in-memory index and the mmap-backed file.
        for h in 0..idx.encoder.table_size() {
            let mem: Vec<_> = idx.lookup_hash(h).to_vec();
            let disk = file.lookup_hash(h);
            assert_eq!(
                mem, disk,
                "kmer hash {h}: in-memory and on-disk lookups diverge",
            );
        }
    }

    #[test]
    fn open_rejects_bad_magic() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("nope.kmi");
        std::fs::write(&path, [0u8; 128]).unwrap();
        let err = KmerIndexFile::open(&path).err().expect("expected BadMagic");
        assert!(
            matches!(err, KmiReaderError::BadMagic(m) if m == [0, 0, 0, 0]),
            "expected BadMagic, got {err}",
        );
    }

    #[test]
    fn open_rejects_future_version() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("v99.kmi");
        let mut hdr = [0u8; 64];
        hdr[..4].copy_from_slice(&KMI_MAGIC);
        hdr[4..8].copy_from_slice(&99u32.to_le_bytes());
        std::fs::write(&path, hdr).unwrap();
        let err = KmerIndexFile::open(&path).err().expect("expected BadVersion");
        assert!(
            matches!(err, KmiReaderError::BadVersion { found: 99, supported: 1 }),
            "expected BadVersion(99), got {err}",
        );
    }

    #[test]
    fn open_rejects_truncated_file() {
        let idx = tiny_index();
        let dir = tempdir().unwrap();
        let path = dir.path().join("short.kmi");
        write_kmi(&idx, &path).unwrap();
        // Truncate a few bytes off the end.
        let full = std::fs::read(&path).unwrap();
        std::fs::write(&path, &full[..full.len() - 4]).unwrap();
        let err = KmerIndexFile::open(&path).err().expect("expected truncation error");
        assert!(
            matches!(err, KmiReaderError::FileSizeMismatch { .. }),
            "expected FileSizeMismatch, got {err}",
        );
    }

    #[test]
    fn empty_buckets_return_empty_vec() {
        // Hash values with zero hits must return an empty slice, not
        // panic on bounds and not return stale data from a neighbor.
        let idx = tiny_index();
        let dir = tempdir().unwrap();
        let path = dir.path().join("tiny.kmi");
        write_kmi(&idx, &path).unwrap();
        let file = KmerIndexFile::open(&path).unwrap();
        // Some hash slots are guaranteed empty at this corpus size
        // (table_size = 13^3 = 2197 slots, only a few are populated).
        let mut empty_seen = 0;
        for h in 0..file.table_size() {
            if file.lookup_hash(h).is_empty() {
                empty_seen += 1;
            }
        }
        assert!(empty_seen > 0, "expected some empty k-mer buckets in a tiny corpus");
    }
}
