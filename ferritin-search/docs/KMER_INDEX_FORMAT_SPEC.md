# Ferritin-native on-disk k-mer index format (`.kmi`)

**Last updated: 2026-04-17**
**Format version: 1**
**Magic: `FKMI` (`0x46 0x4B 0x4D 0x49`)**

This spec defines a ferritin-native on-disk format for the k-mer index
used by `ferritin-search::search::SearchEngine`. The format exists
because the in-memory `KmerIndex` is the dominant RAM consumer at
archive-scale corpora:

- UniRef50 (~50M sequences × avg 245 k-mer starts at k=6) produces
  ~12.3 × 10⁹ (seq_id, pos) entries. At 8 bytes in-memory that's 98 GB
  resident before the engine has done any other work.

The `.kmi` format lays the same CSR-style `(offsets, entries)` structure
out on disk as mmap-ready byte arrays. The engine opens the file, points
two slices into the mmap, and lets the OS page in only the k-mer
postings a given query's prefilter actually touches.

Design decisions:

1. **Not byte-compatible with upstream `mmseqs createindex`.** Upstream
   index format is versioned, platform-dependent, SIMD-layout-tuned, and
   carries distributed-search artifacts we don't need. A ferritin-native
   format keeps the spec small and under our control.
2. **No compression in v1.** Offsets and entries are dense integer
   arrays; zstd on the whole file is a possible v2 layer but reader
   complexity goes up a lot (streaming vs mmap).
3. **Little-endian all the way.** Matches `.dbtype` convention and every
   host ferritin runs on.
4. **Reader is zero-copy.** After mmap + header parse, `offsets` and
   `entries` are bit-cast slices into the mmap — no decode, no copy.
5. **Writer is one-pass over an in-memory `KmerIndex`.** We do not build
   the on-disk index directly from sequences yet; instead we first
   produce an in-memory index (as today), then dump it to disk. A
   later phase may add a two-pass external-memory builder that never
   materializes the full in-memory index, but v1 keeps the write path
   trivial to audit.

## File layout

```
Offset  Size  Field                 Value
------  ----  --------------------  ------------------------------------
0       4     magic                 "FKMI"
4       4     version               u32 LE  (current: 1)
8       4     alphabet_size         u32 LE  (e.g. 21 protein, 13 reduced)
12      4     kmer_size             u32 LE  (k; typically 6)
16      8     table_size            u64 LE  (= alphabet_size^kmer_size)
24      8     n_entries             u64 LE  (total k-mer hits)
32      8     offsets_byte_pos      u64 LE  (= HEADER_SIZE, see below)
40      8     entries_seq_id_pos    u64 LE
48      8     entries_pos_pos       u64 LE
56      8     reserved              0x00 × 8  (pad to 64-byte header)
64      ...   offsets               (table_size + 1) × u64 LE
...     ...   entries_seq_id        n_entries × u32 LE
...     ...   entries_pos           n_entries × u16 LE
```

Header size is fixed at **64 bytes**. The three `*_pos` fields are
absolute file offsets that a reader must cross-check against the
layout's implied positions — they exist so validators can diagnose
truncation without reconstructing the arithmetic.

## On-disk entries: struct-of-arrays

The in-memory `KmerHit` is `{ seq_id: u32, pos: u16 }`; AoS on disk
would require either packed structs (Rust `#[repr(C, packed)]` — hard
to borrow cleanly from an mmap without unaligned access) or a fixed
8-byte wasteful padding. SoA — two parallel arrays of `u32` seq_ids
and `u16` positions — gives the same 6-byte-per-entry storage cost
without alignment headaches. For a hit at index `i`:

```rust
let seq_id = entries_seq_id[i];  // u32 LE, aligned
let pos    = entries_pos[i];     // u16 LE, aligned
```

At 50M UniRef50 targets with ~245 k-mer starts each (k=6), that's
~12.3 × 10⁹ hits × 6 bytes = ~74 GB on disk. The reader mmaps the
whole file; only pages touched by queried k-mers occupy RSS.

Position overflow (`pos > u16::MAX`) is refused at build time exactly
like the in-memory builder does (see `kmer::KmerIndexError::PositionOverflow`).

## `offsets` semantics

`offsets[k..k+1]` delimit the hits for k-mer hash `k` in `entries`:

```rust
let start = offsets[k as usize] as usize;
let end = offsets[k as usize + 1] as usize;
let hits: &[KmerHitOnDisk] = &entries[start..end];
```

`offsets` has `table_size + 1` elements so the last offset doubles as
the total entry count (= `n_entries`). Readers must check
`offsets[table_size] == n_entries` as a sanity guard.

## Validation contract

A reader opens an `.kmi` file and validates, in order:

1. File size ≥ 64 bytes (header).
2. `magic == "FKMI"`.
3. `version == 1` (future versions bump this; a v1 reader refuses v2+).
4. `offsets_byte_pos == 64`.
5. `entries_seq_id_pos == 64 + 8 * (table_size + 1)`.
6. `entries_pos_pos == entries_seq_id_pos + 4 * n_entries`.
7. File size == `entries_pos_pos + 2 * n_entries`.
8. `offsets[table_size] == n_entries`.

Failures raise `KmiReaderError::{BadMagic, BadVersion, TruncatedHeader,
LayoutMismatch, OffsetMismatch}` so callers can distinguish "wrong file
format" from "correct format, wrong version" from "corrupted file."

## Round-trip contract

Building an in-memory `KmerIndex`, writing it to `.kmi`, re-opening via
mmap, and running `lookup_hash(k)` on any k-mer must return byte-identical
postings as the original in-memory lookup. This is the core test in
`ferritin-search/tests/kmer_index_roundtrip.rs`.

## Out of scope for v1

- Compression (zstd, delta encoding on positions within a kmer's posting
  list, etc.). The ~74 GB size on UniRef50 is disk-cheap on monster3;
  compression helps cold-start paging but complicates the reader. Revisit
  if we hit archive-scale corpora where disk pressure matters.
- Distributed indices (shards across multiple files). A single index
  file covers our scale target.
- On-disk similar-k-mer expansion tables. Those are a prefilter-quality
  feature that mmseqs offers; ferritin's prefilter doesn't use them yet.
- Index-of-indices for very sparse tables. At k=6, even the full 21-letter
  alphabet only has 85.8M slots × 8 bytes = 687 MB of offsets, which
  comfortably mmaps.
