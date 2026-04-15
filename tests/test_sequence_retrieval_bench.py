from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "validation" / "bench_sequence_retrieval.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("bench_sequence_retrieval", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_kmer_counts_and_weighted_jaccard():
    bench = _load_module()

    left = bench.kmer_counts("AAAA", 2)
    right = bench.kmer_counts("AAAT", 2)

    assert left == {"AA": 3}
    assert right == {"AA": 2, "AT": 1}
    assert bench.weighted_jaccard(left, right) == 0.5


def test_smith_waterman_normalized_score_orders_similarity():
    bench = _load_module()

    self_score = bench.smith_waterman_score("ACDEFG", "ACDEFG", match=2, mismatch=-1, gap=-2)
    weak_score = bench.smith_waterman_score("ACDEFG", "WWWWWW", match=2, mismatch=-1, gap=-2)

    assert self_score == 1.0
    assert weak_score < self_score


def test_sanitize_sequence_for_hmmer_normalizes_noncanonical_symbols():
    bench = _load_module()

    assert bench.sanitize_sequence_for_hmmer("ACDUOBZ*") == "ACDXXXXX"


def test_rank_sequence_hits_reranks_prefilter_candidates():
    bench = _load_module()
    paths = [Path("q.pdb"), Path("same.pdb"), Path("weak.pdb")]
    sequences = {
        "q.pdb": "ACDEFGHIK",
        "same.pdb": "ACDEFGHIK",
        "weak.pdb": "WWWWWWWWW",
    }
    kmer_index = {key: bench.kmer_counts(value, 3) for key, value in sequences.items()}

    hits = bench.rank_sequence_hits(
        Path("q.pdb"),
        paths,
        sequences,
        kmer_index,
        kmer_prefilter_top_k=3,
        sw_rerank_top_k=3,
        sw_match=2,
        sw_mismatch=-1,
        sw_gap=-2,
    )

    assert hits[0]["source_path"] == "q.pdb"
    assert hits[1]["source_path"] == "same.pdb"
    assert hits[1]["smith_waterman_score"] == 1.0


def test_parse_hmmer_tblout(tmp_path: Path):
    bench = _load_module()
    tblout = tmp_path / "hits.tblout"
    tblout.write_text(
        "\n".join([
            "# target name accession query name accession E-value score bias",
            "seq1 - query - 1e-20 75.5 0.0 1e-20 75.5 0.0 1.0 1 1 1 1 1 1 1 hit one",
            "missing - query - 2e-10 40.0 0.0 2e-10 40.0 0.0 1.0 1 1 1 1 1 1 1 ignored",
        ]),
        encoding="utf-8",
    )

    hits = bench.parse_hmmer_tblout(
        tblout,
        {"seq1": Path("target.pdb")},
        {"target.pdb": "ACDEFG"},
    )

    assert hits == [
        {
            "source_path": "target.pdb",
            "sequence_score": 75.5,
            "hmmer_bits": 75.5,
            "hmmer_evalue": 1e-20,
            "query_length": None,
            "target_length": 6,
        }
    ]


def test_sequence_cache_roundtrip(tmp_path: Path):
    bench = _load_module()
    cache = tmp_path / "seqs.json"
    paths = [Path("a.pdb"), Path("b.pdb")]

    bench.save_sequence_cache(
        cache,
        target_paths=paths,
        sequences={"a.pdb": "ACD"},
        skipped_sequences={"b.pdb": "sequence_too_long:100001"},
    )

    assert bench.load_sequence_cache(cache, paths) == (
        {"a.pdb": "ACD"},
        {"b.pdb": "sequence_too_long:100001"},
    )
    assert bench.load_sequence_cache(cache, [Path("other.pdb")]) is None
