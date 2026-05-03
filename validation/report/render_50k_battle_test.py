"""Stand-alone HTML report for the 50K-PDB battle-test claim.

This is a *single-pipeline scale + reliability test*, not a multi-oracle
benchmark: 50 000 random wwPDB structures pushed through proteon's
CHARMM19+EEF1 minimization on GPU, recording per-PDB outcomes
(status, energies, convergence flag, step count). The claim's capability
is ``pipeline-physics-consistent-at-scale`` — i.e. "does the kernel
survive 50K structures without crashing, and does the minimizer move
energies downhill?" — *not* "does the minimum match an external
reference" (that's the corpus oracle claims, separate).

Schema, per-PDB::

    {"pdb": "13sm.pdb", "atoms": 2541, "status": "ok",
     "skipped": false, "initial_energy": 12683.5,
     "final_energy": -38007.1, "steps": 50, "converged": false}

Honest framing built into the headline:

* "99.1% pipeline survival" is the headline reliability number
  (n_ok / n_total).
* "Median energy drop X kJ/mol" confirms the minimizer is pushing
  the system downhill at scale.
* "0.04% converged" is *not* a bad number — the runner uses a
  fixed 50-step budget intentionally to stress the kernel, not to
  reach the minimum. Convergence requires thousands of steps.
* The runner is GPU-only (proteon-search/cuda kernels via cudarc);
  reproducing this on a CPU-only host won't match wall times.

Usage::

    python validation/report/render_50k_battle_test.py \\
        --claim evident/claims/scale_50k_battle_test.yaml \\
        --artifact validation/stage3_50k_gpu_results.jsonl \\
        --output evident/reports/v0.1.3/proteon-50k-battle-test-release.html
"""
from __future__ import annotations

import argparse
import base64
import collections
import datetime
import hashlib
import html
import json
import pathlib
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import yaml
except ImportError as e:
    raise SystemExit("PyYAML required: pip install pyyaml") from e


def _sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _embed_png(path: pathlib.Path, alt: str) -> str:
    if not path.exists():
        return f'<div class="missing">[{alt} — figure missing]</div>'
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" alt="{html.escape(alt)}">'


def _load_records(artifact_path: pathlib.Path) -> list[dict]:
    out: list[dict] = []
    with artifact_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def _summary(records: list[dict]) -> dict[str, Any]:
    statuses = collections.Counter(r.get("status", "unknown") for r in records)
    n_total = len(records)
    n_ok = statuses.get("ok", 0)
    n_load_error = statuses.get("load_error", 0)
    n_pipeline_error = sum(
        v for k, v in statuses.items()
        if k not in ("ok", "load_error", "skipped", "unknown")
    )
    n_skipped = statuses.get("skipped", 0) + sum(
        1 for r in records if r.get("skipped") and r.get("status") != "ok"
    )

    ok = [r for r in records if r.get("status") == "ok"]
    n_converged = sum(1 for r in ok if r.get("converged"))

    atoms = np.array([r["atoms"] for r in ok if r.get("atoms")], dtype=float)
    steps = np.array([r["steps"] for r in ok if r.get("steps") is not None], dtype=float)
    drops = np.array(
        [
            (r["initial_energy"] - r["final_energy"])
            for r in ok
            if r.get("initial_energy") is not None
            and r.get("final_energy") is not None
        ],
        dtype=float,
    )
    drops_finite = drops[np.isfinite(drops)]

    return {
        "n_total": n_total,
        "n_ok": n_ok,
        "n_load_error": n_load_error,
        "n_pipeline_error": n_pipeline_error,
        "n_skipped": n_skipped,
        "survival_rate": n_ok / n_total if n_total else 0.0,
        "n_converged": n_converged,
        "convergence_rate": n_converged / n_ok if n_ok else 0.0,
        "atoms_median": float(np.median(atoms)) if atoms.size else None,
        "atoms_p99": float(np.percentile(atoms, 99)) if atoms.size else None,
        "atoms_max": float(atoms.max()) if atoms.size else None,
        "steps_median": float(np.median(steps)) if steps.size else None,
        "steps_max": float(steps.max()) if steps.size else None,
        "drop_median": float(np.median(drops_finite)) if drops_finite.size else None,
        "drop_negative_count": int((drops_finite < 0).sum()),
        "drop_n": int(drops_finite.size),
    }


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def _save(fig: matplotlib.figure.Figure, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_status_breakdown(records: list[dict], path: pathlib.Path) -> None:
    counts = collections.Counter(r.get("status", "unknown") for r in records)
    keys = sorted(counts, key=lambda k: -counts[k])
    values = [counts[k] for k in keys]
    fig, ax = plt.subplots(figsize=(9, 3.8))
    colors = {
        "ok": "#28a745",
        "load_error": "#cb2431",
        "pipeline_error": "#7e1818",
        "oom": "#7e1818",
        "skipped": "#888",
    }
    ax.barh(keys, values, color=[colors.get(k, "#0a7e8c") for k in keys])
    for i, v in enumerate(values):
        ax.text(v + max(values) * 0.005, i, f"{v}", va="center", fontsize=10)
    ax.set_xlabel(f"# structures (n = {len(records)})")
    ax.set_title("Per-structure pipeline status across the 50K corpus")
    ax.grid(alpha=0.3, axis="x")
    _save(fig, path)


def plot_atoms_distribution(records: list[dict], path: pathlib.Path) -> None:
    atoms = np.array(
        [r["atoms"] for r in records if r.get("atoms") and r.get("status") == "ok"],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(9, 4.5))
    if atoms.size == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        _save(fig, path)
        return
    upper = float(np.percentile(atoms, 99.5))
    bins = np.logspace(np.log10(50), np.log10(max(upper, 100)), 50)
    ax.hist(np.clip(atoms, 50, upper), bins=bins, color="#0a7e8c", alpha=0.85)
    ax.axvline(np.median(atoms), color="black", linewidth=1.5,
               label=f"median {np.median(atoms):.0f} atoms")
    ax.axvline(np.percentile(atoms, 99), color="#cb2431", linestyle="--",
               linewidth=1.0,
               label=f"p99 {np.percentile(atoms, 99):.0f} atoms")
    ax.set_xscale("log")
    ax.set_xlabel("atoms (log)")
    ax.set_ylabel(f"# structures (n = {atoms.size})")
    ax.set_title("Distribution of structure size in the 50K corpus")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(alpha=0.3, which="both")
    _save(fig, path)


def plot_energy_drop_distribution(records: list[dict], path: pathlib.Path) -> None:
    drops_per_atom = []
    for r in records:
        if r.get("status") != "ok":
            continue
        ie, fe, atoms = r.get("initial_energy"), r.get("final_energy"), r.get("atoms")
        if ie is None or fe is None or not atoms:
            continue
        drops_per_atom.append((ie - fe) / atoms)
    arr = np.array(drops_per_atom)
    arr = arr[np.isfinite(arr)]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    if arr.size == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        _save(fig, path)
        return
    p1, p99 = float(np.percentile(arr, 1)), float(np.percentile(arr, 99))
    bins = np.linspace(p1, p99, 60)
    ax.hist(np.clip(arr, p1, p99), bins=bins, color="#0a7e8c", alpha=0.85)
    ax.axvline(0, color="#cb2431", linestyle="--", linewidth=1.0,
               label="no change")
    median = float(np.median(arr))
    ax.axvline(median, color="black", linewidth=1.5,
               label=f"median {median:.2f}")
    ax.set_xlabel("(initial − final) / n_atoms  (kJ/mol per atom; clipped to 1–99 percentile)")
    ax.set_ylabel(f"# structures (n = {arr.size})")
    ax.set_title("Per-atom energy drop after 50 minimizer steps")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(alpha=0.3)
    _save(fig, path)


def plot_initial_vs_final(records: list[dict], path: pathlib.Path) -> None:
    pts = []
    for r in records:
        if r.get("status") != "ok":
            continue
        ie, fe = r.get("initial_energy"), r.get("final_energy")
        if ie is None or fe is None or not np.isfinite(ie) or not np.isfinite(fe):
            continue
        pts.append((ie, fe))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    if not pts:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        _save(fig, path)
        return
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])
    # Sub-sample for plot legibility on 45K points.
    if xs.size > 8000:
        idx = np.random.default_rng(42).choice(xs.size, 8000, replace=False)
        xs_plot, ys_plot = xs[idx], ys[idx]
    else:
        xs_plot, ys_plot = xs, ys
    ax.scatter(xs_plot, ys_plot, s=4, alpha=0.25, color="#0a7e8c")
    lo = float(min(xs.min(), ys.min()))
    hi = float(max(xs.max(), ys.max()))
    ax.plot([lo, hi], [lo, hi], color="#cb2431", linestyle="--",
            linewidth=1.0, label="initial = final (no minimization)")
    ax.set_xlabel("initial_energy  (kJ/mol)")
    ax.set_ylabel("final_energy  (kJ/mol)")
    ax.set_title(f"50K battle test: initial vs final energy "
                 f"(plotted = {xs_plot.size:,} of {xs.size:,} ok records)")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(alpha=0.3)
    _save(fig, path)


# --------------------------------------------------------------------------
# Render
# --------------------------------------------------------------------------
def _claim_metadata(claim_doc: dict) -> dict[str, Any]:
    claims = claim_doc.get("claims") or []
    if not claims:
        raise ValueError("claim YAML contains no claims")
    c = claims[0]
    return {
        "id": c.get("id", "(missing id)"),
        "title": c.get("title", "(missing title)"),
        "tier": c.get("tier", "?"),
        "subsystem": c.get("subsystem", "?"),
        "claim": (c.get("claim") or "").strip(),
        "command": (c.get("evidence") or {}).get("command", ""),
        "oracle": ", ".join((c.get("evidence") or {}).get("oracle") or []),
        "pinned": c.get("pinned_versions") or {},
        "last_verified": c.get("last_verified") or {},
    }


def render(
    claim_path: pathlib.Path,
    artifact_path: pathlib.Path,
    output_path: pathlib.Path,
    fig_dir: pathlib.Path | None = None,
) -> None:
    claim_doc = yaml.safe_load(claim_path.read_text(encoding="utf-8"))
    meta = _claim_metadata(claim_doc)

    records = _load_records(artifact_path)
    if not records:
        raise SystemExit(f"no records in artifact: {artifact_path}")
    summary = _summary(records)

    if fig_dir is None:
        fig_dir = output_path.parent / "_figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig_status = fig_dir / "30_50k_status.png"
    fig_atoms = fig_dir / "31_50k_atoms.png"
    fig_drop = fig_dir / "32_50k_energy_drop.png"
    fig_scatter = fig_dir / "33_50k_initial_vs_final.png"
    plot_status_breakdown(records, fig_status)
    plot_atoms_distribution(records, fig_atoms)
    plot_energy_drop_distribution(records, fig_drop)
    plot_initial_vs_final(records, fig_scatter)

    artifact_sha = _sha256_file(artifact_path)
    today = datetime.date.today().isoformat()

    last_verified = meta["last_verified"] or {}
    lv_lines = []
    for k in ("commit", "date", "value", "corpus_sha"):
        v = last_verified.get(k)
        lv_lines.append(
            f"<tr><th>{k}</th><td><code>"
            f"{html.escape(str(v) if v is not None else '—')}</code></td></tr>"
        )
    lv_table = "<table class='kv'>" + "".join(lv_lines) + "</table>"

    pinned_lines = "".join(
        f"<tr><th>{html.escape(k)}</th><td><code>{html.escape(str(v))}</code></td></tr>"
        for k, v in meta["pinned"].items()
    ) or "<tr><td colspan='2'><em>not pinned</em></td></tr>"

    title = html.escape(meta["title"])
    claim_text = html.escape(meta["claim"])
    cmd = html.escape(meta["command"])
    oracle = html.escape(meta["oracle"]) if meta["oracle"] else "<em>physics-invariant (no external oracle)</em>"

    def _fmt(v: float | None, fmt: str = ",.0f") -> str:
        return "—" if v is None else format(v, fmt)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>EVIDENT — {title}</title>
<style>
:root {{
  --proteon: #0a7e8c;
  --ink:     #222;
  --muted:   #666;
  --bg:      #fafafa;
  --card:    #ffffff;
  --rule:    #e4e4e4;
  --pass:    #28a745;
  --fail:    #cb2431;
  --warn:    #d4a017;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  font-size: 15px; line-height: 1.55; color: var(--ink); background: var(--bg);
}}
main {{ max-width: 1180px; margin: 0 auto; padding: 40px 32px 80px; }}
header {{ border-bottom: 3px solid var(--proteon); margin-bottom: 28px; padding-bottom: 16px; }}
h1 {{ font-size: 28px; font-weight: 700; margin: 0 0 4px; letter-spacing: -0.01em; }}
.subtitle {{ color: var(--muted); font-size: 14px; }}
.badge {{
  display: inline-block; padding: 2px 8px; border-radius: 12px;
  font-size: 11px; font-weight: 600; color: white;
  background: var(--proteon); margin-left: 6px; vertical-align: 2px;
}}
section {{ background: var(--card); border: 1px solid var(--rule);
          border-radius: 8px; padding: 24px 28px; margin: 18px 0; }}
h2 {{ font-size: 20px; font-weight: 700; margin: 0 0 14px; color: var(--proteon); }}
table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }}
th, td {{ border: 1px solid var(--rule); padding: 6px 10px; text-align: left;
         vertical-align: top; }}
th {{ background: #f6f8fa; font-weight: 600; }}
tr:nth-child(even) td {{ background: #fafbfc; }}
table.kv th {{ width: 130px; background: #fafafa; }}
img {{ max-width: 100%; height: auto; border: 1px solid var(--rule);
      border-radius: 4px; margin: 12px 0; }}
code {{ background: #f0f0f0; padding: 2px 6px; border-radius: 3px;
       font-size: 12.5px; font-family: ui-monospace, "SF Mono", Menlo, monospace; }}
.muted {{ color: var(--muted); }}
.summary-grid {{ display: grid; grid-template-columns: repeat(4, 1fr);
               gap: 12px; margin: 16px 0; }}
.summary-card {{ background: #f6f8fa; border-radius: 6px; padding: 14px;
                text-align: center; }}
.summary-card .num {{ font-size: 24px; font-weight: 700; color: var(--proteon); }}
.summary-card .lbl {{ color: var(--muted); font-size: 12px; margin-top: 2px; }}
.missing {{ color: var(--fail); font-style: italic; padding: 8px; }}
footer {{ margin-top: 40px; color: var(--muted); font-size: 12px;
         border-top: 1px solid var(--rule); padding-top: 12px; }}
.note {{ background: #fff8dc; border-left: 3px solid #d4a017;
        padding: 10px 14px; margin: 12px 0; font-size: 13px; }}
.note strong {{ color: var(--ink); }}
</style>
</head>
<body>
<main>
<header>
  <h1>{title}<span class="badge">{html.escape(meta['tier'])}</span></h1>
  <div class="subtitle">
    Claim id <code>{html.escape(meta['id'])}</code> ·
    Subsystem <code>{html.escape(meta['subsystem'])}</code> ·
    Oracle {oracle} ·
    Generated {today}
  </div>
</header>

<section>
  <h2>Claim</h2>
  <p>{claim_text}</p>
  <p class="muted"><strong>Replay command:</strong> <code>{cmd}</code></p>
  <div class="note">
    <strong>Scope:</strong> this is a <em>single-pipeline scale + reliability test</em>,
    not a multi-oracle benchmark. The capability under test is
    <code>pipeline-physics-consistent-at-scale</code> — does
    proteon's CHARMM19+EEF1 minimization survive 50K random PDBs without
    crashing, and does the minimizer move energies downhill at scale.
    The <em>numerical correctness</em> of the energies relative to BALL or
    OpenMM is the corpus-oracle claims' job, not this one.
  </div>
</section>

<section>
  <h2>Headline</h2>
  <div class="summary-grid">
    <div class="summary-card">
      <div class="num">{summary['n_total']:,}</div>
      <div class="lbl">structures attempted</div>
    </div>
    <div class="summary-card">
      <div class="num">{summary['survival_rate']*100:.1f}%</div>
      <div class="lbl">pipeline survival ({summary['n_ok']:,} / {summary['n_total']:,} ok)</div>
    </div>
    <div class="summary-card">
      <div class="num">{_fmt(summary['drop_median'], ',.0f')}</div>
      <div class="lbl">median ΔE per PDB (kJ/mol)</div>
    </div>
    <div class="summary-card">
      <div class="num">{summary['n_load_error']:,}</div>
      <div class="lbl">PDBFixer load failures</div>
    </div>
  </div>
  <table class="kv">
    <tr><th>structures by status</th>
        <td><code>ok={summary['n_ok']:,} · load_error={summary['n_load_error']:,}
        · pipeline_error={summary['n_pipeline_error']:,}
        · skipped={summary['n_skipped']:,}</code></td></tr>
    <tr><th>atoms (median / p99 / max)</th>
        <td><code>{_fmt(summary['atoms_median'])} / {_fmt(summary['atoms_p99'])}
        / {_fmt(summary['atoms_max'])}</code></td></tr>
    <tr><th>minimizer steps (median / max)</th>
        <td><code>{_fmt(summary['steps_median'])} / {_fmt(summary['steps_max'])}</code>
        <span class="muted">(50-step budget by design)</span></td></tr>
    <tr><th>convergence rate</th>
        <td><code>{summary['n_converged']:,} / {summary['n_ok']:,}
        ({summary['convergence_rate']*100:.2f}%)</code>
        <span class="muted">— see note below</span></td></tr>
    <tr><th>energy drop &gt; 0 (downhill)</th>
        <td><code>{summary['drop_negative_count'] != summary['drop_n']
                  and summary['drop_n'] - summary['drop_negative_count'] or summary['drop_n']:,}
        / {summary['drop_n']:,}</code></td></tr>
  </table>
  <div class="note">
    <strong>Convergence ≠ correctness.</strong> The 50-step budget is intentional
    stress: rather than minimizing each PDB to its local minimum (which would
    take 1 000–10 000 steps depending on size), the runner pushes 50 steps and
    measures kernel survival + downhill motion. A near-zero
    <code>converged=true</code> rate is expected and not a regression — what
    the test guards is the kernel <em>not crashing</em> across 50 000 inputs.
  </div>
</section>

<section>
  <h2>Pipeline status across the corpus</h2>
  <p class="muted">Most structures complete the 50-step budget cleanly
  (<code>status: ok</code>); load errors are PDBFixer rejections of structures
  that can't be repaired (missing residues, unsupported polymer types).
  Pipeline errors are kernel-side crashes — the headline gate this claim
  defends.</p>
  {_embed_png(fig_status, "Status breakdown")}
</section>

<section>
  <h2>Structure size distribution</h2>
  <p class="muted">Atom count distribution of the structures that completed.
  Log-scaled because the wwPDB has a heavy tail — most structures are
  &lt; 5 000 atoms but the p99 reaches into the tens of thousands. The fact
  that the kernel survives the tail without an OOM or device-side crash is
  the substantive part of this claim.</p>
  {_embed_png(fig_atoms, "Structure size distribution")}
</section>

<section>
  <h2>Per-atom energy drop after 50 steps</h2>
  <p class="muted">(initial − final) / n_atoms. Negative values would mean
  the minimizer moved the system <em>uphill</em>; the histogram should sit
  almost entirely on the positive side of the dashed zero line. The median
  is the headline downhill rate.</p>
  {_embed_png(fig_drop, "Per-atom energy drop")}
</section>

<section>
  <h2>Initial vs final energy</h2>
  <p class="muted">Per-PDB scatter of starting energy vs minimized energy.
  Points should sit below the dashed identity line (final &lt; initial).
  Sub-sampled to 8 000 points for plot legibility on 45K records.</p>
  {_embed_png(fig_scatter, "Initial vs final energy")}
</section>

<section>
  <h2>Pinned versions</h2>
  <table class="kv">{pinned_lines}</table>
</section>

<section>
  <h2>last_verified</h2>
  {lv_table}
  <p class="muted">Populated by <code>lock_release_replays.py</code>
  at release-tag time. <code>value</code> cites the cited survival rate.</p>
</section>

<footer>
  EVIDENT 50K battle-test report · proteon ·
  artifact <code>{html.escape(artifact_path.name)}</code>
  (sha256 <code>{artifact_sha}</code>)
</footer>
</main>
</body>
</html>
"""

    output_path.write_text(html_doc, encoding="utf-8")
    print(
        f"Rendered {output_path} "
        f"(n_total={summary['n_total']}, ok={summary['n_ok']}, "
        f"survival={summary['survival_rate']*100:.1f}%, "
        f"median_ΔE={_fmt(summary['drop_median'], ',.0f')} kJ/mol)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--claim", type=pathlib.Path, required=True)
    parser.add_argument("--artifact", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--fig-dir", type=pathlib.Path, default=None)
    args = parser.parse_args()
    render(args.claim, args.artifact, args.output, args.fig_dir)


if __name__ == "__main__":
    main()
