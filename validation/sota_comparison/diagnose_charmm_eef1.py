#!/usr/bin/env python3
"""Diagnose the ferritin CHARMM19+EEF1 wrong-sign-total finding.

Tier-2 weak oracle (commit 8f979f4) flagged that ferritin's
charmm19_eef1 totals are POSITIVE on every v1 PDB (1crn +1849, 1ubq
+5591, ...) while OpenMM CHARMM36+OBC2 totals are NEGATIVE. The bug
appears concentrated in the EEF1 solvation term — fer_solv is positive
on 5/6 v1 PDBs and exactly zero on 1ake.

This script does a stepwise teardown on 1crn to localize the problem:

  Step 1: Load 1crn raw, compute CHARMM19+EEF1 energy WITHOUT
          minimization. Rules out "minimizer is driving us into the
          wrong basin" artifacts.
  Step 2: Print every component, including solvation, with sign.
  Step 3: Print the canonical expectation from the EEF1 parameter
          table — for 1crn (46 residues, all-protein) the
          contribution from peptide N (NH1) atoms ALONE should
          dominate at ≥-1100 kJ/mol from solvation.
  Step 4: Repeat with explicit H placement (no minimization).
  Step 5: Repeat with H + minimization for completeness.
  Step 6: Run the SAME structure through OpenMM CHARMM36+OBC2 for
          direct comparison.

If solvation is positive in step 1 (raw, no H, no minimize), the bug
is purely in the energy kernel — atom typing, parameter lookup, or
sign in the eef1_energy() accumulation. If solvation is negative
in step 1 but positive after H placement, the H placement is
corrupting the atom types. If solvation is negative in step 4 and
positive in step 5, the minimizer's gradient has wrong sign on the
EEF1 contribution.

Output: a single markdown table to stdout, no fancy formatting.

Usage:
    python diagnose_charmm_eef1.py [pdb_path]
"""

from __future__ import annotations

import argparse
import sys
import os

import ferritin

# Default test structure: 1crn from the v1 SOTA reference set.
DEFAULT_PDB = "/globalscratch/dateschn/ferritin-benchmark/sota_pdbs/1crn.pdb"

# Component keys for CHARMM19+EEF1
COMPONENTS = (
    "bond_stretch",
    "angle_bend",
    "torsion",
    "improper_torsion",
    "vdw",
    "electrostatic",
    "solvation",
)


def fmt(v):
    if v is None:
        return "None"
    if isinstance(v, float):
        return f"{v:>+12.3f}"
    return f"{v:>12}"


def print_energy(label: str, e: dict):
    print(f"\n  {label}")
    print(f"    {'component':<20} {'value (kJ/mol)':>14}")
    print(f"    {'-' * 20} {'-' * 14}")
    for k in COMPONENTS:
        v = e.get(k)
        marker = ""
        if k == "solvation" and v is not None:
            marker = "  ← wrong sign?" if v > 0 else "  (negative ✓)"
        print(f"    {k:<20} {fmt(v)}{marker}")
    print(f"    {'total':<20} {fmt(e.get('total'))}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", nargs="?", default=DEFAULT_PDB)
    args = parser.parse_args()

    if not os.path.isfile(args.pdb):
        print(f"ERROR: pdb not found: {args.pdb}", file=sys.stderr)
        return 1

    print("=" * 70)
    print(f"CHARMM19+EEF1 sign diagnostic on {os.path.basename(args.pdb)}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: raw structure, no H, no minimize, compute_energy
    # ------------------------------------------------------------------
    s_raw = ferritin.load(args.pdb)
    print(f"\nLoaded {args.pdb}: {s_raw.atom_count} atoms")
    try:
        n_residues = sum(1 for _ in s_raw.residues)
    except Exception:
        n_residues = "?"
    print(f"  residues: {n_residues}")

    e_raw = ferritin.compute_energy(s_raw, ff="charmm19_eef1", units="kJ/mol")
    print_energy("STEP 1: RAW (no H, no minimization)", e_raw)
    raw_solv_sign = "POSITIVE (BUG)" if e_raw["solvation"] > 0 else "negative (canonical)"
    print(f"\n    => solvation sign on raw structure: {raw_solv_sign}")

    # ------------------------------------------------------------------
    # Canonical expectation (back-of-envelope from EEF1 parameter table)
    # ------------------------------------------------------------------
    print("\n  Canonical expectation for a typical small protein:")
    print(f"    n_residues × NH1 dg_ref ≈ {n_residues} × -5.95 kcal/mol "
          f"= {n_residues * -5.95:.1f} kcal/mol "
          f"= {n_residues * -5.95 * 4.184:.1f} kJ/mol")
    print("    (just from peptide nitrogens — actual sum should be more negative)")

    # ------------------------------------------------------------------
    # Step 2: H placed but no minimization
    # ------------------------------------------------------------------
    s_h = ferritin.load(args.pdb)
    reports_h = ferritin.batch_prepare(
        [s_h],
        reconstruct=False,
        hydrogens="all",
        minimize=False,
        strip_hydrogens=False,
        ff="charmm19_eef1",
    )
    print(f"\n  After H placement: {s_h.atom_count} atoms (was {s_raw.atom_count})")
    e_h = ferritin.compute_energy(s_h, ff="charmm19_eef1", units="kJ/mol")
    print_energy("STEP 2: H placed, NO minimization", e_h)

    # ------------------------------------------------------------------
    # Step 3: H + 200 LBFGS steps
    # ------------------------------------------------------------------
    s_min = ferritin.load(args.pdb)
    reports_min = ferritin.batch_prepare(
        [s_min],
        reconstruct=False,
        hydrogens="all",
        minimize=True,
        minimize_method="lbfgs",
        minimize_steps=200,
        gradient_tolerance=0.1,
        strip_hydrogens=False,
        ff="charmm19_eef1",
    )
    r = reports_min[0]
    print(f"\nSTEP 3: H + LBFGS minimization ({r.minimizer_steps} steps, "
          f"converged={r.converged})")
    print(f"  initial_energy: {r.initial_energy:>+14.3f} kJ/mol")
    print(f"  final_energy:   {r.final_energy:>+14.3f} kJ/mol")
    e_min = ferritin.compute_energy(s_min, ff="charmm19_eef1", units="kJ/mol")
    print_energy("  re-evaluated post-minimization", e_min)

    # ------------------------------------------------------------------
    # Step 4: comparison with raw AMBER96 on the same input (sanity)
    # ------------------------------------------------------------------
    s_amber = ferritin.load(args.pdb)
    e_amber = ferritin.compute_energy(s_amber, ff="amber96", units="kJ/mol")
    print_energy("STEP 4: AMBER96 on the same raw structure (control)", e_amber)
    print("\n  (AMBER96 has no solvation term, so 'solvation' should be 0 or None.)")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Raw CHARMM19+EEF1 solvation: {e_raw['solvation']:>+12.3f} kJ/mol")
    print(f"  H+placement solvation:       {e_h['solvation']:>+12.3f} kJ/mol")
    print(f"  Post-minimization solvation: {e_min['solvation']:>+12.3f} kJ/mol")
    print()
    print(f"  Raw CHARMM19+EEF1 total:     {e_raw['total']:>+12.3f} kJ/mol")
    print(f"  H+placement total:           {e_h['total']:>+12.3f} kJ/mol")
    print(f"  Post-minimization total:     {e_min['total']:>+12.3f} kJ/mol")
    print()
    print(f"  Raw AMBER96 total (control): {e_amber['total']:>+12.3f} kJ/mol")
    print()
    print("Localization key:")
    print("  - solvation > 0 in STEP 1 (raw):")
    print("      ⇒ bug is in eef1_energy() kernel itself")
    print("        (atom typing, dg_ref accumulation, or sign in pair correction)")
    print("  - solvation < 0 in STEP 1, > 0 in STEP 2 (H placed):")
    print("      ⇒ H placement corrupts atom types (e.g. NH1 reclassified as NH2)")
    print("  - solvation < 0 in STEP 2, > 0 in STEP 3 (minimized):")
    print("      ⇒ minimizer gradient on EEF1 has wrong sign — energy and force")
    print("        diverge, structure converges to local max instead of min")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
