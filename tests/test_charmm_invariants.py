"""Internal invariant tests for CHARMM19+EEF1.

Why this file exists
--------------------
CHARMM19+EEF1 has been the production default for `batch_prepare` since
commit 73248f7, but has NO external oracle in the repo. The AMBER96
oracle (tests/oracle/test_ball_energy.py) validates AMBER96 against the
BALL Julia implementation to 0.02% on heavy-only 1crn — but AMBER96 is
not what users actually run. The well-validated FF is the one nobody
uses; the FF that runs on every batch_prepare call is checked only by
"the number is not zero".

This file fills the gap by enforcing INVARIANTS — properties that must
hold regardless of force-field parameter values. Any violation indicates
a real bug. These tests catch:

  - NaN / Inf in any energy component (numerical blowup)
  - Missing components (CHARMM19-EEF1 should produce all seven)
  - Σ(components) ≠ total (component-accounting bug)
  - Negative bond/angle energy (sign error in a harmonic potential)
  - Zero solvation (EEF1 silently fell back to vacuum CHARMM)
  - Non-determinism (same input → different output)
  - Energy increasing under minimization (minimizer broken)
  - CHARMM and AMBER giving the same answer (one is silently falling
    back to the other)

What these tests do NOT do
--------------------------
They do NOT validate parameter correctness — wrong K_b values for a
specific bond type, wrong charges on a specific atom, etc. — because
that requires an external oracle. The "external CHARMM19-EEF1 oracle"
is a known gap (Tier-3 work in the SOTA validation roadmap). These
internal invariants are a strict subset of what an oracle would catch.
"""

import math
import os

import pytest

import ferritin

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_PDBS_DIR = os.path.join(REPO, "test-pdbs")

# v1 SOTA reference set restricted to the PDBs available in test-pdbs/.
# All four are small, clean, single-chain proteins with no nucleic acids
# and no exotic ligands. They cover residue counts from 46 (1crn) to
# 214 (1ake) — small enough that all 36 invariant tests run in a few
# seconds, large enough that "energy is small" doesn't trivially hold.
V1_PDBS = [
    ("1crn", os.path.join(TEST_PDBS_DIR, "1crn.pdb")),
    ("1ubq", os.path.join(TEST_PDBS_DIR, "1ubq.pdb")),
    ("1bpi", os.path.join(TEST_PDBS_DIR, "1bpi.pdb")),
    ("1ake", os.path.join(TEST_PDBS_DIR, "1ake.pdb")),
]

# Force-field name used by ferritin's `compute_energy(ff=...)` API.
CHARMM = "charmm19_eef1"

# All seven CHARMM19-EEF1 component keys. `total` is checked separately.
# Note: CHARMM does not split improper from proper torsion the way some
# AMBER implementations do — but ferritin still returns improper_torsion
# as a separate key (may be 0). Keep it in the list so the
# "components-present" test catches a future schema regression.
CHARMM_COMPONENTS = (
    "bond_stretch",
    "angle_bend",
    "torsion",
    "improper_torsion",
    "vdw",
    "electrostatic",
    "solvation",
)


@pytest.fixture(params=V1_PDBS, ids=[p[0] for p in V1_PDBS])
def charmm_energy(request):
    """(name, energy_dict) tuple, parametrized over the v1 reference set.

    Uses `compute_energy` (not `batch_prepare`) so the input geometry is
    the raw PDB — no H placement, no minimization. This is the cleanest
    "is the energy kernel correct on a known input" measurement.
    """
    name, path = request.param
    s = ferritin.load(path)
    e = ferritin.compute_energy(s, ff=CHARMM, units="kJ/mol")
    return name, e


# =========================================================================
# Numerical sanity
# =========================================================================


class TestNoNanInf:
    """Every energy component is finite (no NaN, no Inf)."""

    def test_total_is_finite(self, charmm_energy):
        name, e = charmm_energy
        assert math.isfinite(e["total"]), f"{name}: total is {e['total']}"

    def test_components_are_finite(self, charmm_energy):
        name, e = charmm_energy
        for comp in CHARMM_COMPONENTS:
            v = e.get(comp)
            assert v is not None, f"{name}: component {comp} is missing"
            assert math.isfinite(v), f"{name}: {comp} is {v}"


class TestComponentsPresent:
    """All seven CHARMM components plus `total` are present in the dict."""

    def test_all_keys_present(self, charmm_energy):
        name, e = charmm_energy
        for comp in CHARMM_COMPONENTS:
            assert comp in e, f"{name}: missing key {comp!r}"
        assert "total" in e, f"{name}: missing key 'total'"


class TestSumsMatchTotal:
    """Σ(components) ≈ total — catches accounting bugs like dropping a
    component from the sum.

    Tolerance: max(1e-3 * |total|, 1e-3) kJ/mol. Relative for large
    energies (a 10000 kJ/mol total can drift ~10 kJ/mol from float
    accumulation), absolute floor for small energies.
    """

    def test_components_sum_to_total(self, charmm_energy):
        name, e = charmm_energy
        component_sum = sum(e[c] for c in CHARMM_COMPONENTS)
        diff = abs(component_sum - e["total"])
        tol = max(1e-3 * abs(e["total"]), 1e-3)
        assert diff < tol, (
            f"{name}: Σ(components) {component_sum:.6f} != total "
            f"{e['total']:.6f} (diff {diff:.4e}, tol {tol:.4e})"
        )


# =========================================================================
# Physical sanity
# =========================================================================


class TestSigns:
    """Harmonic potentials are sums of squares — they MUST be ≥ 0.
    A negative value would indicate a sign error in the energy kernel.
    """

    def test_bond_stretch_nonnegative(self, charmm_energy):
        name, e = charmm_energy
        assert e["bond_stretch"] >= 0, (
            f"{name}: bond_stretch {e['bond_stretch']:.6f} is negative — "
            "harmonic potential should always be ≥ 0"
        )

    def test_angle_bend_nonnegative(self, charmm_energy):
        name, e = charmm_energy
        assert e["angle_bend"] >= 0, (
            f"{name}: angle_bend {e['angle_bend']:.6f} is negative — "
            "harmonic potential should always be ≥ 0"
        )


class TestSolvationActive:
    """EEF1 implicit solvation must contribute non-trivially.

    If solvation == 0 across a real protein, EEF1 silently fell back to
    vacuum CHARMM and the production default is broken. The threshold
    of 1.0 kJ/mol is far below what any real protein produces — buried
    hydrophobics and exposed polars give tens to hundreds of kJ/mol —
    so it only catches the "EEF1 isn't running at all" failure mode,
    not "EEF1 has slightly wrong parameters".
    """

    def test_solvation_nontrivial(self, charmm_energy):
        name, e = charmm_energy
        assert abs(e["solvation"]) > 1.0, (
            f"{name}: solvation = {e['solvation']:.6f} kJ/mol — EEF1 "
            "should produce a non-trivial implicit-solvation contribution"
        )


# =========================================================================
# Determinism
# =========================================================================


class TestDeterminism:
    """Same input → identical output across two compute_energy calls.

    Catches non-deterministic accumulation order, uninitialized memory,
    and parallelism with race conditions in the energy kernel. The
    tolerance is exact equality (abs=0) because compute_energy on the
    same input should be bit-identical, not just numerically close.
    """

    @pytest.mark.parametrize("name,path", V1_PDBS, ids=[p[0] for p in V1_PDBS])
    def test_two_calls_identical(self, name, path):
        s = ferritin.load(path)
        e1 = ferritin.compute_energy(s, ff=CHARMM)
        e2 = ferritin.compute_energy(s, ff=CHARMM)
        for key in CHARMM_COMPONENTS + ("total",):
            assert e1[key] == e2[key], (
                f"{name}: {key} not deterministic — "
                f"call1={e1[key]!r}, call2={e2[key]!r}"
            )


# =========================================================================
# Minimization
# =========================================================================


class TestMinimizationDecreases:
    """LBFGS minimization must end at energy ≤ start.

    If final > initial, the minimizer is broken (wrong gradient sign,
    bad step acceptance, etc.). LBFGS line search can briefly land
    above the start during a search step, but the final accepted state
    should always be at or below the initial energy.
    """

    @pytest.mark.parametrize("name,path", V1_PDBS, ids=[p[0] for p in V1_PDBS])
    def test_final_not_above_initial(self, name, path):
        s = ferritin.load(path)
        reports = ferritin.batch_prepare(
            [s],
            reconstruct=False,
            hydrogens="all",
            minimize=True,
            minimize_method="lbfgs",
            # 200 steps is enough to reach a stable minimum on a small
            # protein and keeps the test suite fast.
            minimize_steps=200,
            gradient_tolerance=0.1,
            strip_hydrogens=False,
            ff=CHARMM,
        )
        r = reports[0]
        # Allow a tiny numerical overshoot: 0.1 kJ/mol absolute or
        # 0.01% relative, whichever is larger.
        tol = max(0.1, 1e-4 * abs(r.initial_energy))
        assert r.final_energy <= r.initial_energy + tol, (
            f"{name}: final_energy {r.final_energy:.4f} > "
            f"initial_energy {r.initial_energy:.4f} (tol {tol:.4f}) — "
            "LBFGS minimization is not decreasing the energy"
        )


# =========================================================================
# CHARMM is actually running (not silently falling back to AMBER)
# =========================================================================


class TestCharmmDistinctFromAmber:
    """If CHARMM and AMBER return the same number, one of them is
    silently falling back to the other. Sanity check that we're
    actually evaluating CHARMM19-EEF1 when we ask for it.
    """

    def test_charmm_total_differs_from_amber(self):
        s = ferritin.load(os.path.join(TEST_PDBS_DIR, "1crn.pdb"))
        e_charmm = ferritin.compute_energy(s, ff=CHARMM)
        e_amber = ferritin.compute_energy(s, ff="amber96")
        assert e_charmm["total"] != e_amber["total"], (
            "CHARMM and AMBER returned identical totals on 1crn — "
            "one force field is silently falling back to the other"
        )

    def test_charmm_has_solvation_amber_does_not(self):
        """The most distinctive signal: AMBER96 in vacuum has no
        solvation term, CHARMM19+EEF1 always does. If both have the
        same solvation value, EEF1 is not active.
        """
        s = ferritin.load(os.path.join(TEST_PDBS_DIR, "1crn.pdb"))
        e_charmm = ferritin.compute_energy(s, ff=CHARMM)
        e_amber = ferritin.compute_energy(s, ff="amber96")
        assert e_charmm["solvation"] != e_amber["solvation"], (
            f"CHARMM solvation ({e_charmm['solvation']}) == AMBER "
            f"solvation ({e_amber['solvation']}) — EEF1 is not running"
        )
