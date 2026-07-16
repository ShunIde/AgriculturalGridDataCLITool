#!/usr/bin/env python3
"""
breeding.py — climate-informed variety recommendation, phase 1.

This is a *from-scratch* reimplementation of the core idea behind PyBrOpS
(Shrote & Beavis, G3, 2024: https://doi.org/10.1093/g3journal/jkae199) —
finding Pareto-optimal trade-offs between competing breeding objectives via
a multi-objective evolutionary algorithm — built directly on `pymoo`
(the same optimization library PyBrOpS itself wraps). It does NOT depend on
the `pybrops` package.

Why not just use `pybrops`: `pybrops` hard-requires `cyvcf2` for genotype
loading, and `cyvcf2` ships no Windows wheel (Linux/macOS only, and building
it from source needs the `htslib` C library). A genuine `pybrops`
integration — real VCF loading via cyvcf2, the full `OptimalContributionSubsetSelection`
protocol, multi-generation breeding simulation — is planned for a later
phase, to be run on Linux.

What this module does instead, staying as close to the paper's actual demo
as possible:
  1. Parses real genotype data (942 lines x 2000 SNPs, the Wisconsin
     Diversity maize panel bundled in PyBrOpS's own examples) using only
     the stdlib `gzip` module — no cyvcf2 needed.
  2. Synthesizes two negatively-correlated trait effects ("cold_tolerance",
     "drought_tolerance") over those markers, mirroring PyBrOpS's own
     tri-objective OCS example almost exactly (same multivariate-normal
     marker-effect trick). These are simulated traits, not real phenotypes
     — the dataset has no public cold/drought phenotyping.
  3. Computes a genomic relationship (kinship) matrix from realized allele
     dosage (VanRaden method 1) as a stand-in for PyBrOpS's molecular
     coancestry matrix.
  4. Runs a genuine multi-objective subset-selection GA (NSGA-II via
     `pymoo`) to find the Pareto frontier trading off mean breeding value
     for each trait against within-subset inbreeding — the same
     conceptual problem as Optimal Contribution Selection.
  5. Given a field's irrigation/frost thresholds as a proxy for how much
     that field leans on a crop's innate cold/drought tolerance, picks the
     frontier point that best matches those priorities and reports which
     real panel lines it corresponds to.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

log = logging.getLogger("farm_report.breeding")

TRAIT_NAMES = ("cold_tolerance", "drought_tolerance")

DEFAULT_VCF_PATH = Path(__file__).parent / "data" / "widiv_2000SNPs.vcf.gz"
DEFAULT_CACHE_PATH = Path(__file__).parent / "data" / "ocs_frontier_cache.npz"

# Synthetic trait model parameters — chosen to mirror PyBrOpS's own
# triobjective_OCS_pareto_frontier.py example (same seed value, same
# negatively-correlated effect covariance, same trait-intercept magnitudes).
TRAIT_SEED = 31621463
TRAIT_COV = [[1.0, -0.4], [-0.4, 1.0]]
TRAIT_INTERCEPTS = (10.0, 25.0)


# ==========================================================================
# Genotype loading (pure stdlib — no cyvcf2)
# ==========================================================================


@dataclass
class GenotypePanel:
    """Real diploid genotypes as additive allele dosage (0, 1, or 2 copies of ALT)."""

    taxa: list[str]
    dosage: np.ndarray  # shape (n_taxa, n_variants), dtype int8


def load_vcf_dosage(path: Path) -> GenotypePanel:
    """
    Parses a (possibly gzipped) biallelic VCF into an additive dosage matrix.
    Only reads the GT field; missing calls ('.') count as 0 ALT copies.
    """
    log.info("Parsing genotype panel from %s (pure-Python VCF reader)...", path)
    taxa: list[str] | None = None
    rows: list[list[int]] = []

    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                taxa = line.rstrip("\n").split("\t")[9:]
                continue
            if taxa is None:
                raise ValueError(f"{path}: no #CHROM header found before data rows")

            cols = line.rstrip("\n").split("\t")
            fmt = cols[8].split(":")
            gt_idx = fmt.index("GT")

            dosages = []
            for sample_field in cols[9:]:
                gt = sample_field.split(":")[gt_idx]
                alleles = gt.replace("|", "/").split("/")
                dosages.append(sum(1 for a in alleles if a == "1"))
            rows.append(dosages)

    if taxa is None or not rows:
        raise ValueError(f"{path}: no genotype data found")

    dosage = np.array(rows, dtype=np.int8).T  # (n_variants, n_taxa) -> (n_taxa, n_variants)
    log.info("Loaded %d taxa x %d variants.", dosage.shape[0], dosage.shape[1])
    return GenotypePanel(taxa=taxa, dosage=dosage)


# ==========================================================================
# Synthetic bi-trait genomic model
# ==========================================================================


def build_synthetic_marker_effects(n_variants: int, seed: int = TRAIT_SEED) -> np.ndarray:
    """
    Returns a (n_variants, 2) array of marker effects for two competing,
    negatively-correlated synthetic traits: cold_tolerance, drought_tolerance.
    """
    rng = np.random.default_rng(seed)
    return rng.multivariate_normal(mean=[0.0, 0.0], cov=TRAIT_COV, size=n_variants)


def compute_gebv(panel: GenotypePanel, marker_effects: np.ndarray) -> np.ndarray:
    """Additive genomic estimated breeding values: intercept + dosage @ effects."""
    return np.array(TRAIT_INTERCEPTS) + panel.dosage.astype(np.float64) @ marker_effects


def compute_grm(dosage: np.ndarray) -> np.ndarray:
    """
    VanRaden method-1 genomic relationship (kinship) matrix from allele dosage.
    Standin for PyBrOpS's DenseMolecularCoancestryMatrixFactory.
    """
    dosage = dosage.astype(np.float64)
    allele_freq = dosage.mean(axis=0) / 2.0
    centered = dosage - 2.0 * allele_freq
    denom = 2.0 * np.sum(allele_freq * (1.0 - allele_freq))
    if denom <= 0:
        raise ValueError("Genotype panel has no polymorphic markers; cannot compute GRM")
    return (centered @ centered.T) / denom


# ==========================================================================
# Multi-objective subset selection (NSGA-II via pymoo)
# ==========================================================================


def _import_pymoo():
    try:
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.core.crossover import Crossover
        from pymoo.core.mutation import Mutation
        from pymoo.core.problem import Problem
        from pymoo.core.sampling import Sampling
        from pymoo.optimize import minimize
    except ImportError as e:
        raise ImportError(
            "The 'variety' subcommand requires pymoo. Install it with:\n"
            "    pip install -r requirements-breeding.txt"
        ) from e
    return NSGA2, Crossover, Mutation, Problem, Sampling, minimize


def _build_subset_ga_classes():
    """
    Builds the fixed-size subset-selection Problem/Sampling/Crossover/Mutation
    classes. Deferred inside a function so importing breeding.py doesn't
    require pymoo unless the variety subcommand is actually invoked.
    """
    NSGA2, Crossover, Mutation, Problem, Sampling, minimize = _import_pymoo()

    class SubsetOCSProblem(Problem):
        """
        Decision variable: a boolean mask of length n_taxa selecting exactly
        `k` lines. Objectives (all minimized):
          f1 = -mean cold_tolerance GEBV of selected subset
          f2 = -mean drought_tolerance GEBV of selected subset
          f3 =  mean pairwise kinship within the selected subset (inbreeding cost)
        This is the same conceptual trade-off as PyBrOpS's Optimal
        Contribution Selection: genetic gain per trait vs. diversity.
        """

        def __init__(self, gebv: np.ndarray, grm: np.ndarray, k: int):
            n = gebv.shape[0]
            super().__init__(n_var=n, n_obj=3, n_ieq_constr=0, xl=0, xu=1, vtype=bool)
            self.gebv = gebv
            self.grm = grm
            self.k = k

        def _evaluate(self, X, out, *args, **kwargs):
            n_pop = X.shape[0]
            f = np.empty((n_pop, 3), dtype=np.float64)
            triu = np.triu_indices(self.k, k=1)
            for i in range(n_pop):
                idx = np.flatnonzero(X[i])
                f[i, 0] = -self.gebv[idx, 0].mean()
                f[i, 1] = -self.gebv[idx, 1].mean()
                sub = self.grm[np.ix_(idx, idx)]
                f[i, 2] = sub[triu].mean()
            out["F"] = f

    class SubsetSampling(Sampling):
        def _do(self, problem, n_samples, **kwargs):
            X = np.zeros((n_samples, problem.n_var), dtype=bool)
            for i in range(n_samples):
                idx = np.random.choice(problem.n_var, problem.k, replace=False)
                X[i, idx] = True
            return X

    class SubsetCrossover(Crossover):
        """Child = a random k-sized pick from the union of both parents' selections."""

        def __init__(self):
            super().__init__(2, 1)

        def _do(self, problem, X, **kwargs):
            _, n_matings, n_var = X.shape
            Y = np.zeros((1, n_matings, n_var), dtype=bool)
            for m in range(n_matings):
                union = np.flatnonzero(X[0, m] | X[1, m])
                k = problem.k
                if len(union) >= k:
                    chosen = np.random.choice(union, k, replace=False)
                else:
                    remaining = np.setdiff1d(np.arange(n_var), union)
                    extra = np.random.choice(remaining, k - len(union), replace=False)
                    chosen = np.concatenate([union, extra])
                Y[0, m, chosen] = True
            return Y

    class SubsetMutation(Mutation):
        """With probability `prob`, swap `nswap` selected lines for unselected ones."""

        def __init__(self, prob: float = 0.3, nswap: int = 2):
            super().__init__()
            self.prob = prob
            self.nswap = nswap

        def _do(self, problem, X, **kwargs):
            Xp = X.copy()
            for i in range(Xp.shape[0]):
                if np.random.random() >= self.prob:
                    continue
                sel = np.flatnonzero(Xp[i])
                unsel = np.flatnonzero(~Xp[i])
                nswap = min(self.nswap, len(sel), len(unsel))
                if nswap == 0:
                    continue
                out_idx = np.random.choice(sel, nswap, replace=False)
                in_idx = np.random.choice(unsel, nswap, replace=False)
                Xp[i, out_idx] = False
                Xp[i, in_idx] = True
            return Xp

    return SubsetOCSProblem, SubsetSampling, SubsetCrossover, SubsetMutation, NSGA2, minimize


@dataclass
class ParetoFrontier:
    taxa: list[str]
    masks: np.ndarray  # (n_solutions, n_taxa) bool
    objectives: np.ndarray  # (n_solutions, 3): [-meanBV_cold, -meanBV_drought, mean_kinship]
    k: int


def _params_fingerprint(vcf_path: Path, k: int, pop_size: int, n_gen: int, seed: int) -> str:
    vcf_stat = vcf_path.stat()
    key = f"{vcf_path.name}:{vcf_stat.st_size}:{k}:{pop_size}:{n_gen}:{seed}:{TRAIT_SEED}"
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def compute_pareto_frontier(
    panel: GenotypePanel,
    gebv: np.ndarray,
    grm: np.ndarray,
    k: int = 20,
    pop_size: int = 80,
    n_gen: int = 120,
    seed: int = 1,
) -> ParetoFrontier:
    (
        SubsetOCSProblem,
        SubsetSampling,
        SubsetCrossover,
        SubsetMutation,
        NSGA2,
        minimize,
    ) = _build_subset_ga_classes()

    log.info(
        "Running NSGA-II subset selection (k=%d parents out of %d lines, "
        "pop_size=%d, n_gen=%d)... this may take a little while.",
        k, len(panel.taxa), pop_size, n_gen,
    )
    problem = SubsetOCSProblem(gebv, grm, k)
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=SubsetSampling(),
        crossover=SubsetCrossover(),
        mutation=SubsetMutation(),
        eliminate_duplicates=False,
    )
    res = minimize(problem, algorithm, ("n_gen", n_gen), seed=seed, verbose=False)
    log.info("Pareto frontier computed: %d non-dominated solutions.", len(res.F))
    return ParetoFrontier(taxa=panel.taxa, masks=res.X.astype(bool), objectives=res.F, k=k)


def compute_or_load_frontier(
    panel: GenotypePanel,
    gebv: np.ndarray,
    grm: np.ndarray,
    vcf_path: Path,
    cache_path: Path = DEFAULT_CACHE_PATH,
    k: int = 20,
    pop_size: int = 80,
    n_gen: int = 120,
    seed: int = 1,
    force_regen: bool = False,
) -> ParetoFrontier:
    """
    NSGA-II over ~1000s of evaluations is the expensive part of this pipeline
    (seconds-to-minutes, not instant) — cache it to disk so repeated CLI runs
    (e.g. one per field) don't recompute it. Cache is invalidated whenever the
    fingerprint of (vcf file, k, pop_size, n_gen, seed, trait seed) changes.
    """
    fingerprint = _params_fingerprint(vcf_path, k, pop_size, n_gen, seed)

    if not force_regen and cache_path.exists():
        with np.load(cache_path, allow_pickle=True) as data:
            if str(data["fingerprint"]) == fingerprint:
                log.info("Loaded cached Pareto frontier from %s", cache_path)
                return ParetoFrontier(
                    taxa=list(data["taxa"]),
                    masks=data["masks"],
                    objectives=data["objectives"],
                    k=int(data["k"]),
                )
            log.info("Cache at %s is stale (params changed); recomputing.", cache_path)

    frontier = compute_pareto_frontier(panel, gebv, grm, k=k, pop_size=pop_size, n_gen=n_gen, seed=seed)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        cache_path,
        fingerprint=fingerprint,
        taxa=np.array(frontier.taxa, dtype=object),
        masks=frontier.masks,
        objectives=frontier.objectives,
        k=frontier.k,
    )
    log.info("Cached Pareto frontier to %s", cache_path)
    return frontier


# ==========================================================================
# Field-driven recommendation
# ==========================================================================


# Kept deliberately small: min-max normalization stretches whichever
# frontier point has the best kinship all the way to 0 penalty and the
# worst all the way to 1, regardless of how large the real gap is. A
# coefficient close to 1 would let that structural swing override a field's
# actual trait preference; keeping it well below the smallest expected
# per-trait weight means it mainly breaks near-ties between otherwise
# similarly-preferred points instead of overriding preference outright.
DIVERSITY_PENALTY_WEIGHT = 0.15


@dataclass
class VarietyRecommendation:
    field_id: str
    field_name: str
    weights: tuple[float, float]  # (cold_weight, drought_weight)
    selected_taxa: list[str]
    mean_bv_cold: float
    mean_bv_drought: float
    mean_kinship: float
    top_lines: list[tuple[str, float, float]]  # (taxon, bv_cold, bv_drought), best few


def field_trait_weights(frost_temp_c: float, soil_moisture_min: float) -> tuple[float, float]:
    """
    Maps a field's configured thresholds to a (cold_weight, drought_weight)
    preference vector. A *lower* frost_temp_c threshold means the grower
    only intervenes at more extreme cold, i.e. is leaning more on the crop's
    innate cold tolerance to survive milder cold unassisted — so it gets
    more weight on cold_tolerance. Same logic for soil_moisture_min and
    drought_tolerance. This is a phase-1 heuristic based on configured
    thresholds only (not live/historical weather) — see README.
    """
    cold_raw = -frost_temp_c
    drought_raw = -soil_moisture_min
    # shift both onto a common [0, positive] scale before normalizing so a
    # single very negative raw value can't zero out the other weight
    cold_shifted = cold_raw + 10.0
    drought_shifted = drought_raw + 10.0
    total = cold_shifted + drought_shifted
    if total <= 0:
        return (0.5, 0.5)
    return (cold_shifted / total, drought_shifted / total)


def recommend_for_field(
    field: Any,  # farm_report.Field — kept as Any to avoid a circular import
    frontier: ParetoFrontier,
    top_n: int = 5,
) -> VarietyRecommendation:
    cold_weight, drought_weight = field_trait_weights(
        field.thresholds["frost_temp_c"], field.thresholds["soil_moisture_min"]
    )

    gain_cold = -frontier.objectives[:, 0]
    gain_drought = -frontier.objectives[:, 1]
    kinship = frontier.objectives[:, 2]

    def normalize(v: np.ndarray) -> np.ndarray:
        span = v.max() - v.min()
        return (v - v.min()) / span if span > 0 else np.zeros_like(v)

    score = (
        cold_weight * normalize(gain_cold)
        + drought_weight * normalize(gain_drought)
        - DIVERSITY_PENALTY_WEIGHT * normalize(kinship)
    )
    best_idx = int(np.argmax(score))

    mask = frontier.masks[best_idx]
    selected_idx = np.flatnonzero(mask)
    selected_taxa = [frontier.taxa[i] for i in selected_idx]

    return VarietyRecommendation(
        field_id=field.id,
        field_name=field.name,
        weights=(cold_weight, drought_weight),
        selected_taxa=selected_taxa,
        mean_bv_cold=float(gain_cold[best_idx]),
        mean_bv_drought=float(gain_drought[best_idx]),
        mean_kinship=float(kinship[best_idx]),
        top_lines=[],  # filled in by build_recommendations, which has GEBVs in scope
    )


def build_recommendations(
    fields: list[Any],
    frontier: ParetoFrontier,
    gebv: np.ndarray,
    taxa: list[str],
    top_n: int = 5,
) -> list[VarietyRecommendation]:
    taxon_to_row = {t: i for i, t in enumerate(taxa)}
    recs = []
    for field in fields:
        rec = recommend_for_field(field, frontier, top_n=top_n)
        scored = sorted(
            (
                (t, float(gebv[taxon_to_row[t], 0]), float(gebv[taxon_to_row[t], 1]))
                for t in rec.selected_taxa
            ),
            key=lambda row: rec.weights[0] * row[1] + rec.weights[1] * row[2],
            reverse=True,
        )
        rec.top_lines = scored[:top_n]
        recs.append(rec)
    return recs


# ==========================================================================
# Report generation
# ==========================================================================


class VarietyReportGenerator:
    def __init__(self, farm_name: str, recommendations: list[VarietyRecommendation]):
        self.farm_name = farm_name
        self.recommendations = recommendations

    def build(self) -> str:
        lines: list[str] = []
        lines.append(f"# Variety Recommendations — {self.farm_name}")
        lines.append(
            "_Climate-informed genomic selection over the Wisconsin Diversity (WiDiv) "
            "maize panel — see caveats below._\n"
        )
        lines.append(
            "> **Simulated data notice:** `cold_tolerance` and `drought_tolerance` are "
            "*synthetic* trait effects generated over real SNP genotypes (see "
            "`breeding.py`), not measured phenotypes — the WiDiv panel has no public "
            "cold/drought phenotyping. This demonstrates the selection methodology, not "
            "an agronomic recommendation.\n"
        )

        if not self.recommendations:
            lines.append("_No corn fields found in this farm layout — nothing to recommend._")
            return "\n".join(lines)

        for rec in self.recommendations:
            lines.append(f"## {rec.field_name} (`{rec.field_id}`)")
            lines.append(
                f"- Trait weights derived from field thresholds: "
                f"cold_tolerance {rec.weights[0]:.2f}, drought_tolerance {rec.weights[1]:.2f}"
            )
            lines.append(
                f"- Selected cross candidates: {len(rec.selected_taxa)} lines "
                f"(mean cold BV {rec.mean_bv_cold:.2f}, mean drought BV {rec.mean_bv_drought:.2f}, "
                f"mean kinship {rec.mean_kinship:.4f})"
            )
            lines.append("")
            lines.append("| Line | Cold-tolerance BV | Drought-tolerance BV |")
            lines.append("|---|---|---|")
            for taxon, bv_cold, bv_drought in rec.top_lines:
                lines.append(f"| {taxon} | {bv_cold:.2f} | {bv_drought:.2f} |")
            lines.append("")

        return "\n".join(lines)
