#!/usr/bin/env python3
"""
breeding_pybrops.py — climate-informed variety recommendation, phase 2.

The real `pybrops` integration promised in breeding.py's phase-1 docstring:
real VCF loading via `cyvcf2` (through pybrops's own `DensePhasedGenotypeMatrix
.from_vcf`), the actual `OptimalContributionSubsetSelection` protocol, and a
multi-generation breeding simulation using pybrops's `TwoWayCross` mating
protocol. Linux/macOS only (`cyvcf2` ships no Windows wheel) — see
breeding.py for the from-scratch Windows-compatible fallback.

Uses the same synthetic bi-trait genomic model as phase 1 (same trait seed,
covariance, intercepts — see breeding.TRAIT_SEED etc.) so the two engines are
directly comparable, and reuses breeding.py's ParetoFrontier / field-scoring /
report-generation code by producing frontiers in the same shape: objectives
columns (-meanBV_cold, -meanBV_drought, mean_kinship), masks as boolean
(n_solutions, n_taxa) arrays.

**Known approximation:** pybrops's mating protocols (`TwoWayCross` etc.)
require a genetic map (`vrnt_xoprob`, crossover probability per marker) to
simulate recombination. No real genetic map exists for this 2000-SNP WiDiv
panel, so `assign_approximate_genetic_map` fabricates one at a flat
1 cM/Mb genome-wide average (a common maize-wide rule of thumb, not this
panel's actual recombination landscape) via `StandardGeneticMap` +
`HaldaneMapFunction`. This only affects the multi-generation breeding
simulation (recombination between markers within a chromosome) — the
single-generation OCS frontier (variety recommendations) doesn't use it at
all, since selecting existing lines doesn't require simulating meiosis.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Any

import numpy as np

import breeding
from breeding import (
    DEFAULT_VCF_PATH,
    TRAIT_COV,
    TRAIT_INTERCEPTS,
    TRAIT_NAMES,
    TRAIT_SEED,
    ParetoFrontier,
)

log = logging.getLogger("farm_report.breeding_pybrops")

DEFAULT_CACHE_PATH = Path(__file__).parent / "data" / "ocs_frontier_cache_pybrops.npz"

# Flat genome-wide recombination rate used to fabricate a genetic map for
# this panel (see module docstring for why). 1 cM/Mb is a commonly-cited
# maize-wide average, not a per-chromosome measurement for these markers.
APPROX_CM_PER_MB = 1.0


def _import_pybrops():
    try:
        from pybrops.breed.prot.mate.TwoWayCross import TwoWayCross
        from pybrops.breed.prot.sel.OptimalContributionSelection import (
            OptimalContributionSubsetSelection,
        )
        from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import (
            DenseAdditiveLinearGenomicModel,
        )
        from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
        from pybrops.popgen.cmat.fcty.DenseVanRadenCoancestryMatrixFactory import (
            DenseVanRadenCoancestryMatrixFactory,
        )
        from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
        from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap
        from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
    except ImportError as e:
        raise ImportError(
            "The 'pybrops' engine requires the real pybrops + cyvcf2 packages "
            "(Linux/macOS only). Install them with:\n"
            "    pip install -r requirements-breeding-pybrops.txt\n"
            "Note: pybrops 1.0.3 is incompatible with numpy>=2.0 (uses a "
            "removed alias), so numpy must stay pinned below 2.0 - the "
            "requirements file already constrains this."
        ) from e
    return (
        DensePhasedGenotypeMatrix,
        DenseAdditiveLinearGenomicModel,
        DenseVanRadenCoancestryMatrixFactory,
        OptimalContributionSubsetSelection,
        NSGA2SubsetGeneticAlgorithm,
        TwoWayCross,
        StandardGeneticMap,
        HaldaneMapFunction,
    )


# ==========================================================================
# Genotype loading + synthetic genomic model (real pybrops objects)
# ==========================================================================


def load_pgmat(vcf_path: Path) -> Any:
    """Loads a phased genotype matrix straight from VCF via cyvcf2."""
    DensePhasedGenotypeMatrix, *_ = _import_pybrops()
    log.info("Parsing genotype panel from %s (pybrops + cyvcf2)...", vcf_path)
    pgmat = DensePhasedGenotypeMatrix.from_vcf(str(vcf_path), auto_group_vrnt=True)
    log.info(
        "Loaded %d taxa x %d variants (ploidy %d) via cyvcf2.",
        pgmat.ntaxa, pgmat.nvrnt, pgmat.ploidy,
    )
    return pgmat


def build_genomic_model(nvrnt: int, seed: int = TRAIT_SEED) -> Any:
    """
    Builds a DenseAdditiveLinearGenomicModel from the same synthetic
    negatively-correlated marker effects phase 1 uses (see
    breeding.build_synthetic_marker_effects) - no fitting step, since
    pybrops accepts hand-specified marker effect coefficients directly.
    """
    _, DenseAdditiveLinearGenomicModel, *_ = _import_pybrops()
    rng = np.random.default_rng(seed)
    u_a = rng.multivariate_normal(mean=[0.0, 0.0], cov=TRAIT_COV, size=nvrnt)
    beta = np.array([list(TRAIT_INTERCEPTS)])  # shape (q=1, t=2)
    return DenseAdditiveLinearGenomicModel(
        beta=beta,
        u_misc=None,
        u_a=u_a,
        trait=np.array(TRAIT_NAMES, dtype=object),
    )


def assign_approximate_genetic_map(pgmat: Any, cm_per_mb: float = APPROX_CM_PER_MB) -> None:
    """
    Fabricates and interpolates a flat-rate genetic map onto `pgmat` in
    place, populating `vrnt_xoprob` (per-marker crossover probability) that
    pybrops's mating protocols require. See module docstring: no real
    genetic map exists for this panel, so this is a known approximation and
    only matters for multi-generation mating simulation, not the
    single-generation OCS frontier.
    """
    *_, StandardGeneticMap, HaldaneMapFunction = _import_pybrops()
    gmap = StandardGeneticMap(
        vrnt_chrgrp=pgmat.vrnt_chrgrp,
        vrnt_phypos=pgmat.vrnt_phypos,
        vrnt_genpos=pgmat.vrnt_phypos * (cm_per_mb * 1e-8),
        vrnt_genpos_units="M",
    )
    pgmat.interp_xoprob(gmap, HaldaneMapFunction())


# ==========================================================================
# Optimal Contribution Selection (real pybrops protocol)
# ==========================================================================


def _reorder_objectives(soln_obj: np.ndarray) -> np.ndarray:
    """
    OptimalContributionSubsetSelectionProblem.latentfn returns columns
    [mean_genomic_relationship, -meanBV_trait0, -meanBV_trait1]. Reorder to
    match breeding.ParetoFrontier's convention:
    [-meanBV_cold, -meanBV_drought, mean_kinship], so phase-1's
    score_frontier/recommend_for_field/build_recommendations/
    VarietyReportGenerator can be reused unmodified.
    """
    return soln_obj[:, [1, 2, 0]]


def compute_ocs_frontier(
    pgmat: Any,
    gpmod: Any,
    k: int = 20,
    pop_size: int = 80,
    n_gen: int = 120,
    seed: int = 1,
) -> ParetoFrontier:
    (
        _,
        _,
        DenseVanRadenCoancestryMatrixFactory,
        OptimalContributionSubsetSelection,
        NSGA2SubsetGeneticAlgorithm,
        *_,
    ) = _import_pybrops()

    cmatfcty = DenseVanRadenCoancestryMatrixFactory()
    bvmat = gpmod.gebv(pgmat)

    log.info(
        "Running pybrops OptimalContributionSubsetSelection (k=%d parents out of "
        "%d taxa, pop_size=%d, n_gen=%d)... this may take a little while.",
        k, pgmat.ntaxa, pop_size, n_gen,
    )
    sel = OptimalContributionSubsetSelection(
        ntrait=len(TRAIT_NAMES),
        cmatfcty=cmatfcty,
        unscale=False,
        ncross=k,
        nparent=1,
        nmating=1,
        nprogeny=1,
        nobj=1 + len(TRAIT_NAMES),
        rng=np.random.default_rng(seed),
        moalgo=NSGA2SubsetGeneticAlgorithm(
            ngen=n_gen, pop_size=pop_size, rng=np.random.default_rng(seed)
        ),
    )
    mosoln = sel.mosolve(
        pgmat=pgmat, gmat=pgmat, ptdf=None, bvmat=bvmat, gpmod=gpmod, t_cur=0, t_max=0,
    )
    log.info("Pareto frontier computed: %d solutions.", len(mosoln.soln_obj))

    objectives = _reorder_objectives(mosoln.soln_obj)
    n_solutions = mosoln.soln_decn.shape[0]
    masks = np.zeros((n_solutions, pgmat.ntaxa), dtype=bool)
    for i, idx in enumerate(mosoln.soln_decn):
        masks[i, idx] = True

    return ParetoFrontier(taxa=list(pgmat.taxa), masks=masks, objectives=objectives, k=k)


def _params_fingerprint(vcf_path: Path, k: int, pop_size: int, n_gen: int, seed: int) -> str:
    vcf_stat = vcf_path.stat()
    key = f"pybrops:{vcf_path.name}:{vcf_stat.st_size}:{k}:{pop_size}:{n_gen}:{seed}:{TRAIT_SEED}"
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def compute_or_load_ocs_frontier(
    pgmat: Any,
    gpmod: Any,
    vcf_path: Path,
    cache_path: Path = DEFAULT_CACHE_PATH,
    k: int = 20,
    pop_size: int = 80,
    n_gen: int = 120,
    seed: int = 1,
    force_regen: bool = False,
) -> ParetoFrontier:
    fingerprint = _params_fingerprint(vcf_path, k, pop_size, n_gen, seed)

    if not force_regen and cache_path.exists():
        with np.load(cache_path, allow_pickle=True) as data:
            if str(data["fingerprint"]) == fingerprint:
                log.info("Loaded cached pybrops OCS frontier from %s", cache_path)
                return ParetoFrontier(
                    taxa=list(data["taxa"]),
                    masks=data["masks"],
                    objectives=data["objectives"],
                    k=int(data["k"]),
                )
            log.info("Cache at %s is stale (params changed); recomputing.", cache_path)

    frontier = compute_ocs_frontier(pgmat, gpmod, k=k, pop_size=pop_size, n_gen=n_gen, seed=seed)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        cache_path,
        fingerprint=fingerprint,
        taxa=np.array(frontier.taxa, dtype=object),
        masks=frontier.masks,
        objectives=frontier.objectives,
        k=frontier.k,
    )
    log.info("Cached pybrops OCS frontier to %s", cache_path)
    return frontier


# ==========================================================================
# Multi-generation breeding simulation (real pybrops mating protocol)
# ==========================================================================


@dataclass
class GenerationRecord:
    generation: int
    n_taxa: int
    n_crosses: int
    mean_bv_cold: float
    mean_bv_drought: float
    mean_kinship: float
    weights: tuple[float, float]


@dataclass
class BreedingSimResult:
    field_id: str
    field_name: str
    generations: list[GenerationRecord] = dataclass_field(default_factory=list)


def simulate_breeding_generations(
    pgmat: Any,
    gpmod: Any,
    field: Any,  # farm_report.Field
    n_generations: int = 5,
    n_crosses: int = 10,
    nmating: int = 1,
    nprogeny: int = 10,
    pop_size: int = 60,
    n_gen_ga: int = 60,
    seed: int = 1,
) -> BreedingSimResult:
    """
    Simulates `n_generations` of: OCS-select the best `n_crosses` biparental
    crosses for this field's trait preference -> mate them via TwoWayCross
    -> re-evaluate GEBVs/kinship on the resulting progeny -> repeat. Tracks
    genetic gain (mean GEBV of the selected crosses' parents) and kinship
    trend across generations, so the report can show whether the population
    is actually improving (and how quickly it's losing diversity) under
    sustained selection - the core question multi-generation breeding
    simulation exists to answer, which a single-generation frontier can't.
    """
    (
        _,
        _,
        DenseVanRadenCoancestryMatrixFactory,
        OptimalContributionSubsetSelection,
        NSGA2SubsetGeneticAlgorithm,
        TwoWayCross,
        *_,
    ) = _import_pybrops()

    cold_weight, drought_weight = breeding.field_trait_weights(
        field.thresholds["frost_temp_c"], field.thresholds["soil_moisture_min"]
    )

    cmatfcty = DenseVanRadenCoancestryMatrixFactory()
    mate_op = TwoWayCross(rng=np.random.default_rng(seed))

    result = BreedingSimResult(field_id=field.id, field_name=field.name)
    pgmat_t = pgmat

    for gen in range(n_generations):
        if pgmat_t.ntaxa < n_crosses * 2:
            log.warning(
                "Field %s: population shrank to %d taxa, too few for %d crosses "
                "of 2 parents each - stopping simulation early at generation %d.",
                field.id, pgmat_t.ntaxa, n_crosses, gen,
            )
            break

        bvmat = gpmod.gebv(pgmat_t)
        sel = OptimalContributionSubsetSelection(
            ntrait=len(TRAIT_NAMES),
            cmatfcty=cmatfcty,
            unscale=False,
            ncross=n_crosses,
            nparent=2,
            nmating=nmating,
            nprogeny=nprogeny,
            nobj=1 + len(TRAIT_NAMES),
            rng=np.random.default_rng(seed + gen),
            moalgo=NSGA2SubsetGeneticAlgorithm(
                ngen=n_gen_ga, pop_size=pop_size, rng=np.random.default_rng(seed + gen)
            ),
        )
        mosoln = sel.mosolve(
            pgmat=pgmat_t, gmat=pgmat_t, ptdf=None, bvmat=bvmat, gpmod=gpmod,
            t_cur=gen, t_max=n_generations,
        )
        objectives = _reorder_objectives(mosoln.soln_obj)
        score = breeding.score_frontier(objectives, cold_weight, drought_weight)
        best_idx = int(np.argmax(score))

        result.generations.append(
            GenerationRecord(
                generation=gen,
                n_taxa=pgmat_t.ntaxa,
                n_crosses=n_crosses,
                mean_bv_cold=float(-objectives[best_idx, 0]),
                mean_bv_drought=float(-objectives[best_idx, 1]),
                mean_kinship=float(objectives[best_idx, 2]),
                weights=(cold_weight, drought_weight),
            )
        )
        log.info(
            "Field %s gen %d: mean BV cold=%.3f drought=%.3f kinship=%.4f (n=%d taxa)",
            field.id, gen,
            result.generations[-1].mean_bv_cold,
            result.generations[-1].mean_bv_drought,
            result.generations[-1].mean_kinship,
            pgmat_t.ntaxa,
        )

        xconfig = mosoln.soln_decn[best_idx].reshape(n_crosses, 2)
        pgmat_t = mate_op.mate(pgmat_t, xconfig, nmating=nmating, nprogeny=nprogeny)

    return result


# ==========================================================================
# Report generation
# ==========================================================================


class BreedingSimReportGenerator:
    def __init__(self, farm_name: str, results: list[BreedingSimResult]):
        self.farm_name = farm_name
        self.results = results

    def build(self) -> str:
        lines: list[str] = []
        lines.append(f"# Multi-Generation Breeding Simulation — {self.farm_name}")
        lines.append(
            "_Real `pybrops` Optimal Contribution Selection + `TwoWayCross` mating, "
            "run generation-over-generation — see caveats below._\n"
        )
        lines.append(
            "> **Simulated data notice:** `cold_tolerance` / `drought_tolerance` are "
            "synthetic trait effects (see `breeding.py`), not measured phenotypes. "
            "**Genetic map notice:** no real genetic map exists for this SNP panel, "
            "so recombination during mating uses a flat 1 cM/Mb genome-wide "
            "approximation (see `breeding_pybrops.py` module docstring) rather than "
            "this panel's true recombination landscape - generation-to-generation "
            "trends here demonstrate the simulation methodology, not a validated "
            "breeding forecast.\n"
        )

        if not self.results:
            lines.append("_No corn fields found in this farm layout — nothing to simulate._")
            return "\n".join(lines)

        for res in self.results:
            lines.append(f"## {res.field_name} (`{res.field_id}`)")
            if not res.generations:
                lines.append("_Simulation produced no generations (see log)._\n")
                continue
            w = res.generations[0].weights
            lines.append(
                f"- Trait weights derived from field thresholds: "
                f"cold_tolerance {w[0]:.2f}, drought_tolerance {w[1]:.2f}"
            )
            lines.append("")
            lines.append("| Gen | Population | Mean cold BV | Mean drought BV | Mean kinship |")
            lines.append("|---|---|---|---|---|")
            for g in res.generations:
                lines.append(
                    f"| {g.generation} | {g.n_taxa} | {g.mean_bv_cold:.2f} | "
                    f"{g.mean_bv_drought:.2f} | {g.mean_kinship:.4f} |"
                )
            lines.append("")

        return "\n".join(lines)
