"""Data types for fpocket results."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FpocketParams:
    """Parameters for fpocket pocket detection.

    All values default to fpocket's built-in defaults when set to 0.
    """

    min_alpha_size: float = 0.0
    max_alpha_size: float = 0.0
    min_spheres_per_pocket: int = 0
    clustering_distance: float = 0.0
    mc_iterations: int = 0


@dataclass(frozen=True)
class Pocket:
    """A detected binding pocket."""

    rank: int
    score: float
    druggability_score: float
    num_alpha_spheres: int
    num_apolar_alpha_spheres: int
    num_polar_alpha_spheres: int
    volume: float
    convex_hull_volume: float
    hydrophobicity_score: float
    polarity_score: float
    charge_score: float
    prop_polar_atoms: float
    mean_alpha_sphere_radius: float
    alpha_sphere_density: float
    mean_local_hydrophobic_density: float
    apolar_alpha_sphere_proportion: float
    max_alpha_sphere_distance: float
    flexibility: float
    surf_vdw: float
    surf_pol_vdw: float
    surf_apol_vdw: float
    center: tuple[float, float, float] = field(default=(0.0, 0.0, 0.0))
    inter_chain: bool = False
    chain1: str = ""
    chain2: str = ""
    lig_tag: str = ""
    aa_composition: tuple[int, ...] = field(default_factory=lambda: (0,) * 20)
    n_abpa: int = 0
