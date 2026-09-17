"""CFFI bindings to the fpocket C wrapper library."""

from __future__ import annotations

import os
from pathlib import Path

import cffi

from fpocket.types import FpocketParams, Pocket

HEADER = """
typedef struct {
    int rank;
    float score;
    float druggability_score;
    int num_alpha_spheres;
    int num_apolar_alpha_spheres;
    int num_polar_alpha_spheres;
    float volume;
    float convex_hull_volume;
    float hydrophobicity_score;
    float polarity_score;
    float charge_score;
    float prop_polar_atoms;
    float mean_alpha_sphere_radius;
    float alpha_sphere_density;
    float mean_local_hydrophobic_density;
    float apolar_alpha_sphere_proportion;
    float max_alpha_sphere_distance;
    float flexibility;
    float surf_vdw;
    float surf_pol_vdw;
    float surf_apol_vdw;
    float bary_x;
    float bary_y;
    float bary_z;
    int inter_chain;
    char chain1[256];
    char chain2[256];
    char lig_tag[8];
    int aa_composition[20];
    int n_abpa;
} fpocket_pocket_result;

typedef struct {
    fpocket_pocket_result *pockets;
    int num_pockets;
    int error_code;
    char error_message[512];
} fpocket_result;

fpocket_result *fpocket_run(const char *pdb_path,
                             float min_alpha_size,
                             float max_alpha_size,
                             int min_spheres_per_pocket,
                             float clustering_distance,
                             int mc_iterations);

void fpocket_free_result(fpocket_result *result);
"""

ffi = cffi.FFI()
ffi.cdef(HEADER)


def _load_lib() -> cffi.FFI:
    """Load the fpocket_wrapper shared library."""
    pkg_dir = Path(__file__).parent
    # Look for the library in common locations
    for suffix in (".so", ".dylib"):
        for search_dir in [pkg_dir, pkg_dir / ".."]:
            lib_path = search_dir / f"libfpocket_wrapper{suffix}"
            if lib_path.exists():
                return ffi.dlopen(str(lib_path))

    # Fall back to system library path
    lib_name = "libfpocket_wrapper.so"
    if os.uname().sysname == "Darwin":
        lib_name = "libfpocket_wrapper.dylib"
    return ffi.dlopen(lib_name)


_lib = _load_lib()


class FpocketError(RuntimeError):
    """Raised when fpocket encounters an error during pocket detection."""


def find_pockets(
    pdb_path: str | Path,
    params: FpocketParams | None = None,
) -> list[Pocket]:
    """Detect binding pockets in a PDB or mmCIF file.

    Args:
        pdb_path: Path to a .pdb or .cif file.
        params: Optional parameters for pocket detection.
            When None, fpocket defaults are used.

    Returns:
        List of detected pockets, sorted by score (best first).

    Raises:
        FileNotFoundError: If the PDB file does not exist.
        FpocketError: If fpocket encounters an internal error.
    """
    pdb_path = Path(pdb_path).resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if params is None:
        params = FpocketParams()

    result = _lib.fpocket_run(
        str(pdb_path).encode(),
        params.min_alpha_size,
        params.max_alpha_size,
        params.min_spheres_per_pocket,
        params.clustering_distance,
        params.mc_iterations,
    )

    if result == ffi.NULL:
        raise FpocketError("fpocket_run returned NULL")

    try:
        if result.error_code != 0:
            msg = ffi.string(result.error_message).decode()
            raise FpocketError(f"fpocket error ({result.error_code}): {msg}")

        pockets = []
        for i in range(result.num_pockets):
            p = result.pockets[i]
            pockets.append(
                Pocket(
                    rank=p.rank,
                    score=p.score,
                    druggability_score=p.druggability_score,
                    num_alpha_spheres=p.num_alpha_spheres,
                    num_apolar_alpha_spheres=p.num_apolar_alpha_spheres,
                    num_polar_alpha_spheres=p.num_polar_alpha_spheres,
                    volume=p.volume,
                    convex_hull_volume=p.convex_hull_volume,
                    hydrophobicity_score=p.hydrophobicity_score,
                    polarity_score=p.polarity_score,
                    charge_score=p.charge_score,
                    prop_polar_atoms=p.prop_polar_atoms,
                    mean_alpha_sphere_radius=p.mean_alpha_sphere_radius,
                    alpha_sphere_density=p.alpha_sphere_density,
                    mean_local_hydrophobic_density=p.mean_local_hydrophobic_density,
                    apolar_alpha_sphere_proportion=p.apolar_alpha_sphere_proportion,
                    max_alpha_sphere_distance=p.max_alpha_sphere_distance,
                    flexibility=p.flexibility,
                    surf_vdw=p.surf_vdw,
                    surf_pol_vdw=p.surf_pol_vdw,
                    surf_apol_vdw=p.surf_apol_vdw,
                    center=(p.bary_x, p.bary_y, p.bary_z),
                    inter_chain=bool(p.inter_chain),
                    chain1=ffi.string(p.chain1).decode(),
                    chain2=ffi.string(p.chain2).decode(),
                    lig_tag=ffi.string(p.lig_tag).decode(),
                    aa_composition=tuple(p.aa_composition[j] for j in range(20)),
                    n_abpa=p.n_abpa,
                )
            )
        return pockets
    finally:
        _lib.fpocket_free_result(result)
