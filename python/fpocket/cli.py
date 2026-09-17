"""CLI interface for fpocket, mirroring the native fpocket command."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from fpocket.binding import FpocketError, find_pockets
from fpocket.types import FpocketParams, Pocket


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fpocket",
        description="Detect binding pockets in PDB/mmCIF files.",
    )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="Input PDB or mmCIF file.",
    )
    parser.add_argument(
        "-m", "--min-alpha-size",
        type=float,
        default=0.0,
        help="Minimum radius of an alpha-sphere (default: fpocket default 3.0).",
    )
    parser.add_argument(
        "-M", "--max-alpha-size",
        type=float,
        default=0.0,
        help="Maximum radius of an alpha-sphere (default: fpocket default 6.0).",
    )
    parser.add_argument(
        "-i", "--min-spheres",
        type=int,
        default=0,
        help="Minimum number of alpha-spheres per pocket (default: fpocket default 30).",
    )
    parser.add_argument(
        "-D", "--clustering-distance",
        type=float,
        default=0.0,
        help="Distance threshold for clustering algorithm.",
    )
    parser.add_argument(
        "-v", "--mc-iterations",
        type=int,
        default=0,
        help="Number of Monte-Carlo iterations for volume calculation (default: 2500).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON.",
    )
    return parser


def _format_pocket_text(pocket: Pocket) -> str:
    lines = [
        f"Pocket {pocket.rank}:",
        f"  Score:                  {pocket.score:.4f}",
        f"  Druggability Score:     {pocket.druggability_score:.4f}",
        f"  Number of Alpha Spheres:{pocket.num_alpha_spheres:>6d}",
        f"    Apolar:               {pocket.num_apolar_alpha_spheres:>6d}",
        f"    Polar:                {pocket.num_polar_alpha_spheres:>6d}",
        f"  Volume:                 {pocket.volume:.4f}",
        f"  Convex Hull Volume:     {pocket.convex_hull_volume:.4f}",
        f"  Hydrophobicity Score:   {pocket.hydrophobicity_score:.4f}",
        f"  Polarity Score:         {pocket.polarity_score:.4f}",
        f"  Charge Score:           {pocket.charge_score:.4f}",
        f"  Flexibility:            {pocket.flexibility:.4f}",
        f"  Surface VdW:            {pocket.surf_vdw:.4f}",
        f"  Center: ({pocket.center[0]:.3f}, {pocket.center[1]:.3f}, {pocket.center[2]:.3f})",
    ]
    if pocket.inter_chain:
        lines.append(f"  Inter-chain: {pocket.chain1} - {pocket.chain2}")
    if pocket.lig_tag:
        lines.append(f"  Ligand: {pocket.lig_tag}")
    return "\n".join(lines)


def _pocket_to_dict(pocket: Pocket) -> dict:
    return {
        "rank": pocket.rank,
        "score": pocket.score,
        "druggability_score": pocket.druggability_score,
        "num_alpha_spheres": pocket.num_alpha_spheres,
        "num_apolar_alpha_spheres": pocket.num_apolar_alpha_spheres,
        "num_polar_alpha_spheres": pocket.num_polar_alpha_spheres,
        "volume": pocket.volume,
        "convex_hull_volume": pocket.convex_hull_volume,
        "hydrophobicity_score": pocket.hydrophobicity_score,
        "polarity_score": pocket.polarity_score,
        "charge_score": pocket.charge_score,
        "prop_polar_atoms": pocket.prop_polar_atoms,
        "mean_alpha_sphere_radius": pocket.mean_alpha_sphere_radius,
        "alpha_sphere_density": pocket.alpha_sphere_density,
        "mean_local_hydrophobic_density": pocket.mean_local_hydrophobic_density,
        "apolar_alpha_sphere_proportion": pocket.apolar_alpha_sphere_proportion,
        "max_alpha_sphere_distance": pocket.max_alpha_sphere_distance,
        "flexibility": pocket.flexibility,
        "surf_vdw": pocket.surf_vdw,
        "surf_pol_vdw": pocket.surf_pol_vdw,
        "surf_apol_vdw": pocket.surf_apol_vdw,
        "center": list(pocket.center),
        "inter_chain": pocket.inter_chain,
        "chain1": pocket.chain1,
        "chain2": pocket.chain2,
        "lig_tag": pocket.lig_tag,
        "aa_composition": list(pocket.aa_composition),
        "n_abpa": pocket.n_abpa,
    }


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    pdb_path = Path(args.file)
    if not pdb_path.exists():
        print(f"Error: file not found: {pdb_path}", file=sys.stderr)
        return 1

    params = FpocketParams(
        min_alpha_size=args.min_alpha_size,
        max_alpha_size=args.max_alpha_size,
        min_spheres_per_pocket=args.min_spheres,
        clustering_distance=args.clustering_distance,
        mc_iterations=args.mc_iterations,
    )

    try:
        pockets = find_pockets(pdb_path, params)
    except FpocketError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if args.json:
        output = {
            "file": str(pdb_path),
            "num_pockets": len(pockets),
            "pockets": [_pocket_to_dict(p) for p in pockets],
        }
        json.dump(output, sys.stdout, indent=2)
        print()
    else:
        print(f"\nfpocket - {pdb_path.name}")
        print(f"Found {len(pockets)} pocket(s)\n")
        for pocket in pockets:
            print(_format_pocket_text(pocket))
            print()

    return 0


def _entry() -> None:
    sys.exit(main())


if __name__ == "__main__":
    _entry()
