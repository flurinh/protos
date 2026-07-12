#!/usr/bin/env python3
"""RFdiffusion2 configuration builder for protein design.

Provides a flexible interface for building RFdiffusion2 configs with:
- Fixed backbone regions (preserve backbone frames - indexed)
- Designable regions (new backbone generation)
- Atomic constraints (preserve specific atoms - can be indexed or unindexed/guidepost)
- Partial diffusion (refinement of existing structures)
- Ligand conditioning with RASA control

Based on: "Atom level enzyme active site scaffolding using RFdiffusion2"
(Ahern et al., 2025)

Key concepts:
- **Backbone Motif**: Fixed backbone frames with known sequence indices
- **Atomic Motif**: Fixed atom positions, model infers rotamer
- **Unindexed/Guidepost**: Model infers both sequence index AND rotamer
- **Partial Diffusion**: Refine existing structure with limited noise
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Literal
from enum import Enum
import json
import subprocess


def cif_to_pdb(cif_path: Path, pdb_path: Optional[Path] = None) -> Path:
    """Convert CIF to PDB format using gemmi.

    RFdiffusion2 requires PDB format input. This converts CIF files
    (e.g., from BoltzGen) to PDB.

    Args:
        cif_path: Path to input CIF file
        pdb_path: Output PDB path (default: same name with .pdb extension)

    Returns:
        Path to output PDB file
    """
    cif_path = Path(cif_path)
    if pdb_path is None:
        pdb_path = cif_path.with_suffix(".pdb")
    else:
        pdb_path = Path(pdb_path)

    try:
        import gemmi
        structure = gemmi.read_structure(str(cif_path))
        structure.write_pdb(str(pdb_path))
    except ImportError:
        # Fallback to command-line gemmi if available
        result = subprocess.run(
            ["gemmi", "convert", str(cif_path), str(pdb_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Failed to convert {cif_path} to PDB. "
                f"Install gemmi: pip install gemmi\n{result.stderr}"
            )

    return pdb_path


class MotifType(Enum):
    """Type of motif representation for RFdiffusion2.

    From the paper (Figure 2B):
    - BACKBONE: Fixed backbone frames with known indices (like RFdiffusion v1)
    - ATOMIC_INDEXED: Fixed atom positions with known indices, model infers rotamer
    - ATOMIC_GUIDEPOST: Fixed atom positions, model infers BOTH index AND rotamer

    For rhodozyme design with guidepost=True:
    - All atom positions are preserved at their 3D coordinates
    - Sequence indices can be inferred by the model (unindexed mode)
    - This allows the model to optimize where in the sequence the motif appears
    """
    BACKBONE = "backbone"           # Indexed backbone frames
    ATOMIC_INDEXED = "atomic"       # Indexed atomic constraints
    ATOMIC_GUIDEPOST = "guidepost"  # Unindexed atomic (model infers index)


@dataclass
class AtomConstraint:
    """Atomic constraint for motif scaffolding.

    Specifies atoms that should be preserved during diffusion.
    """
    chain: str
    resnum: int
    atoms: List[str]

    @property
    def key(self) -> str:
        """Return chain+resnum key (e.g., 'A230')."""
        return f"{self.chain}{self.resnum}"

    @classmethod
    def from_dict(cls, key: str, atoms: List[str]) -> "AtomConstraint":
        """Create from dict key format (e.g., 'A230')."""
        chain = key[0]
        resnum = int(key[1:])
        return cls(chain=chain, resnum=resnum, atoms=atoms)


@dataclass
class Region:
    """A contiguous region in the protein.

    Can be:
    - Fixed backbone (BACKBONE): Keep backbone frames from input, indexed
    - Designable (new backbone generation)

    For the user's rhodozyme case:
    - TM helices are BACKBONE motifs (fixed frames, known indices)
    - Catalytic triad atoms are ATOMIC constraints (separate from regions)
    - Loops between helices are DESIGNABLE regions
    """
    start: Optional[int] = None  # None for designable regions
    end: Optional[int] = None
    chain: str = "A"
    length: Optional[int] = None  # For designable: exact length
    length_range: Optional[tuple] = None  # For designable: (min, max)
    fixed: bool = True  # True = keep from input, False = design new
    motif_type: MotifType = MotifType.BACKBONE  # Only relevant if fixed=True

    def to_contig(self) -> str:
        """Convert to contig string format."""
        if self.fixed:
            if self.start == self.end:
                return f"{self.chain}{self.start}-{self.start}"
            return f"{self.chain}{self.start}-{self.end}"
        else:
            if self.length is not None:
                return str(self.length)
            elif self.length_range is not None:
                return f"{self.length_range[0]}-{self.length_range[1]}"
            else:
                raise ValueError("Designable region needs length or length_range")


@dataclass
class RFD2ConfigBuilder:
    """Flexible configuration builder for RFdiffusion2.

    Example usage:
        builder = RFD2ConfigBuilder(
            input_pdb="structure.pdb",
            output_dir="output/",
        )

        # Define structure: fixed helix, designable loop, fixed helix
        builder.add_fixed_region(chain="A", start=1, end=50)
        builder.add_designable_region(length_range=(10, 20))
        builder.add_fixed_region(chain="A", start=71, end=120)

        # Add atomic constraints for catalytic residues
        builder.add_atom_constraint("A", 45, ["OG", "CB", "CA"])  # SER
        builder.add_atom_constraint("A", 80, ["NE2", "ND1"])      # HIS

        # Build config
        config = builder.build()
    """

    input_pdb: Path
    output_dir: Path
    job_name: str = "rfd2_job"

    # Design parameters
    num_designs: int = 10

    # Regions and constraints
    regions: List[Region] = field(default_factory=list)
    atom_constraints: List[AtomConstraint] = field(default_factory=list)

    # Ligands
    ligands: Optional[List[str]] = None
    partially_fixed_ligand: Optional[Dict[str, List[str]]] = None

    # Diffusion parameters
    partial_T: Optional[int] = None  # None = full diffusion
    total_T: int = 100

    # Scaffolding mode
    guidepost: bool = True  # Unindexed scaffolding (recommended for enzymes)

    # Pipeline control
    stop_step: str = "mpnn"  # 'sweep', 'mpnn', 'thread_mpnn', 'score', 'end'

    # Model selection
    checkpoint: str = "RFD_173"  # 'RFD_173' or 'RFD_140'

    # Additional parameters
    deterministic: bool = True
    seed_offset: int = 0
    extra_params: Dict[str, Any] = field(default_factory=dict)

    def add_fixed_region(
        self,
        start: int,
        end: int,
        chain: str = "A",
    ) -> "RFD2ConfigBuilder":
        """Add a fixed region from the input structure.

        Args:
            start: Start residue number
            end: End residue number
            chain: Chain identifier
        """
        self.regions.append(Region(
            start=start,
            end=end,
            chain=chain,
            fixed=True,
        ))
        return self

    def add_fixed_residue(
        self,
        resnum: int,
        chain: str = "A",
    ) -> "RFD2ConfigBuilder":
        """Add a single fixed residue."""
        return self.add_fixed_region(resnum, resnum, chain)

    def add_designable_region(
        self,
        length: Optional[int] = None,
        length_range: Optional[tuple] = None,
    ) -> "RFD2ConfigBuilder":
        """Add a designable region (new backbone).

        Args:
            length: Exact length of new region
            length_range: (min, max) length range
        """
        if length is None and length_range is None:
            raise ValueError("Must specify length or length_range")

        self.regions.append(Region(
            length=length,
            length_range=length_range,
            fixed=False,
        ))
        return self

    def add_atom_constraint(
        self,
        chain: str,
        resnum: int,
        atoms: List[str],
    ) -> "RFD2ConfigBuilder":
        """Add atomic constraint for motif scaffolding.

        Args:
            chain: Chain identifier
            resnum: Residue number
            atoms: List of atom names to preserve (e.g., ['OG', 'CB', 'CA'])
        """
        self.atom_constraints.append(AtomConstraint(
            chain=chain,
            resnum=resnum,
            atoms=atoms,
        ))
        return self

    def add_catalytic_triad(
        self,
        ser_pos: int,
        his_pos: int,
        asp_pos: int,
        chain: str = "A",
    ) -> "RFD2ConfigBuilder":
        """Add standard catalytic triad constraints.

        Preserves key atoms for SER-HIS-ASP catalytic triad.
        """
        # SER: hydroxyl oxygen + backbone
        self.add_atom_constraint(chain, ser_pos, ["OG", "CB", "CA"])
        # HIS: imidazole ring
        self.add_atom_constraint(chain, his_pos, ["NE2", "ND1", "CE1", "CD2", "CG", "CB", "CA"])
        # ASP: carboxyl group
        self.add_atom_constraint(chain, asp_pos, ["OD1", "OD2", "CG", "CB", "CA"])
        return self

    def add_ligand(
        self,
        name: str,
        fixed_atoms: Optional[List[str]] = None,
    ) -> "RFD2ConfigBuilder":
        """Add a ligand to the design.

        Args:
            name: Ligand name (3-letter code)
            fixed_atoms: Atoms to keep fixed during diffusion (optional)
        """
        if self.ligands is None:
            self.ligands = []
        self.ligands.append(name)

        if fixed_atoms:
            if self.partially_fixed_ligand is None:
                self.partially_fixed_ligand = {}
            self.partially_fixed_ligand[name] = fixed_atoms

        return self

    def set_partial_diffusion(
        self,
        partial_T: int = 20,
        total_T: int = 100,
    ) -> "RFD2ConfigBuilder":
        """Enable partial diffusion (refinement mode).

        Args:
            partial_T: Number of diffusion steps (lower = less change)
            total_T: Total diffusion timesteps
        """
        self.partial_T = partial_T
        self.total_T = total_T
        return self

    def set_full_chain(
        self,
        chain: str = "A",
        length: int = 350,
    ) -> "RFD2ConfigBuilder":
        """Set a single fixed chain (for partial diffusion).

        Use this when doing partial diffusion on an entire structure.
        """
        self.regions = [Region(start=1, end=length, chain=chain, fixed=True)]
        return self

    def _build_contigs(self) -> str:
        """Build contigs string from regions."""
        if not self.regions:
            raise ValueError("No regions defined. Use add_fixed_region/add_designable_region")

        parts = [r.to_contig() for r in self.regions]
        return ",".join(parts)

    def _build_contig_atoms(self) -> Optional[Dict[str, List[str]]]:
        """Build contig_atoms dict from constraints."""
        if not self.atom_constraints:
            return None

        return {c.key: c.atoms for c in self.atom_constraints}

    def _build_ligand_str(self) -> Optional[str]:
        """Build ligand specification string."""
        if not self.ligands:
            return None
        return ",".join(self.ligands)

    def _get_checkpoint_path(self) -> str:
        """Get full checkpoint path."""
        if self.checkpoint.endswith(".pt"):
            return self.checkpoint
        return f"rf_diffusion/model_weights/{self.checkpoint}.pt"

    def build(self) -> Dict[str, Any]:
        """Build the final configuration dict.

        Returns a configuration dictionary for RFdiffusion2 tooling.
        """
        config = {
            "job_name": self.job_name,
            "input_pdb": str(self.input_pdb),
            "output_dir": str(self.output_dir),
            "contigs": self._build_contigs(),
            "num_designs": self.num_designs,
            "guidepost": self.guidepost,
            "stop_step": self.stop_step,
        }

        # Atomic constraints
        contig_atoms = self._build_contig_atoms()
        if contig_atoms:
            config["contig_atoms"] = contig_atoms

        # Ligands
        ligand_str = self._build_ligand_str()
        if ligand_str:
            config["ligands"] = ligand_str

        # Partial diffusion
        if self.partial_T is not None:
            config["partial_T"] = self.partial_T

        # Extra inference parameters
        extra = {
            "inference.deterministic": self.deterministic,
            "inference.seed_offset": self.seed_offset,
            "inference.ckpt_path": self._get_checkpoint_path(),
            "diffuser.T": self.total_T,
        }

        # Partially fixed ligand
        if self.partially_fixed_ligand:
            lig_str = json.dumps(self.partially_fixed_ligand).replace('"', '\\"')
            extra["inference.partially_fixed_ligand"] = f'"{lig_str}"'

        # Merge with user extra params
        extra.update(self.extra_params)
        config["extra_params"] = extra

        return config

    def to_hydra_args(self) -> List[str]:
        """Convert to Hydra command-line arguments.

        Useful for direct command-line invocation.
        """
        config = self.build()
        args = []

        # Input PDB
        args.append(f"inference.input_pdb={config['input_pdb']}")

        # Contigs
        args.append(f"contigmap.contigs=['{config['contigs']}']")

        # Num designs
        args.append(f"inference.num_designs={config['num_designs']}")

        # Guidepost
        args.append(f"inference.contig_as_guidepost={config['guidepost']}")

        # Contig atoms
        if "contig_atoms" in config:
            atoms_str = json.dumps(config["contig_atoms"]).replace('"', "'")
            args.append(f'contigmap.contig_atoms="{atoms_str}"')

        # Ligands
        if "ligands" in config:
            args.append(f"inference.ligand='{config['ligands']}'")

        # Partial T
        if "partial_T" in config:
            args.append(f"diffuser.partial_T={config['partial_T']}")

        # Stop step
        args.append(f"stop_step='{config['stop_step']}'")

        # Extra params
        for key, value in config.get("extra_params", {}).items():
            if isinstance(value, bool):
                args.append(f"{key}={str(value)}")
            elif isinstance(value, str) and not value.startswith('"'):
                args.append(f"{key}='{value}'")
            else:
                args.append(f"{key}={value}")

        return args


# ============================================================
# Convenience functions
# ============================================================

def partial_diffusion_config(
    input_pdb: Path,
    output_dir: Path,
    atom_constraints: Dict[str, List[str]],
    partial_T: int = 20,
    num_designs: int = 10,
    ligands: Optional[List[str]] = None,
    chain: str = "A",
    chain_length: int = 350,
    **kwargs,
) -> Dict[str, Any]:
    """Quick config for partial diffusion refinement.

    Args:
        input_pdb: Input structure
        output_dir: Output directory
        atom_constraints: Dict of {chain+resnum: [atoms]} to preserve
        partial_T: Diffusion steps (lower = more conservative)
        num_designs: Number of designs
        ligands: List of ligand names
        chain: Chain identifier
        chain_length: Length of the chain
        **kwargs: Additional parameters for RFD2ConfigBuilder

    Returns:
        RFdiffusion2 configuration dictionary
    """
    builder = RFD2ConfigBuilder(
        input_pdb=input_pdb,
        output_dir=output_dir,
        num_designs=num_designs,
        **kwargs,
    )

    # Set full chain for partial diffusion
    builder.set_full_chain(chain, chain_length)
    builder.set_partial_diffusion(partial_T)

    # Add atom constraints
    for key, atoms in atom_constraints.items():
        constraint = AtomConstraint.from_dict(key, atoms)
        builder.atom_constraints.append(constraint)

    # Add ligands
    if ligands:
        for lig in ligands:
            builder.add_ligand(lig)

    return builder.build()


def scaffolding_config(
    input_pdb: Path,
    output_dir: Path,
    fixed_residues: List[tuple],  # [(chain, resnum), ...]
    atom_constraints: Dict[str, List[str]],
    target_length: int = 180,
    num_designs: int = 100,
    ligands: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Any]:
    """Quick config for motif scaffolding (de novo design around fixed residues).

    Args:
        input_pdb: Input structure with motif
        output_dir: Output directory
        fixed_residues: List of (chain, resnum) tuples to keep
        atom_constraints: Dict of {chain+resnum: [atoms]} to preserve
        target_length: Target protein length
        num_designs: Number of designs
        ligands: List of ligand names
        **kwargs: Additional parameters

    Returns:
        RFdiffusion2 configuration dictionary
    """
    builder = RFD2ConfigBuilder(
        input_pdb=input_pdb,
        output_dir=output_dir,
        num_designs=num_designs,
        **kwargs,
    )

    # Sort fixed residues by position
    sorted_residues = sorted(fixed_residues, key=lambda x: x[1])

    # Build regions: gap, fixed, gap, fixed, ...
    prev_pos = 0
    remaining = target_length

    for chain, resnum in sorted_residues:
        gap = resnum - prev_pos - 1
        if gap > 0:
            builder.add_designable_region(length=gap)
            remaining -= gap
        builder.add_fixed_residue(resnum, chain)
        remaining -= 1
        prev_pos = resnum

    # Final gap
    if remaining > 0:
        builder.add_designable_region(length=remaining)

    # Add atom constraints
    for key, atoms in atom_constraints.items():
        constraint = AtomConstraint.from_dict(key, atoms)
        builder.atom_constraints.append(constraint)

    # Add ligands
    if ligands:
        for lig in ligands:
            builder.add_ligand(lig)

    return builder.build()


def hybrid_scaffold_config(
    input_pdb: Path,
    output_dir: Path,
    fixed_backbone_regions: List[tuple],  # [(chain, start, end), ...]
    designable_regions: List[tuple],  # [(length,) or (min_len, max_len), ...]
    atom_constraints: Dict[str, List[str]],
    num_designs: int = 100,
    partial_T: Optional[int] = None,
    ligands: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Any]:
    """Config for hybrid scaffolding: fixed backbone + atomic constraints + designable loops.

    This is the recommended approach for the rhodozyme case:
    - TM helices: fixed backbone frames (indexed backbone motifs)
    - Catalytic triad: atomic constraints (preserved atoms)
    - Loops: designable regions (new backbone generation)

    The order of regions should be: fixed1, designable1, fixed2, designable2, ...
    matching the N-to-C sequence of the protein.

    Args:
        input_pdb: Input structure (e.g., BoltzGen placement)
        output_dir: Output directory
        fixed_backbone_regions: List of (chain, start, end) tuples for fixed regions
        designable_regions: List of length specs - (length,) for exact or (min, max) for range
        atom_constraints: Dict of {chain+resnum: [atoms]} for catalytic residues
        num_designs: Number of designs
        partial_T: If set, use partial diffusion (refinement mode)
        ligands: List of ligand names (e.g., ['RET'])
        **kwargs: Additional parameters for RFD2ConfigBuilder

    Returns:
        RFdiffusion2 configuration dictionary

    Example:
        # Rhodozyme with 8 TM helices and connecting loops
        config = hybrid_scaffold_config(
            input_pdb=Path("placement_00.pdb"),
            output_dir=Path("output/"),
            fixed_backbone_regions=[
                ("A", 1, 53),    # TM1
                ("A", 81, 129),  # TM2
                ("A", 162, 217), # TM3-4
                ("A", 261, 301), # TM5
                ("A", 327, 350), # TM6-7
            ],
            designable_regions=[
                (10, 30),  # Loop 1 (variable length)
                (10, 35),  # Loop 2
                (15, 45),  # Loop 3
                (10, 30),  # Loop 4
            ],
            atom_constraints={
                "A230": ["OG", "CB", "CA"],      # SER
                "A245": ["NE2", "ND1", "CB"],    # HIS
                "A250": ["OD1", "OD2", "CB"],    # ASP
            },
            ligands=["RET"],
        )
    """
    builder = RFD2ConfigBuilder(
        input_pdb=input_pdb,
        output_dir=output_dir,
        num_designs=num_designs,
        **kwargs,
    )

    # Interleave fixed and designable regions
    # Expecting: fixed[0], designable[0], fixed[1], designable[1], ...
    n_fixed = len(fixed_backbone_regions)
    n_design = len(designable_regions)

    for i in range(max(n_fixed, n_design)):
        # Add fixed region if available
        if i < n_fixed:
            chain, start, end = fixed_backbone_regions[i]
            builder.add_fixed_region(start, end, chain)

        # Add designable region if available
        if i < n_design:
            spec = designable_regions[i]
            if len(spec) == 1:
                builder.add_designable_region(length=spec[0])
            else:
                builder.add_designable_region(length_range=(spec[0], spec[1]))

    # Add atomic constraints
    for key, atoms in atom_constraints.items():
        constraint = AtomConstraint.from_dict(key, atoms)
        builder.atom_constraints.append(constraint)

    # Add ligands
    if ligands:
        for lig in ligands:
            builder.add_ligand(lig)

    # Set partial diffusion if requested
    if partial_T is not None:
        builder.set_partial_diffusion(partial_T)

    return builder.build()


def rhodozyme_config(
    placement_pdb: Path,
    output_dir: Path,
    triad_positions: Dict[str, int],
    tm_helix_regions: List[tuple],
    loop_length_ranges: List[tuple],
    num_designs: int = 100,
    partial_T: int = 20,
    include_retinal: bool = True,
    **kwargs,
) -> Dict[str, Any]:
    """Quick config builder for rhodozyme design.

    This wraps hybrid_scaffold_config with sensible defaults for rhodozyme.

    NOTE: RFdiffusion2 requires PDB format. If your placement is CIF,
    convert first: pdb_path = cif_to_pdb(cif_path)

    With guidepost=True (default), the model preserves:
    - 3D positions of fixed backbone regions (TM helices)
    - 3D positions of catalytic triad atoms
    But allows the model to infer optimal sequence indices.

    Args:
        placement_pdb: Input structure in PDB format (use cif_to_pdb() if needed)
        output_dir: Output directory
        triad_positions: Dict mapping {'SER': pos, 'HIS': pos, 'ASP': pos}
        tm_helix_regions: List of (chain, start, end) for TM helices
        loop_length_ranges: List of (min_len, max_len) for loops between helices
        num_designs: Number of designs to generate
        partial_T: Partial diffusion steps (lower = more conservative)
        include_retinal: Include retinal ligand
        **kwargs: Additional parameters

    Returns:
        RFdiffusion2 configuration dictionary

    Example:
        # Convert CIF to PDB first if needed
        pdb_path = cif_to_pdb(Path("placement_00.cif"))

        config = rhodozyme_config(
            placement_pdb=pdb_path,
            output_dir=Path("output/"),
            triad_positions={"SER": 230, "HIS": 245, "ASP": 250},
            tm_helix_regions=[
                ("A", 1, 53), ("A", 81, 129), ("A", 162, 217),
                ("A", 261, 301), ("A", 327, 350),
            ],
            loop_length_ranges=[(10, 30), (10, 35), (15, 45), (10, 30)],
        )
    """
    # Build atomic constraints for catalytic triad
    atom_constraints = {}
    if "SER" in triad_positions:
        atom_constraints[f"A{triad_positions['SER']}"] = ["OG", "CB", "CA"]
    if "HIS" in triad_positions:
        atom_constraints[f"A{triad_positions['HIS']}"] = ["NE2", "ND1", "CE1", "CD2", "CG", "CB", "CA"]
    if "ASP" in triad_positions:
        atom_constraints[f"A{triad_positions['ASP']}"] = ["OD1", "OD2", "CG", "CB", "CA"]

    ligands = ["RET"] if include_retinal else None

    return hybrid_scaffold_config(
        input_pdb=placement_pdb,
        output_dir=output_dir,
        fixed_backbone_regions=tm_helix_regions,
        designable_regions=loop_length_ranges,
        atom_constraints=atom_constraints,
        num_designs=num_designs,
        partial_T=partial_T,
        ligands=ligands,
        **kwargs,
    )


# ============================================================
# Usage Examples
# ============================================================

if __name__ == "__main__":
    from pathlib import Path

    print("=" * 60)
    print("Example 1: Partial diffusion refinement")
    print("=" * 60)

    builder = RFD2ConfigBuilder(
        input_pdb=Path("placement_00.pdb"),
        output_dir=Path("output/refinement"),
        job_name="rhodozyme_refine",
        num_designs=100,
    )

    # Full chain with partial diffusion
    builder.set_full_chain("A", 350)
    builder.set_partial_diffusion(partial_T=20)

    # Catalytic triad constraints
    builder.add_catalytic_triad(ser_pos=230, his_pos=245, asp_pos=250)

    # Retinal ligand
    builder.add_ligand("RET")

    config = builder.build()
    print(f"Config: {json.dumps(config, indent=2)}")
    print()

    print("=" * 60)
    print("Example 2: De novo scaffolding around catalytic residues")
    print("=" * 60)

    builder2 = RFD2ConfigBuilder(
        input_pdb=Path("motif.pdb"),
        output_dir=Path("output/scaffold"),
        job_name="enzyme_scaffold",
        num_designs=100,
    )

    # Design new backbone around fixed catalytic residues
    builder2.add_designable_region(length=46)
    builder2.add_fixed_residue(106, "A")
    builder2.add_designable_region(length=59)
    builder2.add_fixed_residue(166, "A")
    builder2.add_designable_region(length_range=(2, 5))
    builder2.add_fixed_residue(169, "A")
    builder2.add_designable_region(length=46)

    # Atom constraints for catalytic residues
    builder2.add_atom_constraint("A", 106, ["NE", "CD", "CZ"])
    builder2.add_atom_constraint("A", 166, ["OD1", "CG"])
    builder2.add_atom_constraint("A", 169, ["NH2", "CZ"])

    # Ligands
    builder2.add_ligand("NAD", fixed_atoms=["O7N", "C7N", "C3N"])
    builder2.add_ligand("OXM")

    config2 = builder2.build()
    print(f"Config: {json.dumps(config2, indent=2)}")
    print()
    print(f"Hydra args: {' '.join(builder2.to_hydra_args())}")

    print()
    print("=" * 60)
    print("Example 3: Quick partial diffusion helper")
    print("=" * 60)

    config3 = partial_diffusion_config(
        input_pdb=Path("structure.pdb"),
        output_dir=Path("output/"),
        atom_constraints={
            "A230": ["OG", "CB", "CA"],
            "A245": ["NE2", "ND1"],
            "A250": ["OD1", "OD2"],
        },
        partial_T=25,
        num_designs=50,
        ligands=["RET"],
        job_name="quick_refine",
    )
    print(f"Config: {json.dumps(config3, indent=2)}")

    print()
    print("=" * 60)
    print("Example 4: Rhodozyme hybrid case")
    print("(Fixed TM helices + catalytic triad + designable loops)")
    print("=" * 60)

    config4 = rhodozyme_config(
        placement_pdb=Path("placement_00.pdb"),
        output_dir=Path("output/rhodozyme"),
        triad_positions={"SER": 230, "HIS": 245, "ASP": 250},
        tm_helix_regions=[
            ("A", 1, 53),     # TM1
            ("A", 81, 129),   # TM2
            ("A", 162, 217),  # TM3-4
            ("A", 261, 301),  # TM5
            ("A", 327, 350),  # TM6-7
        ],
        loop_length_ranges=[
            (10, 30),  # Loop between TM1-TM2
            (10, 35),  # Loop between TM2-TM3/4
            (15, 45),  # Loop between TM3/4-TM5
            (10, 30),  # Loop between TM5-TM6/7
        ],
        num_designs=100,
        partial_T=20,
        include_retinal=True,
        job_name="rhodozyme_hybrid",
    )
    print(f"Config: {json.dumps(config4, indent=2)}")

    print()
    print("=" * 60)
    print("Example 5: Using hybrid_scaffold_config directly")
    print("=" * 60)

    config5 = hybrid_scaffold_config(
        input_pdb=Path("placement_01.pdb"),
        output_dir=Path("output/hybrid"),
        fixed_backbone_regions=[
            ("A", 1, 50),
            ("A", 75, 130),
            ("A", 165, 220),
        ],
        designable_regions=[
            (10, 25),  # Variable length loop
            (15, 40),  # Another variable length loop
        ],
        atom_constraints={
            "A134": ["OG", "CB", "CA"],       # SER
            "A138": ["NE2", "ND1", "CB"],     # HIS
            "A226": ["OD1", "OD2", "CB"],     # ASP
        },
        num_designs=50,
        partial_T=25,
        ligands=["RET"],
        job_name="hybrid_scaffold",
    )
    print(f"Config: {json.dumps(config5, indent=2)}")
