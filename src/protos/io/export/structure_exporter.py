"""Structure exporter implementation.

Handles filesystem exports for StructureProcessor instances while keeping the
processor itself path-agnostic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union, List, Dict

import pandas as pd

from protos.io.core.base_exporter import BaseExporter
from protos.io.formats.cif_utils import write_cif_file
from protos.io.formats.molecule_utils import canonical_dataframe_to_mol
from protos.io.paths import get_structure_entity_path, get_structure_dataset_path, to_data_relative_path
try:  # Optional RDKit dependency
    from rdkit import Chem
except Exception:  # pragma: no cover
    Chem = None


class StructureExporter(BaseExporter):
    """Exporter that serializes structures to CIF files."""

    def export_entity(
        self,
        name: str,
        out_path: Optional[Path] = None,
        format: Optional[str] = None,
        overwrite: bool = True,
        chains: Optional[list[str]] = None,
        # SDF ligand controls
        ligand_group_by: Optional[Union[str, list[str]]] = None,
        include_chains: Optional[list[str]] = None,
        include_comp_ids: Optional[list[str]] = None,
        include_res_ids: Optional[list[str]] = None,
        exclude_res_ids: Optional[list[str]] = None,
        min_atoms: int = 1,
    ) -> Path:
        export_format = (format or "cif").lower()

        df = self._get_structure(name)
        if chains:
            df = self._filter_by_chains(df, chains)

        # Compute default path within Protos paths if not provided
        target_path = Path(out_path) if out_path is not None else get_structure_entity_path(
            self.processor.paths,
            name,
            extension=export_format,
        )
        if target_path.exists() and not overwrite:
            raise FileExistsError(
                f"File {target_path} already exists. Use overwrite=True to replace it."
            )

        target_path.parent.mkdir(parents=True, exist_ok=True)

        if export_format == "cif":
            write_cif_file(str(target_path), df, force_overwrite=overwrite)
        elif export_format == "pdb":
            # Export PDB using gemmi by converting from an existing CIF if available
            try:
                import gemmi  # type: ignore
            except Exception as exc:  # pragma: no cover
                raise RuntimeError("gemmi is required for PDB export") from exc

            # Prefer existing CIF on disk to ensure atom naming/fields
            cif_path = get_structure_entity_path(self.processor.paths, name, extension='cif')
            if not cif_path.exists():
                # As a fallback, write a temporary CIF from DataFrame and convert
                tmp_cif = target_path.with_suffix('.tmp.cif')
                write_cif_file(str(tmp_cif), df, force_overwrite=True)
                cif_path = tmp_cif
                cleanup_tmp = True
            else:
                cleanup_tmp = False

            st = gemmi.read_structure(str(cif_path))
            st.write_pdb(str(target_path))
            if cleanup_tmp:
                try:
                    tmp_cif.unlink()  # type: ignore[name-defined]
                except Exception:
                    pass
        elif export_format == "sdf":
            if Chem is None:
                raise RuntimeError("RDKit required for SDF export")
            # Filter to ligands only when exporting SDF
            df_lig = df[df['group'].str.upper() == 'HETATM'] if 'group' in df.columns else df
            if df_lig.empty:
                raise ValueError("No HETATM atoms available for SDF export")

            # Apply optional filters
            df_lig = df_lig.reset_index()
            if include_chains and 'auth_chain_id' in df_lig.columns:
                df_lig = df_lig[df_lig['auth_chain_id'].isin(include_chains)]
            if include_comp_ids and 'auth_comp_id' in df_lig.columns:
                df_lig = df_lig[df_lig['auth_comp_id'].isin(include_comp_ids)]
            if include_res_ids and 'res_id' in df_lig.columns:
                df_lig = df_lig[df_lig['res_id'].isin(include_res_ids)]
            if exclude_res_ids and 'res_id' in df_lig.columns:
                df_lig = df_lig[~df_lig['res_id'].isin(exclude_res_ids)]

            # Resolve grouping columns
            if ligand_group_by is None:
                if 'res_id' in df_lig.columns:
                    group_cols = ['res_id']
                    key_is_resid = True
                else:
                    cols = [c for c in ['auth_chain_id', 'auth_seq_id', 'auth_comp_id', 'insertion'] if c in df_lig.columns]
                    group_cols = cols if cols else None
                    key_is_resid = False
            else:
                group_cols = [ligand_group_by] if isinstance(ligand_group_by, str) else list(ligand_group_by)
                key_is_resid = ('res_id' in group_cols)

            group_iter = [(name, df_lig)] if not group_cols else df_lig.groupby(group_cols, dropna=False)

            writer = Chem.SDWriter(str(target_path))
            try:
                wrote = 0
                for key, sub in group_iter:
                    if len(sub) < min_atoms:
                        continue
                    mol = canonical_dataframe_to_mol(sub)
                    # Set a meaningful name
                    if key_is_resid:
                        mol_name = f"{name}_{key}"
                    else:
                        vals = key if isinstance(key, tuple) else (key,)
                        mol_name = f"{name}_" + "_".join(str(v) for v in vals if v is not None and str(v) != 'nan')
                    if not mol.HasProp('Name'):
                        mol.SetProp('Name', mol_name)
                    writer.write(mol)
                    wrote += 1
                if wrote == 0:
                    raise ValueError("No ligand groups found for SDF export")
            finally:
                writer.close()
        else:
            raise ValueError("Unsupported export format; use 'cif' or 'sdf'")

        self.logger.info("Exported structure '%s' to %s", name, target_path)
        return target_path

    def export_dataset(
        self,
        dataset_name: str,
        output_dir: Optional[Path] = None,
        format: Optional[str] = None,
        overwrite: bool = False,
        name_pattern: Optional[str] = None,
        # SDF ligand controls
        ligand_group_by: Optional[Union[str, list[str]]] = None,
        include_chains: Optional[list[str]] = None,
        include_comp_ids: Optional[list[str]] = None,
        include_res_ids: Optional[list[str]] = None,
        exclude_res_ids: Optional[list[str]] = None,
        min_atoms: int = 1,
    ) -> Dict[str, Union[str, Dict[str, str]]]:
        export_format = (format or "cif").lower()

        entities = self.processor.get_dataset_entities(dataset_name)
        if not entities:
            raise ValueError(f"Dataset '{dataset_name}' has no entities")
        output_dir = Path(output_dir) if output_dir is not None else (get_structure_dataset_path(
            self.processor.paths,
            dataset_name,
            extension=export_format,
        ).parent)
        output_dir.mkdir(parents=True, exist_ok=True)

        # CIF: write one file per entity
        if export_format == "cif":
            outputs: Dict[str, Path] = {}
            for name in entities:
                target = get_structure_entity_path(self.processor.paths, name, extension='cif') if not name_pattern else (output_dir / name_pattern.format(name=name, ext="cif"))
                self.export_entity(name, target, format="cif", overwrite=overwrite)
                outputs[name] = target
            return {
                "files": {k: str(v) for k, v in outputs.items()},
                "format": "cif",
            }

        # SDF: write one multi-molecule file per dataset
        if export_format == "sdf":
            if Chem is None:
                raise RuntimeError("RDKit required for SDF export")
            target = get_structure_dataset_path(self.processor.paths, dataset_name, extension='sdf') if not name_pattern else (output_dir / name_pattern.format(name=dataset_name, ext="sdf"))
            if target.exists() and not overwrite:
                raise FileExistsError(f"Dataset export already exists at {target}")
            writer = Chem.SDWriter(str(target))
            try:
                wrote = 0
                for ent_name in entities:
                    df = self._get_structure(ent_name)
                    df_lig = df[df['group'].str.upper() == 'HETATM'] if 'group' in df.columns else df
                    if df_lig.empty:
                        continue
                    df_lig = df_lig.reset_index()
                    if include_chains and 'auth_chain_id' in df_lig.columns:
                        df_lig = df_lig[df_lig['auth_chain_id'].isin(include_chains)]
                    if include_comp_ids and 'auth_comp_id' in df_lig.columns:
                        df_lig = df_lig[df_lig['auth_comp_id'].isin(include_comp_ids)]
                    if include_res_ids and 'res_id' in df_lig.columns:
                        df_lig = df_lig[df_lig['res_id'].isin(include_res_ids)]
                    if exclude_res_ids and 'res_id' in df_lig.columns:
                        df_lig = df_lig[~df_lig['res_id'].isin(exclude_res_ids)]

                    # Determine grouping
                    if ligand_group_by is None:
                        if 'res_id' in df_lig.columns:
                            group_cols = ['res_id']
                            key_is_resid = True
                        else:
                            cols = [c for c in ['auth_chain_id', 'auth_seq_id', 'auth_comp_id', 'insertion'] if c in df_lig.columns]
                            group_cols = cols if cols else None
                            key_is_resid = False
                    else:
                        group_cols = [ligand_group_by] if isinstance(ligand_group_by, str) else list(ligand_group_by)
                        key_is_resid = ('res_id' in group_cols)

                    group_iter = [(ent_name, df_lig)] if not group_cols else df_lig.groupby(group_cols, dropna=False)

                    for key, sub in group_iter:
                        if sub.empty:
                            continue
                        if len(sub) < min_atoms:
                            continue
                        mol = canonical_dataframe_to_mol(sub)
                        if key_is_resid:
                            mol_name = f"{ent_name}_{key}"
                        else:
                            vals = key if isinstance(key, tuple) else (key,)
                            mol_name = f"{ent_name}_" + "_".join(str(v) for v in vals if v is not None and str(v) != 'nan')
                        if not mol.HasProp('Name'):
                            mol.SetProp('Name', mol_name)
                        writer.write(mol)
                        wrote += 1
                if wrote == 0:
                    raise ValueError("No ligand groups found for SDF dataset export")
            finally:
                writer.close()
            return {
                "dataset_file": str(target),
                "artifact_path": to_data_relative_path(self.processor.paths, target),
                "format": "sdf",
            }

        raise ValueError("Unsupported export format; use 'cif', 'pdb', or 'sdf'")

    def _get_structure(self, name: str) -> pd.DataFrame:
        df = self.processor.frames.get(name)
        if df is None:
            df = self.processor.load_entity(name)
        if df is None:
            raise ValueError(f"Structure '{name}' not found")
        if not isinstance(df, pd.DataFrame):
            raise TypeError("StructureProcessor must return DataFrame instances")
        return df.reset_index()

    def _filter_by_chains(self, df: pd.DataFrame, chains: list[str]) -> pd.DataFrame:
        chains_set = set(chains)
        mask = df['auth_chain_id'].isin(chains_set)
        if 'label_chain_id' in df.columns:
            mask |= df['label_chain_id'].isin(chains_set)
        filtered = df.loc[mask]
        if filtered.empty:
            self.logger.warning("No atoms found for chains %s; exporting full structure", chains)
            return df
        return filtered
