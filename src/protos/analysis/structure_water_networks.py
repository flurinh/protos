"""Water network analysis utilities built on StructureProcessor outputs.

This module treats water molecules as ligand-like entities and identifies
networks of alternating residue and water nodes that appear in crystal
structures. Results are returned as JSON-friendly dictionaries so other
components (including MCP tooling) can record summaries or render the
networks without touching processor internals.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Any, Deque, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


WATER_CODES = {"HOH", "WAT"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class ResidueKey:
    """Identifier for a residue node."""

    chain_id: str
    seq_id: str
    res_name: str

    @property
    def label(self) -> str:
        return f"{self.chain_id}:{self.res_name}{self.seq_id}"


@dataclass(frozen=True)
class WaterKey:
    """Identifier for a water node."""

    chain_id: str
    seq_id: str

    @property
    def label(self) -> str:
        return f"{self.chain_id}:HOH{self.seq_id}"


def _format_seq_id(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "?"
    try:
        if float(value).is_integer():
            return str(int(float(value)))
    except Exception:
        pass
    return str(value)


class WaterNetworkAnalyzer:
    """Detects water-mediated residue interaction networks."""

    def __init__(self, structure_frame: pd.DataFrame) -> None:
        self.frame = structure_frame.reset_index()
        self.protein_atoms = self.frame[self.frame["group"] == "ATOM"].copy()
        self.water_atoms = self.frame[
            (self.frame["group"] == "HETATM")
            & (self.frame["res_name3l"].isin(WATER_CODES))
        ].copy()

        self._protein_coords = (
            self.protein_atoms[["x", "y", "z"]].values
            if not self.protein_atoms.empty
            else np.empty((0, 3))
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def analyze(
        self,
        *,
        residue_cutoff: float = 3.4,
        water_water_cutoff: float = 3.4,
        hydrogen_bond_cutoff: float = 3.2,
    ) -> Dict[str, Any]:
        """Compute water networks for the supplied structure frame."""

        if self.water_atoms.empty or self.protein_atoms.empty:
            return {
                "networks": [],
                "summary": {
                    "network_count": 0,
                    "residue_count": 0,
                    "water_count": int(self.water_atoms.groupby(
                        ["auth_chain_id", "auth_seq_id"]
                    ).ngroups),
                    "bridging_water_count": 0,
                    "max_residue_path_length": 0,
                    "residues_with_grn": 0,
                },
            }

        water_nodes = self._collect_water_nodes()
        residue_metadata = self._collect_residue_metadata()

        (
            residue_water_edges,
            adjacency,
            involved_nodes,
        ) = self._compute_residue_contacts(
            water_nodes,
            residue_metadata,
            residue_cutoff=residue_cutoff,
            hydrogen_bond_cutoff=hydrogen_bond_cutoff,
        )

        water_water_edges = self._compute_water_water_edges(
            water_nodes, adjacency, cutoff=water_water_cutoff
        )

        networks = self._build_networks(
            involved_nodes,
            adjacency,
            residue_water_edges,
            water_water_edges,
            residue_metadata,
            water_nodes,
        )
        for index, network in enumerate(networks, start=1):
            network["network_id"] = index

        total_residues = sum(n["summary"]["residue_count"] for n in networks)
        total_waters = sum(n["summary"]["water_count"] for n in networks)
        max_path = max((n["summary"].get("max_residue_path_length", 0) for n in networks), default=0)
        residues_with_grn = sum(
            len(
                [r for r in net["residues"] if r.get("grn_labels")]
            )
            for net in networks
        )
        bridging_total = sum(len(net.get("bridging_waters", [])) for net in networks)

        return {
            "networks": networks,
            "summary": {
                "network_count": len(networks),
                "residue_count": total_residues,
                "water_count": total_waters,
                "bridging_water_count": bridging_total,
                "max_residue_path_length": max_path,
                "residues_with_grn": residues_with_grn,
            },
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _collect_water_nodes(self) -> Dict[WaterKey, Dict[str, Any]]:
        nodes: Dict[WaterKey, Dict[str, Any]] = {}
        for (chain_id, seq_id), atoms in self.water_atoms.groupby(
            ["auth_chain_id", "auth_seq_id"]
        ):
            seq_label = _format_seq_id(seq_id)
            key = WaterKey(chain_id=chain_id or "?", seq_id=seq_label)
            oxygen = atoms[atoms["atom_name"].str.upper().str.startswith("O")]
            if not oxygen.empty:
                coords = oxygen[["x", "y", "z"]].values.mean(axis=0)
            else:
                coords = atoms[["x", "y", "z"]].values.mean(axis=0)

            nodes[key] = {
                "type": "water",
                "key": key,
                "label": key.label,
                "chain_id": key.chain_id,
                "seq_id": key.seq_id,
                "coords": coords,
                "atom_ids": atoms["atom_id"].astype(int).tolist(),
                "res_name": atoms["res_name3l"].iloc[0] if not atoms.empty else "HOH",
                "atoms": atoms,
            }
        return nodes

    def _collect_residue_metadata(self) -> Dict[ResidueKey, Dict[str, Any]]:
        metadata: Dict[ResidueKey, Dict[str, Any]] = {}
        if self.protein_atoms.empty:
            return metadata

        for (chain_id, seq_id, res_name), atoms in self.protein_atoms.groupby(
            ["auth_chain_id", "auth_seq_id", "res_name3l"]
        ):
            seq_label = _format_seq_id(seq_id)
            key = ResidueKey(
                chain_id=chain_id or "?",
                seq_id=seq_label,
                res_name=res_name or "UNK",
            )
            grn_values = [
                str(value)
                for value in atoms["grn"].dropna().astype(str).unique()
                if str(value).strip()
            ]
            metadata[key] = {
                "type": "residue",
                "key": key,
                "label": key.label,
                "chain_id": key.chain_id,
                "seq_id": key.seq_id,
                "res_name": key.res_name,
                "grn_labels": grn_values,
                "atoms": atoms,
            }
        return metadata

    def _compute_residue_contacts(
        self,
        water_nodes: Dict[WaterKey, Dict[str, Any]],
        residue_metadata: Dict[ResidueKey, Dict[str, Any]],
        *,
        residue_cutoff: float,
        hydrogen_bond_cutoff: float,
    ) -> Tuple[
        List[Dict[str, Any]],
        Dict[str, Set[str]],
        Set[str],
    ]:
        adjacency: Dict[str, Set[str]] = defaultdict(set)
        edges: List[Dict[str, Any]] = []
        involved_nodes: Set[str] = set()

        if not water_nodes or not residue_metadata or self._protein_coords.size == 0:
            return edges, adjacency, involved_nodes

        # Pre-compute lookup from dataframe index to residue key
        residue_lookup: Dict[int, ResidueKey] = {}
        for key, meta in residue_metadata.items():
            atom_indices = meta["atoms"].index.tolist()
            for idx in atom_indices:
                residue_lookup[idx] = key

        protein_coords = self._protein_coords

        for water_key, water_info in water_nodes.items():
            atoms = water_info["atoms"]
            water_coords = atoms[["x", "y", "z"]].values
            if water_coords.size == 0:
                continue

            distances = (
                np.linalg.norm(
                    protein_coords[None, :, :] - water_coords[:, None, :], axis=-1
                )
                if protein_coords.size
                else np.empty((len(water_coords), 0))
            )

            within_cutoff = np.where(distances <= residue_cutoff)
            if within_cutoff[1].size == 0:
                continue

            residue_contacts: Dict[ResidueKey, Dict[str, Any]] = {}

            for water_idx, prot_idx in zip(*within_cutoff):
                atom_row = self.protein_atoms.iloc[prot_idx]
                residue_key = residue_lookup.get(atom_row.name)
                if residue_key is None:
                    chain = atom_row.get("auth_chain_id") or "?"
                    seq_id = _format_seq_id(atom_row.get("auth_seq_id"))
                    res_name = atom_row.get("res_name3l") or "UNK"
                    residue_key = ResidueKey(chain, seq_id, res_name)
                    if residue_key not in residue_metadata:
                        residue_metadata[residue_key] = {
                            "type": "residue",
                            "key": residue_key,
                            "label": residue_key.label,
                            "chain_id": chain,
                            "seq_id": seq_id,
                            "res_name": res_name,
                            "grn_labels": [],
                            "atoms": self.protein_atoms.iloc[[prot_idx]],
                        }

                min_distance = float(distances[water_idx, prot_idx])
                contact_entry = residue_contacts.setdefault(
                    residue_key,
                    {
                        "residue": residue_metadata[residue_key],
                        "distance": min_distance,
                        "contact_atoms": set(),
                        "hydrogen_bond": False,
                        "backbone": False,
                    },
                )

                contact_entry["distance"] = min(
                    contact_entry["distance"],
                    min_distance,
                )
                atom_name = str(atom_row.get("atom_name") or "")
                contact_entry["contact_atoms"].add(atom_name)

                if atom_name in BACKBONE_ATOMS:
                    contact_entry["backbone"] = True

                if atom_name.startswith("O") or atom_name.startswith("N"):
                    if min_distance <= hydrogen_bond_cutoff:
                        contact_entry["hydrogen_bond"] = True

            for residue_key, info in residue_contacts.items():
                edge_payload = {
                    "type": "residue_water",
                    "residue_key": residue_key.label,
                    "water_key": water_key.label,
                    "distance": float(info["distance"]),
                    "contact_atoms": sorted(info["contact_atoms"]),
                    "hydrogen_bond": info["hydrogen_bond"],
                    "backbone_contact": info["backbone"],
                    "residue": self._serialize_residue(residue_metadata[residue_key]),
                    "water": self._serialize_water(water_info),
                }
                edges.append(edge_payload)

                adjacency[residue_key.label].add(water_key.label)
                adjacency[water_key.label].add(residue_key.label)
                involved_nodes.add(residue_key.label)
                involved_nodes.add(water_key.label)

        return edges, adjacency, involved_nodes

    def _compute_water_water_edges(
        self,
        water_nodes: Dict[WaterKey, Dict[str, Any]],
        adjacency: Dict[str, Set[str]],
        *,
        cutoff: float,
    ) -> List[Dict[str, Any]]:
        edges: List[Dict[str, Any]] = []
        water_items = list(water_nodes.values())
        for idx, node_a in enumerate(water_items):
            coords_a = node_a["coords"]
            for node_b in water_items[idx + 1 :]:
                coords_b = node_b["coords"]
                distance = float(np.linalg.norm(coords_a - coords_b))
                if distance > cutoff:
                    continue
                edge_payload = {
                    "type": "water_water",
                    "water_a": node_a["label"],
                    "water_b": node_b["label"],
                    "distance": distance,
                }
                edges.append(edge_payload)
                adjacency[node_a["label"]].add(node_b["label"])
                adjacency[node_b["label"]].add(node_a["label"])
        return edges

    def _build_networks(
        self,
        involved_nodes: Set[str],
        adjacency: Dict[str, Set[str]],
        residue_water_edges: Sequence[Dict[str, Any]],
        water_water_edges: Sequence[Dict[str, Any]],
        residue_metadata: Dict[ResidueKey, Dict[str, Any]],
        water_nodes: Dict[WaterKey, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        edge_lookup: Dict[frozenset[str], Dict[str, Any]] = {}
        for edge in residue_water_edges:
            key = frozenset([edge["residue_key"], edge["water_key"]])
            edge_lookup[key] = edge
        for edge in water_water_edges:
            key = frozenset([edge["water_a"], edge["water_b"]])
            edge_lookup[key] = edge

        networks: List[Dict[str, Any]] = []
        visited: Set[str] = set()
        water_by_label = {info["label"]: info for info in water_nodes.values()}
        residue_by_label = {meta["label"]: meta for meta in residue_metadata.values()}

        for node in involved_nodes:
            if node in visited:
                continue
            component_nodes: Set[str] = set()
            queue: Deque[str] = deque([node])
            visited.add(node)

            while queue:
                current = queue.popleft()
                component_nodes.add(current)
                for neighbour in adjacency.get(current, set()):
                    if neighbour not in visited:
                        visited.add(neighbour)
                        queue.append(neighbour)

            network = self._serialize_component(
                component_nodes,
                adjacency,
                edge_lookup,
                residue_by_label,
                water_by_label,
            )
            networks.append(network)

        return networks

    def _serialize_component(
        self,
        component_nodes: Set[str],
        adjacency: Dict[str, Set[str]],
        edge_lookup: Dict[frozenset[str], Dict[str, Any]],
        residue_by_label: Dict[str, Dict[str, Any]],
        water_by_label: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        residue_nodes: List[Dict[str, Any]] = []
        water_node_payloads: List[Dict[str, Any]] = []
        node_payload_map: Dict[str, Dict[str, Any]] = {}

        total_waters = 0
        total_residues = 0

        for node in component_nodes:
            if node in water_by_label:
                total_waters += 1
                payload = self._serialize_water(water_by_label[node])
                node_payload_map[node] = payload
                water_node_payloads.append(payload)
                continue

            residue_meta = residue_by_label.get(node)
            if residue_meta is not None:
                total_residues += 1
                payload = self._serialize_residue(residue_meta)
                node_payload_map[node] = payload
                residue_nodes.append(payload)

        residue_nodes.sort(key=lambda item: (item["chain_id"], item["seq_id"]))
        water_node_payloads.sort(key=lambda item: (item["chain_id"], item["seq_id"]))

        residue_keys = [node["label"] for node in residue_nodes]
        water_keys = [node["label"] for node in water_node_payloads]

        bridging_waters = [
            node_payload_map[label]
            for label in water_keys
            if len(
                {
                    neighbour
                    for neighbour in adjacency.get(label, set())
                    if neighbour in node_payload_map
                    and node_payload_map[neighbour]["type"] == "residue"
                }
            )
            >= 2
        ]

        residue_pair_paths = self._compute_residue_paths(
            residue_keys,
            adjacency,
            edge_lookup,
            node_payload_map,
        )

        summary = {
            "residue_count": total_residues,
            "water_count": total_waters,
            "chains": sorted({node["chain_id"] for node in residue_nodes}),
            "grn_labels": sorted({grn for node in residue_nodes for grn in node["grn_labels"]}),
            "bridging_water_count": len(bridging_waters),
            "max_residue_path_length": max(
                (path["length"] for path in residue_pair_paths),
                default=0,
            ),
        }

        return {
            "residues": residue_nodes,
            "waters": water_node_payloads,
            "residue_water_edges": list(residue_water_edges_from_component(
                component_nodes, edge_lookup
            )),
            "water_water_edges": list(water_water_edges_from_component(
                component_nodes, edge_lookup
            )),
            "bridging_waters": bridging_waters,
            "residue_pair_paths": residue_pair_paths,
            "summary": summary,
        }

    def _compute_residue_paths(
        self,
        residue_labels: Sequence[str],
        adjacency: Dict[str, Set[str]],
        edge_lookup: Dict[frozenset[str], Dict[str, Any]],
        node_payload_map: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        from itertools import combinations

        paths: List[Dict[str, Any]] = []
        for source, target in combinations(residue_labels, 2):
            path = self._shortest_path(adjacency, source, target)
            if not path:
                continue
            serialized_path = [node_payload_map[n] for n in path]
            edge_details = []
            for start, end in zip(path, path[1:]):
                key = frozenset([start, end])
                edge_info = edge_lookup.get(key, {})
                edge_details.append({
                    "type": edge_info.get("type"),
                    "distance": edge_info.get("distance"),
                    "between": [start, end],
                })

            paths.append(
                {
                    "source": node_payload_map[source]["label"],
                    "target": node_payload_map[target]["label"],
                    "source_grn": node_payload_map[source].get("grn_labels", []),
                    "target_grn": node_payload_map[target].get("grn_labels", []),
                    "path_nodes": serialized_path,
                    "edges": edge_details,
                    "length": len(path) - 1,
                }
            )
        return paths

    def _shortest_path(
        self, adjacency: Dict[str, Set[str]], source: str, target: str
    ) -> Optional[List[str]]:
        queue: Deque[str] = deque([source])
        parents: Dict[str, Optional[str]] = {source: None}

        while queue:
            node = queue.popleft()
            if node == target:
                break
            for neighbour in adjacency.get(node, set()):
                if neighbour not in parents:
                    parents[neighbour] = node
                    queue.append(neighbour)

        if target not in parents:
            return None

        path: List[str] = []
        current: Optional[str] = target
        while current is not None:
            path.append(current)
            current = parents[current]
        return list(reversed(path))

    def _serialize_residue(self, meta: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "type": "residue",
            "label": meta["label"],
            "chain_id": meta["chain_id"],
            "seq_id": meta["seq_id"],
            "res_name": meta["res_name"],
            "grn_labels": meta.get("grn_labels", []),
        }

    def _serialize_water(self, meta: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "type": "water",
            "label": meta["label"],
            "chain_id": meta["chain_id"],
            "seq_id": meta["seq_id"],
            "res_name": meta.get("res_name", "HOH"),
        }


def residue_water_edges_from_component(
    component_nodes: Iterable[str],
    edge_lookup: Dict[frozenset[str], Dict[str, Any]],
) -> Iterable[Dict[str, Any]]:
    for edge_key, payload in edge_lookup.items():
        nodes = list(edge_key)
        if payload.get("type") != "residue_water":
            continue
        if all(node in component_nodes for node in nodes):
            yield payload


def water_water_edges_from_component(
    component_nodes: Iterable[str],
    edge_lookup: Dict[frozenset[str], Dict[str, Any]],
) -> Iterable[Dict[str, Any]]:
    for edge_key, payload in edge_lookup.items():
        if payload.get("type") != "water_water":
            continue
        nodes = list(edge_key)
        if all(node in component_nodes for node in nodes):
            yield payload


def analyze_water_networks(
    structure_frame: pd.DataFrame,
    *,
    residue_cutoff: float = 3.4,
    water_water_cutoff: float = 3.4,
    hydrogen_bond_cutoff: float = 3.2,
) -> Dict[str, Any]:
    """Convenience wrapper around :class:`WaterNetworkAnalyzer`."""

    analyzer = WaterNetworkAnalyzer(structure_frame)
    return analyzer.analyze(
        residue_cutoff=residue_cutoff,
        water_water_cutoff=water_water_cutoff,
        hydrogen_bond_cutoff=hydrogen_bond_cutoff,
    )


def _primary_label(node: Dict[str, Any]) -> str:
    if node.get("type") == "residue":
        labels = node.get("grn_labels") or []
        if labels:
            return labels[0]
    return str(node.get("label", "?"))


def _condense_paths(paths: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    best: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for path in paths:
        source = path.get("source")
        target = path.get("target")
        if not source or not target:
            continue
        key = tuple(sorted((source, target)))
        length = path.get("length")
        if key not in best or (length is not None and length < best[key].get("length", float("inf"))):
            best[key] = path
    return list(best.values())


def summarize_structure_water_network(structure_result: Dict[str, Any]) -> Dict[str, Any]:
    """Produce a condensed summary of a single structure's water networks."""

    summary = structure_result.get("summary", {})
    networks: List[Dict[str, Any]] = []

    for network in structure_result.get("networks", []):
        net_summary = network.get("summary", {})
        residues = sorted({
            label for label in (_primary_label(node) for node in network.get("residues", []))
            if label and label != "?"
        })

        condensed_paths: List[Dict[str, Any]] = []
        for path in _condense_paths(network.get("residue_pair_paths", [])):
            sequence_nodes = [_primary_label(node) for node in path.get("path_nodes", [])]
            condensed_paths.append(
                {
                    "source": path.get("source"),
                    "target": path.get("target"),
                    "source_grn": path.get("source_grn", []),
                    "target_grn": path.get("target_grn", []),
                    "length": path.get("length"),
                    "sequence": sequence_nodes,
                    "sequence_str": " → ".join(sequence_nodes),
                }
            )

        networks.append(
            {
                "network_id": network.get("network_id"),
                "residues": residues,
                "waters": net_summary.get("water_count", 0),
                "bridging_waters": net_summary.get("bridging_water_count", 0),
                "max_path_length": net_summary.get("max_residue_path_length", 0),
                "chains": net_summary.get("chains", []),
                "paths": condensed_paths,
            }
        )

    return {"summary": summary, "networks": networks}


def summarize_water_networks(result: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Summarize water network analysis for every structure in a result map."""

    structures = result.get("structures") if "structures" in result else result
    return {
        structure_id: summarize_structure_water_network(structure_data)
        for structure_id, structure_data in structures.items()
    }
