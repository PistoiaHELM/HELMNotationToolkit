"""Phase 4 (incremental): port of org.helm.notation.peptide.PeptideStructureParser.

Current scope:
- Convert peptide structure (SMILES/MolBlock/RDKit Mol) to HELM notation.
- Supports linear peptides, cyclic backbone closure, and side-chain crosslinks
    represented as standard HELM edges.

Still out-of-scope for this phase:
- generic/unknown residue fallback inline notation.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem

from helm_notation.exceptions import StructureException
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.chem import ChemEngine


@dataclass
class _FragInfo:
    frag: Chem.Mol
    monomer_id: str
    connectors: List[Tuple[int, str]]


class PeptideStructureParser:
    _instance: Optional["PeptideStructureParser"] = None

    @classmethod
    def get_instance(cls) -> "PeptideStructureParser":
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def molfile2notation(self, molfile: str) -> str:
        mol = ChemEngine.get_molecule(molfile)
        return self.molecule2notation(mol)

    def smiles2notation(self, smiles: str) -> str:
        mol = ChemEngine.get_molecule(smiles)
        return self.molecule2notation(mol)

    def molecule2notation(self, molecule: Chem.Mol) -> str:
        tmp = Chem.Mol(molecule)
        frags = Chem.GetMolFrags(tmp, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            raise StructureException("Input structure can't be a mixture")

        cut_bonds = _find_backbone_amide_bonds(tmp) + _find_disulfide_bonds(tmp)
        if not cut_bonds:
            # Single residue edge case: try direct monomer match
            mid = _match_peptide_monomer_id(tmp)
            if mid is None:
                raise StructureException("Unable to map structure to peptide monomer")
            return mid

        # Label each cut pair with a deterministic id (1..n) on both sides
        labels = [(idx + 1, idx + 1) for idx in range(len(cut_bonds))]
        broken = Chem.FragmentOnBonds(tmp, cut_bonds, addDummies=True, dummyLabels=labels)
        pieces = Chem.GetMolFrags(broken, asMols=True, sanitizeFrags=True)

        infos: List[_FragInfo] = []
        for piece in pieces:
            monomer_id = _match_peptide_monomer_id(piece)
            if monomer_id is None:
                raise StructureException("Unable to map a peptide fragment to monomer database")
            connectors = _classify_connectors(piece)
            infos.append(
                _FragInfo(
                    frag=piece,
                    monomer_id=monomer_id,
                    connectors=connectors,
                )
            )

        return _assemble_helm_from_fragments(infos)


def _find_backbone_amide_bonds(mol: Chem.Mol) -> List[int]:
    # N-C(=O) peptide bond pattern
    pattern = Chem.MolFromSmarts("[NX3]-[CX3](=O)")
    matches = mol.GetSubstructMatches(pattern)
    bond_ids = set()
    for n_idx, c_idx, _ in matches:
        bond = mol.GetBondBetweenAtoms(n_idx, c_idx)
        if bond is not None:
            bond_ids.add(bond.GetIdx())
    return sorted(bond_ids)


def _find_disulfide_bonds(mol: Chem.Mol) -> List[int]:
    pattern = Chem.MolFromSmarts("[S]-[S]")
    matches = mol.GetSubstructMatches(pattern)
    bond_ids = set()
    for s1, s2 in matches:
        bond = mol.GetBondBetweenAtoms(s1, s2)
        if bond is not None:
            bond_ids.add(bond.GetIdx())
    return sorted(bond_ids)


def _classify_connectors(frag: Chem.Mol) -> List[Tuple[int, str]]:
    connectors: List[Tuple[int, str]] = []

    for atom in frag.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso <= 0:
            continue
        if atom.GetDegree() != 1:
            continue
        nbr = atom.GetNeighbors()[0]
        sym = nbr.GetSymbol()
        if sym == "N":
            connectors.append((iso, "R1"))
        elif sym == "C":
            connectors.append((iso, "R2"))
        elif sym == "S":
            connectors.append((iso, "R3"))

    connectors.sort(key=lambda x: (x[0], x[1]))
    return connectors


def _assemble_helm_from_fragments(infos: List[_FragInfo]) -> str:
    if not infos:
        raise StructureException("No fragments found for peptide chain")

    conn_map: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
    for frag_idx, info in enumerate(infos):
        for conn_id, rlabel in info.connectors:
            conn_map[conn_id].append((frag_idx, rlabel))

    backbone_directed: List[Tuple[int, int, int]] = []
    side_edges_raw: List[Tuple[Tuple[int, str], Tuple[int, str], int]] = []
    for conn_id, occ in conn_map.items():
        if len(occ) != 2:
            continue
        (f1, r1), (f2, r2) = occ
        labels = {r1, r2}
        if labels == {"R1", "R2"}:
            if r1 == "R2":
                backbone_directed.append((f1, f2, conn_id))
            else:
                backbone_directed.append((f2, f1, conn_id))
        else:
            side_edges_raw.append(((f1, r1), (f2, r2), conn_id))

    chains, chain_is_cycle = _build_backbone_chains(len(infos), backbone_directed)

    frag_to_chain_pos: Dict[int, Tuple[int, int]] = {}
    for chain_idx, chain in enumerate(chains):
        for pos, frag_idx in enumerate(chain, start=1):
            frag_to_chain_pos[frag_idx] = (chain_idx, pos)

    chain_notations = [".".join(infos[f].monomer_id for f in chain) for chain in chains]

    cycle_edges: List[Tuple[int, int, int, str, int, str]] = []
    for chain_idx, chain in enumerate(chains):
        if not chain_is_cycle[chain_idx]:
            continue
        last_pos = len(chain)
        cycle_edges.append((chain_idx, last_pos, chain_idx, "R2", 1, "R1"))

    side_edges: List[Tuple[int, int, str, int, int, str]] = []
    for (f1, r1), (f2, r2), _conn_id in side_edges_raw:
        if f1 not in frag_to_chain_pos or f2 not in frag_to_chain_pos:
            continue
        c1, p1 = frag_to_chain_pos[f1]
        c2, p2 = frag_to_chain_pos[f2]
        side_edges.append((c1, p1, r1, c2, p2, r2))

    # Simple notation return only for one linear chain and no extra edges.
    if len(chains) == 1 and not chain_is_cycle[0] and not side_edges:
        return chain_notations[0]

    node_strings = [f"PEPTIDE{i+1}{{{notation}}}" for i, notation in enumerate(chain_notations)]

    edge_strings: List[str] = []
    for cidx, src_pos, cidx2, src_r, tgt_pos, tgt_r in cycle_edges:
        edge_strings.append(
            f"PEPTIDE{cidx+1},PEPTIDE{cidx2+1},{src_pos}:{src_r}-{tgt_pos}:{tgt_r}"
        )

    for c1, p1, r1, c2, p2, r2 in side_edges:
        left = (c1, p1, r1)
        right = (c2, p2, r2)
        if right < left:
            (c1, p1, r1), (c2, p2, r2) = (c2, p2, r2), (c1, p1, r1)
        edge_strings.append(f"PEPTIDE{c1+1},PEPTIDE{c2+1},{p1}:{r1}-{p2}:{r2}")

    edge_strings = sorted(set(edge_strings))

    complex_notation = "|".join(node_strings) + "$" + "|".join(edge_strings) + "$$$"

    # Canonicalize ordering using complex parser logic (already implemented in Phase 3).
    from helm_notation.complex_notation_parser import ComplexNotationParser

    return ComplexNotationParser.get_canonical_notation(complex_notation)


def _build_backbone_chains(
    n_frags: int,
    directed_edges: List[Tuple[int, int, int]],
) -> Tuple[List[List[int]], List[bool]]:
    out_adj: Dict[int, List[int]] = defaultdict(list)
    in_adj: Dict[int, List[int]] = defaultdict(list)
    undirected: Dict[int, Set[int]] = defaultdict(set)

    for src, tgt, _cid in directed_edges:
        out_adj[src].append(tgt)
        in_adj[tgt].append(src)
        undirected[src].add(tgt)
        undirected[tgt].add(src)

    for i in range(n_frags):
        out_adj.setdefault(i, [])
        in_adj.setdefault(i, [])
        undirected.setdefault(i, set())

    # Connected components in backbone graph (isolated residues become solo chains).
    visited: Set[int] = set()
    components: List[List[int]] = []
    for i in range(n_frags):
        if i in visited:
            continue
        q = deque([i])
        visited.add(i)
        comp: List[int] = []
        while q:
            cur = q.popleft()
            comp.append(cur)
            for nb in undirected[cur]:
                if nb not in visited:
                    visited.add(nb)
                    q.append(nb)
        components.append(sorted(comp))

    chains: List[List[int]] = []
    is_cycle: List[bool] = []

    for comp in components:
        starts = [node for node in comp if len(in_adj[node]) == 0 and len(out_adj[node]) <= 1]
        if len(comp) == 1 and len(out_adj[comp[0]]) == 0 and len(in_adj[comp[0]]) == 0:
            chains.append([comp[0]])
            is_cycle.append(False)
            continue

        if len(starts) == 1:
            # Linear chain
            order: List[int] = []
            cur = starts[0]
            seen_local: Set[int] = set()
            while True:
                if cur in seen_local:
                    raise StructureException("Cycle detected in linear peptide ordering")
                seen_local.add(cur)
                order.append(cur)
                outs = out_adj[cur]
                if len(outs) == 0:
                    break
                if len(outs) != 1:
                    raise StructureException("Branching peptide backbone is not supported")
                cur = outs[0]

            if set(order) != set(comp):
                raise StructureException("Disconnected peptide backbone component detected")

            chains.append(order)
            is_cycle.append(False)
            continue

        if len(starts) == 0:
            # Candidate cycle: each node should have in=1 and out=1
            if any(len(in_adj[n]) != 1 or len(out_adj[n]) != 1 for n in comp):
                raise StructureException("Unsupported peptide topology")
            start = min(comp)
            order: List[int] = []
            cur = start
            seen_local: Set[int] = set()
            while cur not in seen_local:
                seen_local.add(cur)
                order.append(cur)
                cur = out_adj[cur][0]

            if cur != start or set(order) != set(comp):
                raise StructureException("Invalid cyclic peptide backbone topology")

            chains.append(order)
            is_cycle.append(True)
            continue

        raise StructureException("Branching peptide backbone is not supported")

    return chains, is_cycle


@lru_cache(maxsize=1)
def _peptide_signature_index() -> Dict[str, List[str]]:
    factory = MonomerFactory.get_instance()
    peptide_db = factory.get_monomer_db().get("PEPTIDE", {})

    idx: Dict[str, List[str]] = {}
    for monomer_id, monomer in peptide_db.items():
        smi = monomer.can_smiles
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        sig = _normalized_signature(mol)
        idx.setdefault(sig, []).append(monomer_id)

    # deterministic ordering for ambiguous signatures
    for sig in idx:
        idx[sig] = sorted(idx[sig], key=lambda x: (len(x), x))

    return idx


def _normalized_signature(mol: Chem.Mol) -> str:
    clone = Chem.Mol(mol)
    for atom in clone.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        atom.SetIsotope(0)
        atom.SetAtomMapNum(0)
        for key in list(atom.GetPropsAsDict().keys()):
            if atom.HasProp(key):
                atom.ClearProp(key)
    return Chem.MolToSmiles(clone, canonical=True, isomericSmiles=True)


def _match_peptide_monomer_id(frag: Chem.Mol) -> Optional[str]:
    sig = _normalized_signature(frag)
    hits = _peptide_signature_index().get(sig)
    if not hits:
        return None
    return hits[0]
