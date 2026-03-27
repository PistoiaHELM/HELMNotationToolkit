"""Parser for polymer edge notation in HELM."""

from __future__ import annotations

from helm_notation.exceptions import NotationException
from helm_notation.model.polymer_edge import (
    CONNECTION_SEPARATOR,
    EDGE_COMPONENT_DELIMITER,
    GENERIC_EDGE,
    GENERIC_EDGE_ATTACHMENT,
    GENERIC_EDGE_KEY,
    MONOMER_ATTACHMENT_SEPARATOR,
    NODE_CONCATENATOR,
    PAIR_EDGE,
    PAIR_EDGE_ATTACHMENT,
    PAIR_EDGE_KEY,
    STANDARD_EDGE,
    STANDARD_EDGE_ATTACHMENT,
    PolymerEdge,
)
from helm_notation.model.polymer_node import get_polymer_type


def parse(edge_notation: str) -> PolymerEdge:
    """Parse edge notation string into a PolymerEdge object."""
    if not edge_notation:
        raise NotationException("Empty polymer edge notation")

    edge_components = edge_notation.split(EDGE_COMPONENT_DELIMITER)
    if len(edge_components) != 3:
        raise NotationException(
            "Invalid polymer edge notation: must contain three components"
        )

    edge = PolymerEdge()

    source_node = edge_components[0]
    _check_node(source_node)
    edge.source_node = source_node

    target_node = edge_components[1]
    _check_node(target_node)
    edge.target_node = target_node

    connection = edge_components[2]
    tmp = _parse_connection(connection)
    edge.connection = connection
    edge.edge_type = tmp.edge_type
    edge.source_attachment_type = tmp.source_attachment_type
    edge.target_attachment_type = tmp.target_attachment_type
    edge.edge_notation = edge_notation

    return edge


def canonicalize(polymer_edge: PolymerEdge, sorted_node_ids: list[str]) -> None:
    """Canonicalize a polymer edge based on sorted node IDs (modifies in place)."""
    source_nodes = polymer_edge.source_nodes
    if len(source_nodes) > 1:
        new_source = _canonicalize_node(source_nodes, sorted_node_ids)
        polymer_edge.source_node = new_source

    target_nodes = polymer_edge.target_nodes
    if len(target_nodes) > 1:
        new_target = _canonicalize_node(target_nodes, sorted_node_ids)
        polymer_edge.target_node = new_target

    source_index = _get_node_index(polymer_edge.source_nodes, sorted_node_ids)
    target_index = _get_node_index(polymer_edge.target_nodes, sorted_node_ids)

    if source_index == target_index:
        if len(source_nodes) > 1 or len(target_nodes) > 1:
            source_index = _get_min_node_index(
                polymer_edge.source_nodes, sorted_node_ids
            )
            target_index = _get_min_node_index(
                polymer_edge.target_nodes, sorted_node_ids
            )

    # Reverse if necessary
    if source_index > target_index:
        old_source = polymer_edge.source_node
        old_target = polymer_edge.target_node
        polymer_edge.source_node = old_target
        polymer_edge.target_node = old_source
        polymer_edge.connection = polymer_edge.reverse_connection


def sort(
    canonical_polymer_edges: list[PolymerEdge], sorted_node_ids: list[str]
) -> list[PolymerEdge]:
    """Sort canonical polymer edges based on sorted nodes."""
    sorted_edge_node_list = _get_sorted_edge_node_list(
        canonical_polymer_edges, sorted_node_ids
    )

    sorted_edge_list: list[PolymerEdge] = []
    sorted_edge_key_list: list[str] = []
    group_edge_map: dict[str, list[PolymerEdge]] = {}

    for i, id_ in enumerate(sorted_edge_node_list):
        for j in range(i, len(sorted_edge_node_list)):
            next_id = sorted_edge_node_list[j]
            key = id_ + next_id

            for edge in canonical_polymer_edges:
                if edge.source_node == id_ and edge.target_node == next_id:
                    tmp_edge = PolymerEdge(
                        source_node=edge.source_node,
                        target_node=edge.target_node,
                        connection=edge.connection,
                    )
                    if key in group_edge_map:
                        group_edge_map[key].append(tmp_edge)
                    else:
                        group_edge_map[key] = [tmp_edge]
                        sorted_edge_key_list.append(key)

    for key in sorted_edge_key_list:
        edges = group_edge_map[key]
        edges.sort(key=lambda e: str(e))
        sorted_edge_list.extend(edges)

    return sorted_edge_list


def validate_attachment_label(label: str) -> None:
    """Validate an attachment label (e.g., R1, R2, R3)."""
    if not label or not label.startswith("R"):
        raise NotationException(f"Invalid attachment label: {label}")
    try:
        num = int(label[1:])
        if num <= 0:
            raise NotationException(f"Invalid attachment label: {label}")
    except ValueError:
        raise NotationException(f"Invalid attachment label: {label}")


def _check_node(node: str) -> bool:
    if not node:
        raise NotationException("Node is empty")

    node_ids = node.split("+")
    if len(node_ids) > 1:
        type_ = get_polymer_type(node_ids[0])
        for nid in node_ids[1:]:
            tmp = get_polymer_type(nid)
            if tmp != type_:
                raise NotationException("Group polymers must of the same type")
    return True


def _parse_connection(connection: str) -> PolymerEdge:
    pe = PolymerEdge()
    pe.connection = connection

    if not connection:
        raise NotationException("Empty connection notation")

    attachments = connection.split(CONNECTION_SEPARATOR)
    if len(attachments) != 2:
        raise NotationException(
            f"Connection must contain two parts separated with {CONNECTION_SEPARATOR}"
        )

    src_att = _parse_attachment(attachments[0])
    pe.source_attachment_type = src_att

    tgt_att = _parse_attachment(attachments[1])
    pe.target_attachment_type = tgt_att

    if (
        src_att == STANDARD_EDGE_ATTACHMENT
        and tgt_att == STANDARD_EDGE_ATTACHMENT
    ):
        edge_type = STANDARD_EDGE
    elif (
        src_att == PAIR_EDGE_ATTACHMENT and tgt_att == PAIR_EDGE_ATTACHMENT
    ):
        edge_type = PAIR_EDGE
    elif (
        src_att == GENERIC_EDGE_ATTACHMENT
        or tgt_att == GENERIC_EDGE_ATTACHMENT
    ):
        edge_type = GENERIC_EDGE
    else:
        raise NotationException(f"Invalid connection notation: {connection}")

    pe.edge_type = edge_type
    return pe


def _parse_attachment(attachment: str) -> int:
    if not attachment:
        raise NotationException("Empty attachment notation")

    parts = attachment.split(MONOMER_ATTACHMENT_SEPARATOR)
    if len(parts) != 2:
        raise NotationException(
            f"Attachment must contain two parts separated with {MONOMER_ATTACHMENT_SEPARATOR}"
        )

    if parts[0].isdigit():
        if parts[1].startswith("R"):
            try:
                validate_attachment_label(parts[1])
                return STANDARD_EDGE_ATTACHMENT
            except Exception:
                pass
        elif parts[1].lower() == PAIR_EDGE_KEY:
            return PAIR_EDGE_ATTACHMENT
    elif parts[0].lower() == GENERIC_EDGE_KEY:
        return GENERIC_EDGE_ATTACHMENT

    raise NotationException(f"Invalid attachment notation: {attachment}")


def _canonicalize_node(
    node_ids: list[str], sorted_node_ids: list[str]
) -> str:
    results = [nid for nid in sorted_node_ids if nid in node_ids]

    if len(results) != len(node_ids):
        raise NotationException(
            "Invalid polymer node ID found in Generic connection"
        )

    return NODE_CONCATENATOR.join(results)


def _get_node_index(nodes: list[str], sorted_node_ids: list[str]) -> int:
    index = 0
    for node in nodes:
        found = False
        for i, nid in enumerate(sorted_node_ids):
            if nid == node:
                index += i
                found = True
                break
        if not found:
            raise NotationException(
                "Invalid polymer node ID found in Generic connection"
            )
    return index


def _get_min_node_index(
    nodes: list[str], sorted_node_ids: list[str]
) -> int:
    index = 1000
    for node in nodes:
        found = False
        for i, nid in enumerate(sorted_node_ids):
            if nid == node:
                if i < index:
                    index = i
                found = True
                break
        if not found:
            raise NotationException(
                "Invalid polymer node ID found in Generic connection"
            )
    return index


def _get_sorted_edge_node_list(
    polymer_edges: list[PolymerEdge], sorted_node_ids: list[str]
) -> list[str]:
    from collections import OrderedDict

    index_nodes_map: dict[int, list[str]] = {}

    for pe in polymer_edges:
        src_node = pe.source_node
        src_index = _get_node_index(pe.source_nodes, sorted_node_ids)
        if src_index in index_nodes_map:
            nodes = index_nodes_map[src_index]
            if src_node not in nodes:
                nodes.append(src_node)
        else:
            index_nodes_map[src_index] = [src_node]

        tar_node = pe.target_node
        tar_index = _get_node_index(pe.target_nodes, sorted_node_ids)
        if tar_index in index_nodes_map:
            nodes = index_nodes_map[tar_index]
            if tar_node not in nodes:
                nodes.append(tar_node)
        else:
            index_nodes_map[tar_index] = [tar_node]

    edge_node_list: list[str] = []
    for index in sorted(index_nodes_map.keys()):
        nodes = index_nodes_map[index]
        if len(nodes) > 1:
            min_index_map: dict[int, str] = {}
            for node in nodes:
                min_idx = _get_min_node_index(node.split("+"), sorted_node_ids)
                min_index_map[min_idx] = node
            for min_idx in sorted(min_index_map.keys()):
                edge_node_list.append(min_index_map[min_idx])
        else:
            edge_node_list.extend(nodes)

    return edge_node_list
