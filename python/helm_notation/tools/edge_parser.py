"""Port of org.helm.notation.tools.EdgeParser."""

from __future__ import annotations

from typing import List, Optional

from helm_notation.exceptions import NotationException
from helm_notation.model import (
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
    PolymerEdge,
    STANDARD_EDGE,
    STANDARD_EDGE_ATTACHMENT,
    get_polymer_type,
)


def parse(edge_notation: str) -> PolymerEdge:
    """Parse input edgeNotation into a PolymerEdge object.
    
    Args:
        edge_notation: Edge notation in format "SourceNode,TargetNode,SourceConnection-TargetConnection"
    
    Returns:
        PolymerEdge object with full edge type information
    
    Raises:
        NotationException: If edge notation is invalid
    """
    edge = PolymerEdge()
    if not edge_notation or len(edge_notation) == 0:
        raise NotationException("Empty polymer edge notation")

    edge_components = edge_notation.split(EDGE_COMPONENT_DELIMITER, -1)
    if len(edge_components) != 3:
        raise NotationException(
            "Invalid polymer edge notation: must contain three components"
        )

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


def canonicalize(polymer_edge: PolymerEdge, sorted_node_ids: List[str]) -> None:
    """Canonicalize input polymerEdge based on sortedNodeIDs.
    
    Modifies the polymer_edge in place.
    
    Args:
        polymer_edge: The edge to canonicalize
        sorted_node_ids: List of node IDs in sorted order
    
    Raises:
        NotationException: If node IDs are invalid
    """
    source_nodes = polymer_edge.source_nodes
    if len(source_nodes) > 1:
        new_source_node = _canonicalize_node(source_nodes, sorted_node_ids)
        polymer_edge.source_node = new_source_node

    target_nodes = polymer_edge.target_nodes
    if len(target_nodes) > 1:
        new_target_node = _canonicalize_node(target_nodes, sorted_node_ids)
        polymer_edge.target_node = new_target_node

    source_index = _get_source_node_index(polymer_edge, sorted_node_ids)
    target_index = _get_target_node_index(polymer_edge, sorted_node_ids)

    if source_index == target_index:
        if len(source_nodes) > 1 or len(target_nodes) > 1:
            source_index = _get_source_node_min_index(polymer_edge, sorted_node_ids)
            target_index = _get_target_node_min_index(polymer_edge, sorted_node_ids)

    if source_index > target_index:
        old_source_node = polymer_edge.source_node
        old_target_node = polymer_edge.target_node
        polymer_edge.source_node = old_target_node
        polymer_edge.target_node = old_source_node
        polymer_edge.connection = polymer_edge.reverse_connection


def sort(
    canonical_polymer_edges: List[PolymerEdge],
    sorted_node_ids: List[str],
) -> List[PolymerEdge]:
    """Sort canonical polymer edges based on sorted nodes.
    
    Args:
        canonical_polymer_edges: List of canonicalized PolymerEdge objects
        sorted_node_ids: List of node IDs in sorted order
    
    Returns:
        Sorted list of PolymerEdge objects
    
    Raises:
        NotationException: If node IDs are invalid
    """
    sorted_edge_node_list = _get_sorted_edge_node_list(
        canonical_polymer_edges, sorted_node_ids
    )

    sorted_edge_list: List[PolymerEdge] = []
    sorted_edge_key_list: List[str] = []
    group_edge_map: dict[str, List[PolymerEdge]] = {}

    for i, id in enumerate(sorted_edge_node_list):
        for next_id in sorted_edge_node_list[i:]:
            key = id + next_id
            for edge in canonical_polymer_edges:
                if edge.source_node == id and edge.target_node == next_id:
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
        l = group_edge_map[key]
        l.sort(key=lambda e: str(e))
        sorted_edge_list.extend(l)

    return sorted_edge_list


def _check_node(node: str) -> None:
    """Validate node notation."""
    if not node or len(node) == 0:
        raise NotationException("Node is empty")

    node_ids = node.split(NODE_CONCATENATOR, -1)
    if len(node_ids) > 1:
        poly_type = get_polymer_type(node_ids[0])
        for i in range(1, len(node_ids)):
            tmp = get_polymer_type(node_ids[i])
            if tmp != poly_type:
                raise NotationException("Group polymers must be of the same type")


def _parse_connection(connection: str) -> PolymerEdge:
    """Parse connection string into edge type info."""
    pe = PolymerEdge(connection=connection)

    if not connection or len(connection) == 0:
        raise NotationException("Empty connection notation")

    attachments = connection.split(CONNECTION_SEPARATOR, -1)
    if len(attachments) != 2:
        raise NotationException(
            f"Connection must contain two parts separated with {CONNECTION_SEPARATOR}"
        )

    src_att = _parse_attachment(attachments[0])
    pe.source_attachment_type = src_att

    tgt_att = _parse_attachment(attachments[1])
    pe.target_attachment_type = tgt_att

    edge_type = 0
    if src_att == STANDARD_EDGE_ATTACHMENT and tgt_att == STANDARD_EDGE_ATTACHMENT:
        edge_type = STANDARD_EDGE
    elif src_att == PAIR_EDGE_ATTACHMENT and tgt_att == PAIR_EDGE_ATTACHMENT:
        edge_type = PAIR_EDGE
    elif src_att == GENERIC_EDGE_ATTACHMENT or tgt_att == GENERIC_EDGE_ATTACHMENT:
        edge_type = GENERIC_EDGE
    else:
        raise NotationException(f"Invalid connection notation: {connection}")
    pe.edge_type = edge_type

    return pe


def _parse_attachment(attachment: str) -> int:
    """Parse attachment string to determine attachment type."""
    if not attachment or len(attachment) == 0:
        raise NotationException("Empty attachment notation")

    parts = attachment.split(MONOMER_ATTACHMENT_SEPARATOR, -1)
    if len(parts) != 2:
        raise NotationException(
            f"Attachment must contain two parts separated with {MONOMER_ATTACHMENT_SEPARATOR}"
        )

    attachment_type = 0

    if parts[0].isdigit():
        if parts[1].startswith("R"):
            try:
                _validate_attachment_label(parts[1])
                attachment_type = STANDARD_EDGE_ATTACHMENT
            except Exception:
                pass
        elif parts[1].lower() == PAIR_EDGE_KEY:
            attachment_type = PAIR_EDGE_ATTACHMENT
    elif parts[0].lower() == GENERIC_EDGE_KEY:
        attachment_type = GENERIC_EDGE_ATTACHMENT
    else:
        raise NotationException(f"Invalid attachment notation: {attachment}")

    return attachment_type


def _validate_attachment_label(label: str) -> None:
    """Validate attachment label format."""
    if not label.startswith("R"):
        raise NotationException(f"Invalid attachment label: {label}")
    num = label[1:]
    if not num.isdigit():
        raise NotationException(f"Invalid attachment label: {label}")


def _canonicalize_node(node_ids: List[str], sorted_node_ids: List[str]) -> str:
    """Canonicalize a list of node IDs based on sorted order."""
    results = [id for id in sorted_node_ids if id in node_ids]

    if len(node_ids) != len(results):
        raise NotationException("Invalid polymer node ID found in Generic connection")

    return NODE_CONCATENATOR.join(results)


def _get_target_node_index(polymer_edge: PolymerEdge, sorted_node_ids: List[str]) -> int:
    """Get the target node index."""
    nodes = polymer_edge.target_nodes
    return _get_node_index(nodes, sorted_node_ids)


def _get_source_node_index(polymer_edge: PolymerEdge, sorted_node_ids: List[str]) -> int:
    """Get the source node index."""
    nodes = polymer_edge.source_nodes
    return _get_node_index(nodes, sorted_node_ids)


def _get_node_index(nodes: List[str], sorted_node_ids: List[str]) -> int:
    """Calculate node index from list of nodes."""
    index = 0

    for node in nodes:
        found = False
        for i, tmp in enumerate(sorted_node_ids):
            if tmp == node:
                index = index + i
                found = True
                break
        if not found:
            raise NotationException("Invalid polymer node ID found in Generic connection")

    return index


def _get_target_node_min_index(polymer_edge: PolymerEdge, sorted_node_ids: List[str]) -> int:
    """Get the minimum target node index."""
    nodes = polymer_edge.target_nodes
    return _get_min_node_index(nodes, sorted_node_ids)


def _get_source_node_min_index(polymer_edge: PolymerEdge, sorted_node_ids: List[str]) -> int:
    """Get the minimum source node index."""
    nodes = polymer_edge.source_nodes
    return _get_min_node_index(nodes, sorted_node_ids)


def _get_min_node_index(nodes: List[str], sorted_node_ids: List[str]) -> int:
    """Get minimum index for a list of nodes."""
    index = 1000

    for node in nodes:
        found = False
        for i, tmp in enumerate(sorted_node_ids):
            if tmp == node:
                if i < index:
                    index = i
                found = True
                break
        if not found:
            raise NotationException("Invalid polymer node ID found in Generic connection")

    return index


def _get_sorted_edge_node_list(
    polymer_edges: List[PolymerEdge],
    sorted_node_ids: List[str],
) -> List[str]:
    """Get sorted list of edge nodes."""
    edge_node_list: List[str] = []

    index_nodes_map: dict[int, List[str]] = {}
    for pe in polymer_edges:
        src_node = pe.source_node
        src_index = _get_source_node_index(pe, sorted_node_ids)
        if src_index in index_nodes_map:
            nodes = index_nodes_map[src_index]
            if src_node not in nodes:
                nodes.append(src_node)
        else:
            index_nodes_map[src_index] = [src_node]

        tar_node = pe.target_node
        tar_index = _get_target_node_index(pe, sorted_node_ids)
        if tar_index in index_nodes_map:
            nodes = index_nodes_map[tar_index]
            if tar_node not in nodes:
                nodes.append(tar_node)
        else:
            index_nodes_map[tar_index] = [tar_node]

    for index in sorted(index_nodes_map.keys()):
        nodes = index_nodes_map[index]
        if len(nodes) > 1:
            min_index_node_map: dict[int, str] = {}
            for node in nodes:
                min_index = _get_node_min_index(node, sorted_node_ids)
                min_index_node_map[min_index] = node

            for min_index in sorted(min_index_node_map.keys()):
                node = min_index_node_map[min_index]
                edge_node_list.append(node)
        else:
            edge_node_list.extend(nodes)

    return edge_node_list


def _get_node_min_index(node: str, sorted_node_ids: List[str]) -> int:
    """Get minimum index for a single node (splits by + if needed)."""
    nodes = node.split(NODE_CONCATENATOR, -1)
    return _get_min_node_index(nodes, sorted_node_ids)