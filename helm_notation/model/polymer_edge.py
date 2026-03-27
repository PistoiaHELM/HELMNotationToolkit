"""Data model for polymer edges in complex polymer notation."""

from __future__ import annotations

from dataclasses import dataclass, field


# Delimiters
EDGE_COMPONENT_DELIMITER = ","
CONNECTION_SEPARATOR = "-"
MONOMER_ATTACHMENT_SEPARATOR = ":"
NODE_CONCATENATOR = "+"
NODE_CONCATENATOR_REGEX = r"\+"

# Edge types
STANDARD_EDGE = 1
PAIR_EDGE = 2
GENERIC_EDGE = 3

# Attachment types
STANDARD_EDGE_ATTACHMENT = 1
PAIR_EDGE_ATTACHMENT = 2
GENERIC_EDGE_ATTACHMENT = 3

# Identifiers
PAIR_EDGE_KEY = "pair"
GENERIC_EDGE_KEY = "generic"
UNKNOWN_ATTACHMENT = "?"


@dataclass
class PolymerEdge:
    """Connection between polymer nodes in the HELM graph."""

    edge_notation: str = ""
    source_node: str = ""
    target_node: str = ""
    connection: str = ""
    edge_type: int = 0
    source_attachment_type: int = 0
    target_attachment_type: int = 0

    @property
    def source_nodes(self) -> list[str]:
        if self.source_node:
            return self.source_node.split("+")
        return []

    @property
    def target_nodes(self) -> list[str]:
        if self.target_node:
            return self.target_node.split("+")
        return []

    @property
    def source_connection(self) -> str | None:
        if not self.connection:
            return None
        tokens = self.connection.split(CONNECTION_SEPARATOR)
        return tokens[0]

    @property
    def target_connection(self) -> str | None:
        if not self.connection:
            return None
        tokens = self.connection.split(CONNECTION_SEPARATOR)
        return tokens[1] if len(tokens) == 2 else None

    @property
    def reverse_connection(self) -> str:
        return f"{self.target_connection}{CONNECTION_SEPARATOR}{self.source_connection}"

    @property
    def source_uid(self) -> str:
        return f"{self.source_node}{MONOMER_ATTACHMENT_SEPARATOR}{self.source_connection}"

    @property
    def target_uid(self) -> str:
        return f"{self.target_node}{MONOMER_ATTACHMENT_SEPARATOR}{self.target_connection}"

    def to_reverse_string(self) -> str:
        return (
            f"{self.target_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.source_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.target_connection}{CONNECTION_SEPARATOR}{self.source_connection}"
        )

    def __str__(self) -> str:
        return (
            f"{self.source_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.target_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.source_connection}{CONNECTION_SEPARATOR}{self.target_connection}"
        )

    @property
    def source_r(self) -> str | None:
        src_conn = self.source_connection
        if src_conn is None:
            return None
        if self.source_attachment_type == STANDARD_EDGE_ATTACHMENT:
            tokens = src_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            return tokens[1]
        return None

    @source_r.setter
    def source_r(self, value: str) -> None:
        src_conn = self.source_connection
        tar_conn = self.target_connection
        if src_conn is not None:
            tokens = src_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            if len(tokens) == 2:
                new_conn = f"{tokens[0]}{MONOMER_ATTACHMENT_SEPARATOR}{value}"
                self.connection = f"{new_conn}{CONNECTION_SEPARATOR}{tar_conn}"

    @property
    def target_r(self) -> str | None:
        tar_conn = self.target_connection
        if tar_conn is None:
            return None
        if self.target_attachment_type == STANDARD_EDGE_ATTACHMENT:
            tokens = tar_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            return tokens[1]
        return None

    @target_r.setter
    def target_r(self, value: str) -> None:
        src_conn = self.source_connection
        tar_conn = self.target_connection
        if tar_conn is not None:
            tokens = tar_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            if len(tokens) == 2:
                new_conn = f"{tokens[0]}{MONOMER_ATTACHMENT_SEPARATOR}{value}"
                self.connection = f"{src_conn}{CONNECTION_SEPARATOR}{new_conn}"

    @property
    def source_monomer_number(self) -> int:
        src_conn = self.source_connection
        if src_conn is None:
            return 0
        if self.source_attachment_type == STANDARD_EDGE_ATTACHMENT:
            tokens = src_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            return int(tokens[0])
        return 0

    @source_monomer_number.setter
    def source_monomer_number(self, num: int) -> None:
        src_conn = self.source_connection
        tar_conn = self.target_connection
        if src_conn is not None:
            tokens = src_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            if len(tokens) == 2:
                new_conn = f"{num}{MONOMER_ATTACHMENT_SEPARATOR}{tokens[1]}"
                self.connection = f"{new_conn}{CONNECTION_SEPARATOR}{tar_conn}"

    @property
    def target_monomer_number(self) -> int:
        tar_conn = self.target_connection
        if tar_conn is None:
            return 0
        if self.target_attachment_type == STANDARD_EDGE_ATTACHMENT:
            tokens = tar_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            return int(tokens[0])
        return 0

    @target_monomer_number.setter
    def target_monomer_number(self, num: int) -> None:
        src_conn = self.source_connection
        tar_conn = self.target_connection
        if tar_conn is not None:
            tokens = tar_conn.split(MONOMER_ATTACHMENT_SEPARATOR)
            if len(tokens) == 2:
                new_conn = f"{num}{MONOMER_ATTACHMENT_SEPARATOR}{tokens[1]}"
                self.connection = f"{src_conn}{CONNECTION_SEPARATOR}{new_conn}"

    @property
    def is_self_cycle(self) -> bool:
        return bool(self.source_node and self.source_node == self.target_node)

    @property
    def is_backbone_self_cycle(self) -> bool:
        if self.is_self_cycle and self.connection:
            return "R1" in self.connection and "R2" in self.connection
        return False

    @property
    def source_generic_descriptor(self) -> str | None:
        src_conn = self.source_connection
        if self.source_attachment_type == GENERIC_EDGE_ATTACHMENT and src_conn:
            return src_conn.split(MONOMER_ATTACHMENT_SEPARATOR)[1]
        return None

    @property
    def target_generic_descriptor(self) -> str | None:
        tar_conn = self.target_connection
        if self.target_attachment_type == GENERIC_EDGE_ATTACHMENT and tar_conn:
            return tar_conn.split(MONOMER_ATTACHMENT_SEPARATOR)[1]
        return None
