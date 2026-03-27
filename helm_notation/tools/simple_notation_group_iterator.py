"""Iterator class for groups within a simple notation string.

Iterates over monomer groups in a HELM simple notation, correctly
handling bracket and parenthesis delimiters.
"""

from __future__ import annotations


class SimpleNotationGroupIterator:
    """Iterator for groups within a simple notation string.

    Each call to next_group() returns the next monomer group token,
    properly handling brackets ([...]) and parentheses ((...)).
    """

    def __init__(self, notation_string: str) -> None:
        """Construct a group iterator for a notation string.

        Args:
            notation_string: The HELM simple notation string.
        """
        self.characters: list[str] = list(notation_string)
        self.notation_string: str = notation_string
        self.position: int = 0

    def __iter__(self):
        return self

    def __next__(self) -> str:
        if not self.has_next_group():
            raise StopIteration
        return self.next_group()

    def has_next_group(self) -> bool:
        """Check if there is a next group.

        Returns:
            True if the iterator has another group, otherwise False.
        """
        return self.position < len(self.characters)

    def next_group(self) -> str:
        """Return the next group.

        Returns:
            Notation group string.
        """
        current_position = self.position

        while True:
            current_character = self.characters[current_position]
            if current_character == "[":
                current_position = _get_matching_bracket_position(
                    self.characters, current_position, "[", "]"
                )
            elif current_character == "(":
                current_position = _get_matching_bracket_position(
                    self.characters, current_position, "(", ")"
                )
            elif current_character != ".":
                current_position += 1

            if current_position < 0:
                current_position = len(self.characters)

            if current_position >= len(self.characters):
                break
            if self.characters[current_position] == ".":
                break

        token = self.notation_string[self.position : current_position]
        self.position = current_position + 1
        return token


def _get_matching_bracket_position(
    characters: list[str],
    start: int,
    open_bracket: str,
    close_bracket: str,
) -> int:
    """Find the position after the matching closing bracket.

    This is a local helper that mirrors SimpleNotationParser.getMatchingBracketPosition.

    Args:
        characters: List of characters.
        start: Starting position of the open bracket.
        open_bracket: The opening bracket character.
        close_bracket: The closing bracket character.

    Returns:
        Position after the closing bracket, or -1 if not found.
    """
    depth = 0
    for i in range(start, len(characters)):
        if characters[i] == open_bracket:
            depth += 1
        elif characters[i] == close_bracket:
            depth -= 1
            if depth == 0:
                return i + 1
    return -1
