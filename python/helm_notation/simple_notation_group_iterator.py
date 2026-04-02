from __future__ import annotations


def get_matching_bracket_position(
    characters: list[str],
    position: int,
    opening_bracket: str,
    closing_bracket: str,
) -> int:
    if position < (len(characters) - 1) and characters[position] == opening_bracket:
        current_position = position
        opening_bracket_count = 1

        while opening_bracket_count > 0 and current_position < (len(characters) - 1):
            current_position += 1
            current_character = characters[current_position]
            if current_character == opening_bracket:
                opening_bracket_count += 1
            elif current_character == closing_bracket:
                opening_bracket_count -= 1

        if current_position < len(characters) and characters[current_position] == closing_bracket:
            return current_position
        return -1
    return -1


class SimpleNotationGroupIterator:
    def __init__(self, notation_string: str):
        self.characters = list(notation_string)
        self.notation_string = notation_string
        self.position = 0

    def has_next_group(self) -> bool:
        return self.position < len(self.characters)

    def next_group(self) -> str:
        current_position = self.position

        while True:
            current_character = self.characters[current_position]
            if current_character == "[":
                current_position = get_matching_bracket_position(
                    self.characters,
                    current_position,
                    "[",
                    "]",
                )
            elif current_character == "(":
                current_position = get_matching_bracket_position(
                    self.characters,
                    current_position,
                    "(",
                    ")",
                )
            elif current_character != ".":
                current_position += 1

            if current_position < 0:
                current_position = len(self.characters)

            if current_position >= len(self.characters) or self.characters[current_position] == ".":
                break

        token = self.notation_string[self.position : current_position]
        self.position = current_position + 1
        return token
