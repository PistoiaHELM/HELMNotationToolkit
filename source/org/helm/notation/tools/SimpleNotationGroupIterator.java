package org.helm.notation.tools;

/**
 * Iterator class for groups within a simple notation string
 * 
 * @author maisel
 * 
 */
public class SimpleNotationGroupIterator {

	protected char[] characters;
	protected String notationString;
	protected int position;

	/**
	 * Constructs a group iterator for a notation string
	 * 
	 * @param notationString
	 */
	public SimpleNotationGroupIterator(String notationString) {
		this.characters = notationString.toCharArray();
		this.notationString = notationString;
		this.position = 0;
	}

	/**
	 * Checks if there is a next group
	 * 
	 * @return true if the iterator has another group otherwise false
	 */
	public boolean hasNextGroup() {
		return this.position < this.characters.length;
	}

	/**
	 * Returns the next group
	 * 
	 * @return notation group
	 */
	public String nextGroup() {
		int currentPosition = this.position;

		do {
			char currentCharacter = this.characters[currentPosition];
			if (currentCharacter == '[') {
				currentPosition = SimpleNotationParser
						.getMatchingBracketPosition(this.characters,
								currentPosition, '[', ']');
			} else if (currentCharacter == '(') {
				currentPosition = SimpleNotationParser
						.getMatchingBracketPosition(this.characters,
								currentPosition, '(', ')');
			} else if (currentCharacter != '.') {
				currentPosition++;
			}

			if (currentPosition < 0) {
				currentPosition = this.characters.length;
			}

		} while ((currentPosition < this.characters.length)
				&& (this.characters[currentPosition] != '.'));

		String token = this.notationString.substring(this.position,
				currentPosition);

		this.position = currentPosition + 1;

		return token;
	}
}
