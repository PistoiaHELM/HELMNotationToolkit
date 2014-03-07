package org.helm.notation.tools;

public class SimpleNotationGroupIterator {
	protected char[] characters;
	protected String notationString;
	protected int position;
	
	public SimpleNotationGroupIterator( String notationString) {
		this.characters = notationString.toCharArray();
		this.notationString = notationString;
		this.position = 0;
	}

	public boolean hasNextGroup() {
		return this.position < this.characters.length;
	}

	
	
	public String nextGroup() {		
		int currentPosition = this.position;
		
		do {
			char currentCharacter = this.characters[currentPosition];
			if ( currentCharacter == '[') {
				currentPosition = SimpleNotationParser.getMatchingBracketPosition(this.characters, currentPosition, '[', ']');
			}
			else if ( currentCharacter == '(') {
				currentPosition = SimpleNotationParser.getMatchingBracketPosition(this.characters, currentPosition, '(', ')');
			}
			else if ( currentCharacter != '.') {
				currentPosition++;
			}
			
			
			if ( currentPosition < 0) {
				currentPosition = this.characters.length;
			}
			
			
		} while ( (currentPosition < this.characters.length)
			&& (this.characters[currentPosition] != '.'));
		
		String token = this.notationString.substring(this.position, currentPosition);
		
		this.position = currentPosition + 1;
				
		return token;
	}
}
