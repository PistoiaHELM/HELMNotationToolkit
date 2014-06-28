/*******************************************************************************
 * Copyright C 2012, The Pistoia Alliance
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/
package org.helm.notation.model;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class RNAPolymerNode extends PolymerNode {
	private String sequence; // natural sequence
	private String modifiedSequence; // modified sequence

	public RNAPolymerNode() {
	}

	public RNAPolymerNode(PolymerNode node) {
		super.setId(node.getId());
		super.setLabel(node.getLabel());
		super.setAnotation(node.getAnotation());
	}

	public String getType() {
		return Monomer.NUCLIEC_ACID_POLYMER_TYPE;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(getId() + "\t");
		sb.append(getLabel() + "\t");
		sb.append(getAnotation() + "\t");
		sb.append(getSequence() + "\t");
		sb.append(getModifiedSequence());
		return sb.toString();
	}

	/**
	 * @return the modifiedSequence
	 */
	public String getModifiedSequence() {
		return modifiedSequence;
	}

	/**
	 * @param modifiedSequence
	 *            the modifiedSequence to set
	 */
	public void setModifiedSequence(String modifiedSequence) {
		this.modifiedSequence = modifiedSequence;
	}
}
