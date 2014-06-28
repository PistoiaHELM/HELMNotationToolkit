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
 * This is the data model for PolymerEdge in complex polymer notation.
 * edgeNotation is the complete edge description, sourceNode is the source
 * PolymerNode id, targetNode is the target PolymerNode id, and connection is
 * the connection property between sourceNode and targetNode. Connection could
 * take three different forms: 1. standard, as MonomerPosition:Rgroup 2:R1 2.
 * pair, as MonomerPosition:pair 5:pair 3. generic, as generic:connectionID
 * generic:C114+C288 For generic edge, source node and/or target node could be
 * groups of polymer nodes as PEPTIDE1+PEPTIDE2
 * 
 * @author zhangtianhong
 */
public class PolymerEdge {

	public static final String EDGE_COMPONENT_DELIMITER = ",";
	public static final String CONNECTION_SEPARATOR = "-";
	public static final String MONOMER_ATTACHEMENT_SEPARATOR = ":";
	public static final String NODE_CONCATENATOR = "+";
	public static final String NODE_CONCATENATOR_REGEX = "\\+";
	// Edge Types
	public static final int STANDARD_EDGE = 1;
	public static final int PAIR_EDGE = 2;
	public static final int GENERIC_EDGE = 3;
	// Attachment Types
	public static final int STANDARD_EDGE_ATTACHMENT = 1;
	public static final int PAIR_EDGE_ATTACHMENT = 2;
	public static final int GENERIC_EDGE_ATTACHMENT = 3;
	// Identifiers
	public static final String PAIR_EDGE_KEY = "pair";
	public static final String GENERIC_EDGE_KEY = "generic";
	public static final String UNKNOWN_ATTACHMENT = "?";

	private String edgeNotation; // RNA1, CHEM2, 5:R2-1:R3
	private String sourceNode; // RNA1
	private String targetNode; // RNA2
	private String connection; // 1:R2-34:R1
	private int edgeType;
	private int sourceAttachmentType;
	private int targetAttachmentType;

	public String getEdgeNotation() {
		return edgeNotation;
	}

	public void setEdgeNotation(String edgeString) {
		edgeNotation = edgeString;
	}

	public String getSourceNode() {
		return sourceNode;
	}

	public String[] getSourceNodes() {
		if (null != sourceNode) {
			return sourceNode.split(NODE_CONCATENATOR_REGEX, -1);
		}
		return new String[0];
	}

	public void setSourceNode(String sourceNode) {
		this.sourceNode = sourceNode;
	}

	public String getTargetNode() {
		return targetNode;
	}

	public String[] getTargetNodes() {
		if (null != targetNode) {
			return targetNode.split(NODE_CONCATENATOR_REGEX, -1);
		}
		return new String[0];
	}

	public void setTargetNode(String targetNode) {
		this.targetNode = targetNode;
	}

	public String getConnection() {
		return connection;
	}

	public void setConnection(String connection) {
		this.connection = connection;
	}

	public int getEdgeType() {
		return edgeType;
	}

	public void setEdgeType(int edgeType) {
		this.edgeType = edgeType;
	}

	public int getSourceAttachmentType() {
		return sourceAttachmentType;
	}

	public void setSourceAttachmentType(int sourceAttachmentType) {
		this.sourceAttachmentType = sourceAttachmentType;
	}

	public int getTargetAttachmentType() {
		return targetAttachmentType;
	}

	public void setTargetAttachmentType(int targetAttachmentType) {
		this.targetAttachmentType = targetAttachmentType;
	}

	public String getSourceConnection() {
		if (null == getConnection()) {
			return null;
		}

		String[] conTokens = getConnection().split(CONNECTION_SEPARATOR);
		return conTokens[0];
	}

	public String getTargetConnection() {
		if (null == getConnection()) {
			return null;
		}

		String[] tokens = getConnection().split(CONNECTION_SEPARATOR);
		if (tokens.length == 2) {
			return tokens[1];
		} else {
			return null;
		}
	}

	public String getReverseConnection() {
		return getTargetConnection() + CONNECTION_SEPARATOR
				+ getSourceConnection();
	}

	public String getSourceUID() {
		return getSourceNode() + MONOMER_ATTACHEMENT_SEPARATOR
				+ getSourceConnection();
	}

	public String getTargetUID() {
		return getTargetNode() + MONOMER_ATTACHEMENT_SEPARATOR
				+ getTargetConnection();
	}

	public String toReverseString() {
		return getTargetNode() + EDGE_COMPONENT_DELIMITER + getSourceNode()
				+ EDGE_COMPONENT_DELIMITER + getTargetConnection()
				+ CONNECTION_SEPARATOR + getSourceConnection();
	}

	public String toString() {
		return getSourceNode() + EDGE_COMPONENT_DELIMITER + getTargetNode()
				+ EDGE_COMPONENT_DELIMITER + getSourceConnection()
				+ CONNECTION_SEPARATOR + getTargetConnection();
	}

	public String getSourceR() {
		if (null == getSourceConnection()) {
			return null;
		}

		if (getSourceAttachmentType() == STANDARD_EDGE_ATTACHMENT) {
			String[] tokens = getSourceConnection().split(
					MONOMER_ATTACHEMENT_SEPARATOR);
			return tokens[1];
		}

		return null;
	}

	public void setSourceR(String sourceR) {
		String srcCon = getSourceConnection();
		String tarCon = getTargetConnection();

		if (null != srcCon) {
			String[] tokens = srcCon.split(MONOMER_ATTACHEMENT_SEPARATOR);
			if (tokens.length == 2) {
				String newCon = tokens[0] + MONOMER_ATTACHEMENT_SEPARATOR
						+ sourceR;
				setConnection(newCon + CONNECTION_SEPARATOR + tarCon);
			}
		}
	}

	public String getTargetR() {
		if (null == getTargetConnection()) {
			return null;
		}

		if (getTargetAttachmentType() == STANDARD_EDGE_ATTACHMENT) {
			String[] tokens = getTargetConnection().split(
					MONOMER_ATTACHEMENT_SEPARATOR);
			return tokens[1];
		}

		return null;
	}

	public void setTargetR(String targetR) {
		String srcCon = getSourceConnection();
		String tarCon = getTargetConnection();

		if (null != tarCon) {
			String[] tokens = tarCon.split(MONOMER_ATTACHEMENT_SEPARATOR);
			if (tokens.length == 2) {
				String newCon = tokens[0] + MONOMER_ATTACHEMENT_SEPARATOR
						+ targetR;
				setConnection(srcCon + CONNECTION_SEPARATOR + newCon);
			}
		}
	}

	public int getSourceMonomerNumber() {
		if (null == getSourceConnection()) {
			return 0;
		}

		if (getSourceAttachmentType() == STANDARD_EDGE_ATTACHMENT) {
			String[] tokens = getSourceConnection().split(
					MONOMER_ATTACHEMENT_SEPARATOR);
			return Integer.parseInt(tokens[0]);
		}

		return 0;
	}

	public int getTargetMonomerNumber() {
		if (null == getTargetConnection()) {
			return 0;
		}

		if (getTargetAttachmentType() == STANDARD_EDGE_ATTACHMENT) {
			String[] tokens = getTargetConnection().split(
					MONOMER_ATTACHEMENT_SEPARATOR);
			return Integer.parseInt(tokens[0]);
		}

		return 0;
	}

	public void setSourceMonomerNumber(int num) {
		String srcCon = getSourceConnection();
		String tarCon = getTargetConnection();

		if (null != srcCon) {
			String[] tokens = srcCon.split(MONOMER_ATTACHEMENT_SEPARATOR);
			if (tokens.length == 2) {
				String newCon = "" + num + MONOMER_ATTACHEMENT_SEPARATOR
						+ tokens[1];
				setConnection(newCon + CONNECTION_SEPARATOR + tarCon);
			}
		}
	}

	public void setTargetMonomerNumber(int num) {
		String srcCon = getSourceConnection();
		String tarCon = getTargetConnection();

		if (null != tarCon) {
			String[] tokens = tarCon.split(MONOMER_ATTACHEMENT_SEPARATOR);
			if (tokens.length == 2) {
				String newCon = "" + num + MONOMER_ATTACHEMENT_SEPARATOR
						+ tokens[1];
				setConnection(srcCon + CONNECTION_SEPARATOR + newCon);
			}
		}
	}

	public boolean isSelfCycle() {
		if (null != getSourceNode() && getSourceNode().equals(getTargetNode())) {
			return true;
		}

		return false;
	}

	public boolean isBackboneSelfCycle() {
		if (isSelfCycle()) {
			if (null != getConnection() && getConnection().contains("R1")
					&& getConnection().contains("R2")) {
				return true;
			}
		}

		return false;
	}

	public String getSourceGenericDescriptor() {
		String sourceConnection = getSourceConnection();
		if (getSourceAttachmentType() == GENERIC_EDGE_ATTACHMENT) {
			return sourceConnection.split(MONOMER_ATTACHEMENT_SEPARATOR)[1];
		}
		return null;
	}

	public String getTargetGenericDescriptor() {
		String targetConnection = getTargetConnection();
		if (getTargetAttachmentType() == GENERIC_EDGE_ATTACHMENT) {
			return targetConnection.split(MONOMER_ATTACHEMENT_SEPARATOR)[1];
		}
		return null;
	}
}
