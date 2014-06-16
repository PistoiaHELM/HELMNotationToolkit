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

import java.io.Serializable;

/**
 * This is a data model for Attachment point on monomer
 * 
 * @author zhangtianhong
 */
public class Attachment implements Serializable {

	public static final String BACKBONE_MONOMER_LEFT_ATTACHEMENT = "R1";
	public static final String BACKBONE_MONOMER_RIGHT_ATTACHEMENT = "R2";
	public static final String BACKBONE_MONOMER_BRANCH_ATTACHEMENT = "R3";
	public static final String BRANCH_MONOMER_ATTACHEMENT = "R1";
	public static final String STARTING_MONOMER_ATTACHEMENT = "Starting";

	public static final String CAP_GROUP_OH = "OH";
	public static final String CAP_GROUP_H = "H";

	public static final String PAIR_ATTACHMENT = "pair";

	// surrogate database ID, users may never see it
	private int id;

	// String ID for this attachment
	private String alternateId;

	// R group used for this attachment, such as R1
	private String label;

	// The name for the capping group, such as hydroxyl or OH, does not specify
	// connection atom for multi-atom groups.
	private String capGroupName;

	// The canonical SMILES for the capping group, include the R group to
	// indicate how connection should be made with monomer
	private String capGroupSMILES;

	// private boolean connected;

	public Attachment() {
		// connected = false;
	}

	public Attachment(String label, String capGroupName) {
		this.label = label;
		this.capGroupName = capGroupName;
		// connected = false;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public String getAlternateId() {
		return alternateId;
	}

	public void setAlternateId(String alternateId) {
		this.alternateId = alternateId;
	}

	public String getLabel() {
		if (label == null) {
			return " ";
		}
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public String getCapGroupName() {
		return capGroupName;
	}

	public void setCapGroupName(String capGroupName) {
		this.capGroupName = capGroupName;
	}

	public String getCapGroupSMILES() {
		return capGroupSMILES;
	}

	public void setCapGroupSMILES(String capGroupSMILES) {
		this.capGroupSMILES = capGroupSMILES;
	}

	// public boolean isConnected() {
	// return connected;
	// }
	//
	// public void setConnected(boolean connected) {
	// this.connected = connected;
	// }

	@Override
	public String toString() {
		return (label + ": " + capGroupName);
	}
}
