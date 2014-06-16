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

import chemaxon.marvin.plugin.PluginException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import org.helm.notation.tools.StructureParser;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This is a data model for Monomer. alternateId is used in polymer notation.
 * 
 * @author zhangtianhong
 */
public class Monomer implements Serializable {

	public static final String NUCLIEC_ACID_POLYMER_TYPE = "RNA";
	public static final String PEPTIDE_POLYMER_TYPE = "PEPTIDE";
	public static final String CHEMICAL_POLYMER_TYPE = "CHEM";
	public static final String[] SUPPORTED_POLYMER_TYPES = {
			NUCLIEC_ACID_POLYMER_TYPE, PEPTIDE_POLYMER_TYPE,
			CHEMICAL_POLYMER_TYPE };
	public static final String BACKBONE_MOMONER_TYPE = "Backbone";
	public static final String BRANCH_MOMONER_TYPE = "Branch";
	public static final String UNDEFINED_MOMONER_TYPE = "Undefined";
	public static final String STARTING_NAME = "5";
	public static final String ATTACHMENT_LIST_DELIMITER = "$"; // surrogate
																// database ID,
																// users may
																// never see it
	private int id; // unique text ID for the monomer, will be used in polymer
					// notation
	private String alternateId; // the ID of its closest natural analog
	private String naturalAnalog; // the long name of the monomer
	private String name; // canonical SMILES that represents the monomer
	private String canSMILES;
	private String molfile; // monomer type, Backone, Branch, UnDefined
	private String monomerType; // polymer type, NucleicAcid, Peptide,
								// ChemicalStructure
	private String polymerType; // list of attachments in the monomer
	private List<Attachment> attachmentList; // mark monomer as new
	private boolean newMonomer;
	private boolean adHocMonomer;
	public static final String ID_A = "A";
	public static final String ID_G = "G";
	public static final String ID_C = "C";
	public static final String ID_U = "U";
	public static final String ID_T = "T";
	public static final String ID_R = "R";
	public static final String ID_dR = "dR";
	public static final String ID_P = "P";
	public static final String ID_X = "X";
	public static final String ID_ALA = "Ala";
	public static final String ID_ARG = "Arg";
	public static final String ID_ASN = "Asn";
	public static final String ID_ASP = "Asp";
	public static final String ID_CYS = "Cys";
	public static final String ID_GLU = "Glu";
	public static final String ID_GLN = "Gln";
	public static final String ID_GLY = "Gly";
	public static final String ID_HIS = "His";
	public static final String ID_ILE = "Ile";
	public static final String ID_LEU = "Leu";
	public static final String ID_LYS = "Lys";
	public static final String ID_MET = "Met";
	public static final String ID_PHE = "Phe";
	public static final String ID_PRO = "Pro";
	public static final String ID_SER = "Ser";
	public static final String ID_THR = "Thr";
	public static final String ID_TRP = "Trp";
	public static final String ID_TYR = "Tyr";
	public static final String ID_VAL = "Val";
	public static final String ID_CHEMICAL_STRUCTURE = "chemical structure";

	/**
	 * constructor
	 */
	public Monomer() {
		attachmentList = new ArrayList<Attachment>();
		name = "";
	}

	/**
	 * Create a new monomer.
	 * 
	 * @param polymerType
	 *            : unique text ID for the monomer, will be used in polymer
	 *            notation
	 * @param monomerType
	 *            : monomer type, Backbone, Branch, UnDefined
	 * @param naturalAnalog
	 *            : the ID of its closest natural analog
	 * @param alternateId
	 *            : unique text ID for the monomer, will be used in polymer
	 *            notation
	 */
	public Monomer(String polymerType, String monomerType,
			String naturalAnalog, String alternateId) {
		attachmentList = new ArrayList<Attachment>();
		name = "";
		setPolymerType(polymerType);
		setMonomerType(monomerType);
		setNaturalAnalog(naturalAnalog);
		setAlternateId(alternateId);
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

	/**
	 * get the natural analog of this monomer. For standard, unmodified monomer,
	 * the natural analog should be itself
	 */
	public String getNaturalAnalog() {
		if (alternateId.length() == 1) {
			return alternateId;
		} else {
			return naturalAnalog;
		}
	}

	public void setNaturalAnalog(String naturalAnalog) {
		this.naturalAnalog = naturalAnalog;
	}

	public void setAdHocMonomer(boolean adHocMonomer) {
		this.adHocMonomer = adHocMonomer;
	}

	public boolean isAdHocMonomer() {
		return this.adHocMonomer;
	}

	public String getCanSMILES() {
		return canSMILES;
	}

	public void setCanSMILES(String canSMILES) {
		this.canSMILES = canSMILES;
	}

	public String getMonomerType() {
		return monomerType;
	}

	public void setMonomerType(String monomerType) {
		this.monomerType = monomerType;
	}

	public String getPolymerType() {
		return polymerType;
	}

	public void setPolymerType(String polymerType) {
		this.polymerType = polymerType;
	}

	public List<Attachment> getAttachmentList() {
		return attachmentList;
	}

	public void setAttachmentList(List<Attachment> attachmentList) {
		this.attachmentList = attachmentList;
	}

	public String getMolfile() {
		return molfile;
	}

	public void setMolfile(String molfile) {
		this.molfile = molfile;
	}

	/**
	 * get a specific attachment by passing in a label
	 * 
	 * @param label
	 *            : unique for each attach point
	 * @return Attachment or null if there is no such attach point
	 */
	public Attachment getAttachment(String label) {
		for (Attachment attachment : attachmentList) {
			if (attachment.getLabel().equalsIgnoreCase(label)) {
				return attachment;
			}
		}
		return null;
	}

	/**
	 * This method returns the MoleculeInfo for the input R group label of this
	 * monomer
	 * 
	 * @param label
	 *            - R1, R2...
	 * @return MoleculeInfo for the cap group, R group will contribute nothing
	 * @throws IOException
	 * @throws PluginException
	 */
	public MoleculeInfo getCapMoleculeInfo(String label) throws IOException,
			PluginException {
		for (Attachment attachment : attachmentList) {
			if (attachment.getLabel().equalsIgnoreCase(label)) {
				String capSmi = attachment.getCapGroupSMILES();
				return StructureParser.getMoleculeInfo(capSmi);
			}
		}
		return null;
	}

	/**
	 * Try to add a new attachment to this monomer
	 * 
	 * @param attachment
	 *            -- new attachment to be add in
	 * @return true for success and false if there is one such attach point
	 *         exist
	 */
	public boolean addAttachment(Attachment attachment) {
		boolean isExist = false;
		for (Attachment a : attachmentList) {
			if (a.getLabel().equalsIgnoreCase(attachment.getLabel())) {
				isExist = true;
			}
		}
		if (!isExist) {
			return attachmentList.add(attachment);
		}
		return false;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Test if this monomer has been modified. The length of the alternatedId of
	 * an modified monomer should be greater than one.
	 * 
	 * @return true/false
	 */
	public boolean isModified() {
		return (alternateId.length() > 1);
	}

	/**
	 * Compare this momoner with another one, polymerType, monomerType and
	 * naturalAnalog (could be derived) must be the same to return true
	 * 
	 * @param m
	 * @return true or false
	 */
	public boolean isSameType(Monomer m) {
		String analog1 = getNaturalAnalog();
		if (null == analog1) {
			analog1 = getAlternateId();
		}
		String analog2 = m.getNaturalAnalog();
		if (null == analog2) {
			analog2 = m.getAlternateId();
		}
		if (getMonomerType().equals(m.getMonomerType())
				&& getPolymerType().equals(m.getPolymerType())
				&& analog1.equals(analog2)) {
			return true;
		} else {
			return false;
		}
	}

	public String getAttachmentListString() {
		StringBuilder sb = new StringBuilder();
		List<Attachment> al = this.getAttachmentList();
		List<String> l = new ArrayList<String>();
		for (Attachment a : al) {
			l.add(a.getAlternateId());
		}
		Collections.sort(l);

		for (int i = 0; i < l.size(); i++) {
			if (sb.length() > 0) {
				sb.append(ATTACHMENT_LIST_DELIMITER);
			}
			sb.append(l.get(i));
		}
		return sb.toString();
	}

	public boolean isNewMonomer() {
		return newMonomer;
	}

	public void setNewMonomer(boolean newMonomer) {
		this.newMonomer = newMonomer;
	}

	public boolean attachmentEquals(Monomer monomer) {
		String tmpListString = this.getAttachmentListString();
		String monomerListString = monomer.getAttachmentListString();
		return tmpListString.equals(monomerListString);
	}

	public boolean attachmentContains(Monomer monomer) {
		String tmpListString = this.getAttachmentListString();
		String monomerListString = monomer.getAttachmentListString();
		int index = tmpListString.indexOf(monomerListString);
		if (index >= 0) {
			return true;
		} else {
			return false;
		}
	}

	public boolean containAnyAtom() throws IOException {
		boolean containsA = false;
		String smiles = getCanSMILES();
		if (null != smiles && smiles.length() > 0) {
			Molecule mol = StructureParser.getMolecule(smiles);
			MolAtom[] atoms = mol.getAtomArray();
			for (MolAtom atom : atoms) {
				String symbol = atom.getSymbol();
				if ("A".equals(symbol)) {
					containsA = true;
					break;
				}
			}
		}
		return containsA;
	}
}
