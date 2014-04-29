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
package org.helm.notation.tools;

import chemaxon.marvin.plugin.PluginException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.MonomerStore;
import org.helm.notation.NotationConstant;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import org.helm.notation.model.RgroupStructure;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import org.jdom.JDOMException;

/**
 * 2/24/2011 by Tianhong
 * extends to support ad hoc CHEM monomers, monomer could be added as extended ChemAxon SMILES, instead of ID
 * general process is to add it into monomer database factory first, then treat as normal monomer. The default attachment group will be Hydrogen, as R1-H, R2-H
 * 
 */

/**
 * This class provides methods that handle simple polymer notation. It is
 * unlikely that they will be called by clients directly.
 * 
 * @author zhangtianhong
 */
public class SimpleNotationParser {

	public static final String GROUP_LEVEL_DELIMITER_REGEX = "\\.";
	public static final String GROUP_LEVEL_DELIMITER = ".";
	public static final char MODIFICATION_START_SYMBOL = '[';
	public static final char MODIFICATION_END_SYMBOL = ']';

	public static final String MODIFICATION_DELIMITER_REGEX = "[\\[|\\]]";
	public static final char BRANCH_START_SYMBOL = '(';
	public static final char BRANCH_END_SYMBOL = ')';

	public static final String BRANCH_DELIMITER_REGEX = "[(|)]";
	public static final String RNA_SENSE_STRAND_ANNOTATION = "ss";
	public static final String RNA_ANTISENSE_STRAND_ANNOTATION = "as";

	/**
	 * This method generates the SMILES string for simple polymer notation
	 * 
	 * @param polymerNotation
	 * @param polymerType
	 * @return SMILES string for the simple polymer notation, still contains R
	 *         groups
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.StructureException
	 * @throws org.helm.notation.MonomerException
	 */
	public static String getSimplePolymerSMILES(String polymerNotation,
			String polymerType) throws IOException, NotationException,
			StructureException, MonomerException, JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimplePolymerSMILES(polymerNotation, polymerType,
				factory.getMonomerStore());
	}

	public static String getSimplePolymerSMILES(String polymerNotation,
			String polymerType, MonomerStore store) throws IOException,
			NotationException, StructureException, MonomerException,
			JDOMException {
		RgroupStructure struc = getSimplePolymerStructure(polymerNotation,
				polymerType, store);
		if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			if (struc.getMolecule() == null) {
				throw new NotationException(
						"Polymer notation contain non-specifc monomer structure");
			}
		}
		Molecule mol = struc.getMolecule();
		return mol.toFormat("cxsmiles:u");
	}

	/**
	 * validate chem simple notation
	 * 
	 * @param polymerNotation
	 * @return true or throw exception
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean validateSimpleNotationForChem(String polymerNotation)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return validateSimpleNotationForChem(polymerNotation,
				factory.getMonomerStore());
	}

	public static boolean validateSimpleNotationForChem(String polymerNotation,
			MonomerStore monomerStore) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {
		return validateSimpleNotation(polymerNotation,
				Monomer.CHEMICAL_POLYMER_TYPE, monomerStore);
	}

	/**
	 * validate peptide simple notation
	 * 
	 * @param polymerNotation
	 * @return true or throws exceptions
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */

	public static boolean validateSimpleNotationForPeptide(
			String polymerNotation) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return validateSimpleNotationForPeptide(polymerNotation,
				factory.getMonomerStore());
	}

	public static boolean validateSimpleNotationForPeptide(
			String polymerNotation, MonomerStore monomerStore)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {
		return validateSimpleNotation(polymerNotation,
				Monomer.PEPTIDE_POLYMER_TYPE, monomerStore);
	}

	/**
	 * validate RNA simple notation
	 * 
	 * @param polymerNotation
	 * @return true or exception
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean validateSimpleNotationForRNA(String polymerNotation)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return validateSimpleNotationForRNA(polymerNotation,
				factory.getMonomerStore());
	}

	public static boolean validateSimpleNotationForRNA(String polymerNotation,
			MonomerStore monomerStore) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {
		return validateSimpleNotation(polymerNotation,
				Monomer.NUCLIEC_ACID_POLYMER_TYPE, monomerStore);
	}

	/**
	 * This method checks if simple polymer notation is valid
	 * 
	 * @param polymerNotation
	 * @param polymerType
	 * @return true or exception
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean validateSimpleNotation(String polymerNotation,
			String polymerType) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return validateSimpleNotation(polymerNotation, polymerType,
				factory.getMonomerStore());
	}

	public static boolean validateSimpleNotation(String polymerNotation,
			String polymerType, MonomerStore monomerStore) throws IOException,
			NotationException, MonomerException, StructureException,
			JDOMException {
		getSimplePolymerStructure(polymerNotation, polymerType, monomerStore);
		return true;
	}

	/**
	 * This method generates the RgroupStructure of simple polymer notation
	 * 
	 * @param polymerNotation
	 * @param polymerType
	 * @return RgroupStructure
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 */
	public static RgroupStructure getSimplePolymerStructure(
			String polymerNotation, String polymerType) throws IOException,
			NotationException, MonomerException, StructureException,
			JDOMException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimplePolymerStructure(polymerNotation, polymerType,
				factory.getMonomerStore());
	}

	public static RgroupStructure getSimplePolymerStructure(
			String polymerNotation, String polymerType,
			MonomerStore monomerStore) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {
		List<Monomer> monomerList = getMonomerList(polymerNotation,
				polymerType, monomerStore);
		if (monomerList == null || monomerList.size() == 0) {
			throw new NotationException("Polymer notation contains no monomer");
		}
		List<RgroupStructure> structureList = getMonomerStructureList(monomerList);

		if (monomerList.size() == structureList.size()) {

			Molecule molecule = null;
			Map<String, MolAtom> rmap = new HashMap<String, MolAtom>();
			if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
				if (null != monomerList.get(0).getCanSMILES()) {

					molecule = structureList.get(0).getMolecule();
					Map rgroupMap = structureList.get(0).getRgroupMap();
					Set keyset = rgroupMap.keySet();
					for (Iterator it = keyset.iterator(); it.hasNext();) {
						String key = (String) it.next();
						rmap.put("1:" + key, (MolAtom) rgroupMap.get(key));
					}
				}

			} else {

				int prevMonomerIndex = -1;
				// Map<String, MolAtom> rgroupMap = null;

				for (int i = 0; i < monomerList.size(); i++) {
					Monomer m = monomerList.get(i);
					RgroupStructure ms = structureList.get(i);
					Molecule mol = ms.getMolecule();
					Map<String, MolAtom> rgroupMap = ms.getRgroupMap(); // modify
																		// this
																		// map
																		// when
																		// element
																		// used

					if (null != molecule) {
						RgroupStructure prevMonomerStructure = structureList
								.get(prevMonomerIndex);
						Map prevMonomerRgroupMap = prevMonomerStructure
								.getRgroupMap();

						if (m.getMonomerType().equals(
								Monomer.BACKBONE_MOMONER_TYPE)) {
							StructureParser
									.merge(molecule,
											(MolAtom) prevMonomerRgroupMap
													.get(Attachment.BACKBONE_MONOMER_RIGHT_ATTACHEMENT),
											mol,
											(MolAtom) rgroupMap
													.get(Attachment.BACKBONE_MONOMER_LEFT_ATTACHEMENT));

							prevMonomerRgroupMap
									.remove(Attachment.BACKBONE_MONOMER_RIGHT_ATTACHEMENT);
							rgroupMap
									.remove(Attachment.BACKBONE_MONOMER_LEFT_ATTACHEMENT);

							prevMonomerIndex = i;

							// possible unused R groups on previous backbone
							// monomer
							Set keySet = prevMonomerRgroupMap.keySet();
							for (Iterator it = keySet.iterator(); it.hasNext();) {
								String key = (String) it.next();
								int monomerCount = i;
								rmap.put("" + monomerCount + ":" + key,
										(MolAtom) prevMonomerRgroupMap.get(key));
							}

						} else if (m.getMonomerType().equals(
								Monomer.BRANCH_MOMONER_TYPE)) {
							StructureParser
									.merge(molecule,
											(MolAtom) prevMonomerRgroupMap
													.get(Attachment.BACKBONE_MONOMER_BRANCH_ATTACHEMENT),
											mol,
											(MolAtom) rgroupMap
													.get(Attachment.BRANCH_MONOMER_ATTACHEMENT));

							prevMonomerRgroupMap
									.remove(Attachment.BACKBONE_MONOMER_BRANCH_ATTACHEMENT);
							rgroupMap
									.remove(Attachment.BRANCH_MONOMER_ATTACHEMENT);
							// possible unused R groups on branch monomer
							Set keySet = rgroupMap.keySet();
							for (Iterator it = keySet.iterator(); it.hasNext();) {
								String key = (String) it.next();
								if (!(key
										.equals(Attachment.BRANCH_MONOMER_ATTACHEMENT))) {
									int monomerCount = i + 1;
									rmap.put("" + monomerCount + ":" + key,
											(MolAtom) rgroupMap.get(key));
								}
							}
						} else {
							throw new NotationException(
									"Undefined Monomer Type is not supported in simple polymer");
						}

					} else {
						// first monomer
						molecule = mol;
						prevMonomerIndex = i;
						rmap.put(
								"1:"
										+ Attachment.BACKBONE_MONOMER_LEFT_ATTACHEMENT,
								(MolAtom) rgroupMap
										.get(Attachment.BACKBONE_MONOMER_LEFT_ATTACHEMENT));
						rgroupMap
								.remove(Attachment.BACKBONE_MONOMER_LEFT_ATTACHEMENT);
					}
				}

				// check unused R group on the last backbone monomer
				int monomerCount = prevMonomerIndex + 1;
				Map map = structureList.get(prevMonomerIndex).getRgroupMap();
				Set keySet = map.keySet();
				for (Iterator it = keySet.iterator(); it.hasNext();) {
					String key = (String) it.next();
					rmap.put(monomerCount + ":" + key, (MolAtom) map.get(key));
				}

				// last monomer
				// if (null != rgroupMap) {
				// rmap.put("" + monomerList.size() + ":" +
				// Attachment.BACKBONE_MONOMER_RIGHT_ATTACHEMENT, (MolAtom)
				// rgroupMap.get(Attachment.BACKBONE_MONOMER_RIGHT_ATTACHEMENT));
				// }
			}

			RgroupStructure structure = new RgroupStructure();
			structure.setMolecule(molecule);
			structure.setRgroupMap(rmap);

			return structure;

		} else {
			throw new NotationException(
					"The number of monomers and structures do not match");
		}

	}

	/**
	 * This method returns Monomer object based on monomer ID and polymer type
	 * 
	 * @param monomerID
	 * @param polymerType
	 * @return Monomer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static Monomer getMonomer(String monomerID, String polymerType)
			throws NotationException {

		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMonomer(monomerID, polymerType, factory.getMonomerStore());
	}

	public static Monomer getMonomer(String monomerID, String polymerType,
			MonomerStore monomerStore) throws NotationException {

		Map<String, Monomer> monomers = monomerStore.getMonomers(polymerType);

		if (!monomers.containsKey(monomerID)) {
			throw new NotationException("Unknow monomerID '" + monomerID
					+ "' for " + polymerType);
		}
		return monomers.get(monomerID);
	}

	/**
	 * Convert list of Monomer to list of RgroupStructure
	 * 
	 * @param monomerList
	 * @return list of RgroupStructure
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 */
	private static List<RgroupStructure> getMonomerStructureList(
			List<Monomer> monomerList) throws NotationException, IOException {
		List<RgroupStructure> list = new ArrayList<RgroupStructure>();

		for (int i = 0; i < monomerList.size(); i++) {
			RgroupStructure ms = getMonomerStructure(monomerList.get(i));
			list.add(ms);
		}
		return list;
	}

	/**
	 * This methods generates the list of Monomer from polymer notation
	 * 
	 * @param polymerNotation
	 * @param polymerType
	 * @return list of Monomer
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.jdom.JDOMException
	 */
	public static List<Monomer> getMonomerList(String polymerNotation,
			String polymerType) throws IOException, NotationException,
			MonomerException, JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMonomerList(polymerNotation, polymerType,
				factory.getMonomerStore());
	}

	public static List<Monomer> getMonomerList(String polymerNotation,
			String polymerType, MonomerStore monomerStore) throws IOException,
			NotationException, MonomerException, JDOMException,
			StructureException {
		List<Monomer> list = new ArrayList<Monomer>();
		List<String> ids = getMonomerIDList(polymerNotation, polymerType,
				monomerStore);

		for (int i = 0; i < ids.size(); i++) {
			Monomer m = getMonomer(ids.get(i), polymerType, monomerStore);
			list.add(m);
		}
		return list;
	}

	/**
	 * This methods returns the nucleotide list for RNA polymer type, with
	 * structure validation on.
	 * 
	 * @param polymerNotation
	 *            - simple RNA polymer notation
	 * @return Nucleotide list
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static List<Nucleotide> getNucleotideList(String polymerNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getNucleotideList(polymerNotation, true,
				factory.getMonomerStore());
	}

	public static List<Nucleotide> getNucleotideList(String polymerNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		return getNucleotideList(polymerNotation, true, monomerStore);
	}

	/**
	 * This methods returns the nucleotide list for RNA polymer type. Validation
	 * of polymer notation could be slow
	 * 
	 * @param polymerNotation
	 *            - simple RNA polymer notation
	 * @param validate
	 *            - true will run structure validation, false will not run
	 *            structure validation.
	 * @return Nucleotide list
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static List<Nucleotide> getNucleotideList(String polymerNotation,
			boolean validate) throws NotationException, MonomerException,
			IOException, JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getNucleotideList(polymerNotation, validate,
				factory.getMonomerStore());
	}

	public static List<Nucleotide> getNucleotideList(String polymerNotation,
			boolean validate, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		if (validate) {
			validateSimpleNotationForRNA(polymerNotation, monomerStore);
		}

		List<Nucleotide> ids = new ArrayList<Nucleotide>();
		Map<String, String> reverseNucMap = NucleotideFactory.getInstance()
				.getReverseNucleotideTemplateMap();

		String[] notations = polymerNotation.split(GROUP_LEVEL_DELIMITER_REGEX);
		for (int i = 0; i < notations.length; i++) {
			String notation = notations[i];

			String symbol = null;

			// last nucleotide will be handled differently
			String tmpNotation = notation;
			if (i == (notations.length - 1) && notation.endsWith(")")) {
				tmpNotation = notation + "P";
			}

			if (reverseNucMap.containsKey(tmpNotation)) {
				symbol = reverseNucMap.get(tmpNotation);
			} else {
				// boolean modified =
				// (notation.indexOf(MODIFICATION_START_SYMBOL) >=0);
				String[] tokens = notation.split(BRANCH_DELIMITER_REGEX);
				String base = null;
				if (tokens.length >= 2) {
					base = tokens[1].replaceAll(MODIFICATION_DELIMITER_REGEX,
							"");
				}
				if (null != base) {
					if (base.length() == 1) {
						symbol = base;
					} else {
						Monomer monomer = SimpleNotationParser
								.getMonomer(base,
										Monomer.NUCLIEC_ACID_POLYMER_TYPE,
										monomerStore);
						if (null == monomer.getNaturalAnalog()) {
							symbol = "X";
						} else {
							symbol = monomer.getNaturalAnalog();
						}
					}
				} else {
					symbol = "X";
				}
			}
			Nucleotide nuc = new Nucleotide(symbol, notation);
			ids.add(nuc);
		}

		if (ids.size() > 0) {
			ids.get(0).setPositionType(Nucleotide.STARTING_POSITION_TYPE);
		}

		if (ids.size() > 1) {
			ids.get(ids.size() - 1).setPositionType(
					Nucleotide.ENDING_POSITION_TYPE);
		}

		return ids;
	}

	/**
	 * getNucleotideList is a more forgiving function that will try to
	 * substitute the natural sequence where possible to provide information for
	 * queries. This getStrictNucleotideList routine is more strict in that it
	 * will utilize modX for a natural seuqence X or include ? for unknown
	 * bases.
	 * 
	 * @param polymerNotation
	 * @param validate
	 * @return list of nucleotide
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */
	public static List<Nucleotide> getStrictNucleotideList(
			String polymerNotation, boolean validate) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getStrictNucleotideList(polymerNotation, validate,
				factory.getMonomerStore());
	}

	public static List<Nucleotide> getStrictNucleotideList(
			String polymerNotation, boolean validate, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		if (validate) {
			validateSimpleNotationForRNA(polymerNotation, monomerStore);
		}

		List<Nucleotide> ids = new ArrayList<Nucleotide>();
		Map<String, String> reverseNucMap = NucleotideFactory.getInstance()
				.getReverseNucleotideTemplateMap();
		String[] notations = polymerNotation.split(GROUP_LEVEL_DELIMITER_REGEX);
		for (int i = 0; i < notations.length; i++) {
			String notation = notations[i];
			String symbol = null;

			// last nucleotide will be handled differently
			Nucleotide nuc = null;
			if (reverseNucMap.containsKey(notation)) {
				symbol = reverseNucMap.get(notation);
				nuc = new Nucleotide(symbol, notation);
			} else {
				nuc = new Nucleotide("temp", notation);
				String naturalAnalog = nuc.getNaturalAnalog();
				if (nuc.isModified()) {
					if (i == (notations.length - 1)
							&& nuc.unmodifiedWithoutPhosphate()) {
						nuc.setSymbol("end" + naturalAnalog);
					} else {
						nuc.setSymbol("mod" + naturalAnalog);
					}
				} else {
					nuc.setSymbol(naturalAnalog);
				}
			}
			ids.add(nuc);
		}

		if (ids.size() > 0) {
			ids.get(0).setPositionType(Nucleotide.STARTING_POSITION_TYPE);
		}

		if (ids.size() > 1) {
			ids.get(ids.size() - 1).setPositionType(
					Nucleotide.ENDING_POSITION_TYPE);
		}

		return ids;
	}

	/**
	 * This method returns the list of Monomer IDs from simple polymer notation
	 * 
	 * @param polymerNotation
	 * @param polymerType
	 * @return list of Monomer IDs
	 * @throws org.helm.notation.NotationException
	 */

	public static List<String> getMonomerIDList(String polymerNotation,
			String polymerType) throws NotationException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMonomerIDList(polymerNotation, polymerType,
				factory.getMonomerStore());
	}

	public static List<String> getMonomerIDList(String polymerNotation,
			String polymerType, MonomerStore monomerStore)
			throws NotationException {
		List<String> ids = new ArrayList<String>();

		
		// CHEMICAL can have only one monomer
		if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			String id = processNode(polymerNotation, polymerType, monomerStore);
			ids.add(id);
		} else {
			/*
			 * if (polymerType.equals(Monomer.PEPTIDE_POLYMER_TYPE) &&
			 * polymerNotation.indexOf(BRANCH_START_SYMBOL) >= 0) { throw new
			 * NotationException
			 * ("Invalid Polymer Notation: Peptide polymer cannot have branches"
			 * ); }
			 */

			SimpleNotationGroupIterator groupIterator = new SimpleNotationGroupIterator(
					polymerNotation);

			while (groupIterator.hasNextGroup()) {
				String notation = groupIterator.nextGroup();
				char[] chars = notation.toCharArray();

				String curId = null;
				char prevLetter = 0;

				for (int i = 0; i < chars.length; i++) {
					char letter = chars[i];
					if (letter == MODIFICATION_START_SYMBOL) {
						int matchingPos = getMatchingBracketPosition(chars, i,
								MODIFICATION_START_SYMBOL,
								MODIFICATION_END_SYMBOL);
						i++;

						if (matchingPos == -1) {
							throw new NotationException(
									"Invalid Polymer Notation: modified monomer must be enclosed by square brackets");
						} else {
							curId = notation.substring(i, matchingPos);
						}
						i = matchingPos;
						/*
						 * while (!(String.valueOf(chars[i]).equals(
						 * MODIFICATION_END_SYMBOL))) {
						 * sb.append(String.valueOf(chars[i])); i++; }
						 */

					} else if (letter == BRANCH_START_SYMBOL) {
						if (i == 0) {
							throw new NotationException(
									"Invalid Polymer Notation: branch monomer is not allowed at the beginnig of notation");
						}

						if (prevLetter == BRANCH_END_SYMBOL) {
							throw new NotationException(
									"Invalid Polymer Notation: branch monomers cannot be connected with each other");
						}

						String string = null;

						int matchingPos = getMatchingBracketPosition(chars, i,
								BRANCH_START_SYMBOL, BRANCH_END_SYMBOL);
						i++;
						if (matchingPos == -1) {
							throw new NotationException(
									"Invalid Polymer Notation: modified monomer must be enclosed by brackets");
						} else {
							string = notation.substring(i, matchingPos);
						}
						i = matchingPos;

						/*
						 * while
						 * (!(String.valueOf(chars[i]).equals(BRANCH_END_SYMBOL
						 * ))) { sb.append(String.valueOf(chars[i])); i++; }
						 */

						if (string.length() > 1) {
							if ((string.charAt(0) == MODIFICATION_START_SYMBOL)
									&& (string.charAt(string.length() - 1) == MODIFICATION_END_SYMBOL)) {
								curId = string
										.substring(1, string.length() - 1);
							} else {
								throw new NotationException(
										"Invalid Polymer Notation: modified monomer must be enclosed by square brackets");
							}

						} else if (string.length() == 1) {
							curId = string;
						} else {
							throw new NotationException(
									"Invalid Polymer Notation: branch monomer not found");
						}
					} else {
						curId = String.valueOf(letter);
					}

					curId = processNode(curId, polymerType, monomerStore);

					ids.add(curId);
					prevLetter = letter;
				}
			}
		}
		return ids;
	}

	protected static String processNode(String nodeDesc, String polymerType,
			MonomerStore monomerStore) throws NotationException {

		//remove brackets 
		if ((nodeDesc.charAt(0)==MODIFICATION_START_SYMBOL) && (nodeDesc.charAt(nodeDesc.length()-1)==MODIFICATION_END_SYMBOL)){
			nodeDesc=nodeDesc.substring(1, nodeDesc.length()-1);
		}
		
		boolean isSmilesCode = Pattern.matches(".*\\$\\|$", nodeDesc);

		if (!isSmilesCode) {
			Map<String, Monomer> monomers = monomerStore.getMonomers(polymerType);
			if (monomers.containsKey(nodeDesc)){
				return nodeDesc;
			
			}
		}

	Map<String, Monomer> smilesDB = monomerStore.getSmilesMonomerDB();

		String alternateId = null;

		if (smilesDB.containsKey(nodeDesc)) {
			alternateId = smilesDB.get(nodeDesc).getAlternateId();

		} else {

			MonomerFactory factory;
			try {
				factory = MonomerFactory.getInstance();
			} catch (Exception ex) {
				throw new NotationException(
						"Unable to initialize monomer factory", ex);
			}

			alternateId = generateNextAdHocMonomerID(polymerType, monomerStore);

			Map<String, Attachment> ids = factory.getAttachmentDB();
			Attachment R1HAtt = ids.get("R1-H");

			Monomer m = null;
			if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
				m = new Monomer(polymerType, Monomer.UNDEFINED_MOMONER_TYPE,
						null, alternateId);
			}
			else if (polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
				m = new Monomer(polymerType, Monomer.BRANCH_MOMONER_TYPE,
						"X", alternateId);
			} 
			else {
				m = new Monomer(polymerType, Monomer.BACKBONE_MOMONER_TYPE,
						"X", alternateId);
			}
			m.setAdHocMonomer( true);
			m.setCanSMILES(nodeDesc);

			List<Attachment> al = new ArrayList<Attachment>();
			int start = 0;
			int pos = nodeDesc.indexOf("R", start);
			String number = "";
			while (pos >= 0) {
				pos++;
				String letter = nodeDesc.substring(pos, pos + 1);
				while (letter.matches("\\d")) {
					number = number + letter;
					pos++;
					letter = nodeDesc.substring(pos, pos + 1);
				}

				try {
					Attachment tmpAtt = DeepCopy.copy(R1HAtt);
					tmpAtt.setLabel("R" + number);
					tmpAtt.setAlternateId("R" + number + "-H");
					String oldSmi = tmpAtt.getCapGroupSMILES();
					String newSmi = oldSmi.replace("R1", "R" + number);
					tmpAtt.setCapGroupSMILES(newSmi);
					al.add(tmpAtt);
				} catch (Exception ex) {
					throw new NotationException(
							"Unable to create attachment by copying from attachment database",
							ex);
				}

				start = pos;
				pos = nodeDesc.indexOf("R", start);
				number = "";
			}

			m.setAttachmentList(al);
			try {
				monomerStore.addNewMonomer(m);
			} catch (Exception ex) {
				throw new NotationException(
						"Unable to add adhoc new monomer into monomer databse",
						ex);
			}
		}
		return alternateId;
	}

	/**
	 * This method returns the RgroupStructure object for the monomer. R group
	 * ID must be unique in each monomer
	 * 
	 * @param monomer
	 * @return RgroupStructure
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 */
	private static RgroupStructure getMonomerStructure(Monomer monomer)
			throws IOException, NotationException {
		RgroupStructure ms = null;
		if (null != monomer.getCanSMILES()) {
			ms = new RgroupStructure();
			Molecule mol = StructureParser.getMolecule(monomer.getCanSMILES());
			mol.dearomatize();
			Map<String, MolAtom> rgroupMap = new HashMap<String, MolAtom>();

			MolAtom[] atoms = mol.getAtomArray();
			for (int i = 0; i < atoms.length; i++) {
				int rId = atoms[i].getRgroup();

				// found R group
				if (rId > 0) {
					rgroupMap.put("R" + rId, atoms[i]);
				}
			}
			ms.setMolecule(mol);
			ms.setRgroupMap(rgroupMap);
		}

		return ms;
	}

	/**
	 * generate single letter (natural analog) nucleotide sequence for simple
	 * RNA notation, remove nucleotide without base from both ends
	 * 
	 * @param polymerNotation
	 * @return single letter (natural analog) sequence
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	
	public static String getTrimmedNucleotideSequence(String polymerNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getTrimmedNucleotideSequence(polymerNotation,factory.getMonomerStore());
	}
	
		
	public static String getTrimmedNucleotideSequence(String polymerNotation,MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		List<Nucleotide> list = getNucleotideList(polymerNotation,monomerStore);

		int start = 0;
		Nucleotide na = list.get(start);
		while (null == na.getBaseMonomer(monomerStore)) {
			start++;
			na = list.get(start);
		}

		int end = list.size() - 1;
		na = list.get(end);
		while (null == na.getBaseMonomer(monomerStore)) {
			end--;
			na = list.get(end);
		}

		StringBuffer sb = new StringBuffer();
		for (int i = start; i <= end; i++) {
			sb.append(list.get(i).getNaturalAnalog());
		}
		return sb.toString();
	}

	/**
	 * generate single letter (natural analog) nucleotide sequence for simple
	 * RNA notatoin
	 * 
	 * @param polymerNotation
	 * @return single letter (natural analog) sequence
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	
	public static String getNucleotideSequence(String polymerNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		MonomerFactory factory = null;
    	try {
    		factory = MonomerFactory.getInstance();
    	} catch (Exception ex) {
    		throw new NotationException("Unable to initialize monomer factory",
    				ex);
    	}
    	return getNucleotideSequence(polymerNotation,factory.getMonomerStore());
	}
	
	
	public static String getNucleotideSequence(String polymerNotation,MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		List<Nucleotide> list = getNucleotideList(polymerNotation,monomerStore);
		return getNucleotideSequence(list);
	}

	protected static String getNucleotideSequence(
			List<Nucleotide> nucleotideList) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < nucleotideList.size(); i++) {
			sb.append(nucleotideList.get(i).getNaturalAnalog());
		}
		return sb.toString();
	}

	/**
	 * convert RNA simple notation to modified nucleotide sequence
	 * 
	 * @param polymerNotation
	 * @return modified nucleotide sequence string
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static String getModifiedNucleotideSequence(String polymerNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		
		MonomerFactory factory = null;
    	try {
    		factory = MonomerFactory.getInstance();
    	} catch (Exception ex) {
    		throw new NotationException("Unable to initialize monomer factory",
    				ex);
    	}
    	return getModifiedNucleotideSequence(polymerNotation,factory.getMonomerStore());
		
	}
	
	public static String getModifiedNucleotideSequence(String polymerNotation,MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		List<Nucleotide> list = getNucleotideList(polymerNotation,monomerStore);
		return getModifiedNucleotideSequence(list);
	}

	protected static String getModifiedNucleotideSequence(
			List<Nucleotide> nucleotideList) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < nucleotideList.size(); i++) {
			sb.append(nucleotideList.get(i).getSymbol());
		}
		return sb.toString();
	}

	/**
	 * convert chemical simple notation to complex notation
	 * 
	 * @param simpleNotation
	 * @return complext notation for CHEM polymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplextNotationForChem(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getComplextNotationForChem(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getComplextNotationForChem(String simpleNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, StructureException, JDOMException, IOException {
		return getComplexNotation(simpleNotation,
				Monomer.CHEMICAL_POLYMER_TYPE, monomerStore);
	}

	/**
	 * convert Peptide simple notation to complex notation
	 * 
	 * @param simpleNotation
	 * @return complext notation for PEPTIDE polymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplextNotationForPeptide(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getComplextNotationForPeptide(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getComplextNotationForPeptide(String simpleNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, StructureException, JDOMException, IOException {
		return getComplexNotation(simpleNotation, Monomer.PEPTIDE_POLYMER_TYPE,
				monomerStore);
	}

	/**
	 * convert RNA simple notation to complex notation and annotate as antisense
	 * strand 'as'
	 * 
	 * @param simpleNotation
	 * @return complext notation for RNA polymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplextNotationForAntisenseRNA(
			String simpleNotation) throws NotationException, MonomerException,
			StructureException, JDOMException, IOException {
		MonomerFactory factory = null;
    	try {
    		factory = MonomerFactory.getInstance();
    	} catch (Exception ex) {
    		throw new NotationException("Unable to initialize monomer factory",
    				ex);
    	}
    	return getComplextNotationForAntisenseRNA(simpleNotation,factory.getMonomerStore());
	}
	
	
	public static String getComplextNotationForAntisenseRNA(
			String simpleNotation,MonomerStore monomerStore) throws NotationException, MonomerException,
			StructureException, JDOMException, IOException {
		String notation = getComplextNotationForRNA(simpleNotation,monomerStore);
		return notation.substring(0, notation.length() - 1)
				+ Monomer.NUCLIEC_ACID_POLYMER_TYPE + "1{"
				+ RNA_ANTISENSE_STRAND_ANNOTATION + "}$";
	}

	/**
	 * convert RNA simple notation to complex notation and annotate as sense
	 * strand 'ss'
	 * 
	 * @param simpleNotation
	 * @return complext notation for RNA polymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplextNotationForSenseRNA(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
    	try {
    		factory = MonomerFactory.getInstance();
    	} catch (Exception ex) {
    		throw new NotationException("Unable to initialize monomer factory",
    				ex);
    	}
		return getComplextNotationForSenseRNA(simpleNotation,factory.getMonomerStore());
		
	}
	
	
	public static String getComplextNotationForSenseRNA(String simpleNotation,MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		String notation = getComplextNotationForRNA(simpleNotation,monomerStore);
		return notation.substring(0, notation.length() - 1)
				+ Monomer.NUCLIEC_ACID_POLYMER_TYPE + "1{"
				+ RNA_SENSE_STRAND_ANNOTATION + "}$";
	}

	/**
	 * convert RNA simple notation to complex notation
	 * 
	 * @param simpleNotation
	 * @return complext notation for RNA polymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplextNotationForRNA(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getComplextNotationForRNA(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getComplextNotationForRNA(String simpleNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, StructureException, JDOMException, IOException {
		return getComplexNotation(simpleNotation,
				Monomer.NUCLIEC_ACID_POLYMER_TYPE, monomerStore);
	}

	/**
	 * This method converts simple notation for a given polymer type to complex
	 * notation. Complex notation is needed to calculate molecule info such as
	 * MW and formula
	 * 
	 * @param simpleNotation
	 * @param polymerType
	 * @return complex notation
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String getComplexNotation(String simpleNotation,
			String polymerType) throws NotationException, MonomerException,
			StructureException, JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getComplexNotation(simpleNotation, polymerType,
				factory.getMonomerStore());
	}

	public static String getComplexNotation(String simpleNotation,
			String polymerType, MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		validateSimpleNotation(simpleNotation, polymerType, monomerStore);
		return polymerType + "1{" + simpleNotation + "}$$$$";
	}

	/**
	 * This method returns the total number of monomers in peptide simple
	 * notation
	 * 
	 * @param simpleNotation
	 * @return the total number of monomers in the notation
	 * @throws org.helm.notation.NotationException
	 */
	public static int getMonomerCountForPeptide(String simpleNotation)
			throws NotationException {

		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMonomerCountForPeptide(simpleNotation,
				factory.getMonomerStore());
	}

	public static int getMonomerCountForPeptide(String simpleNotation,
			MonomerStore monomerStore) throws NotationException {
		return getMonomerCount(simpleNotation, Monomer.PEPTIDE_POLYMER_TYPE,
				monomerStore);
	}

	/**
	 * This method returns the total number of monomers in RNA simple notation
	 * 
	 * @param simpleNotation
	 * @return the total number of monomers in the notation
	 * @throws org.helm.notation.NotationException
	 */
	public static int getMonomerCountForRNA(String simpleNotation)
			throws NotationException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMonomerCountForRNA(simpleNotation, factory.getMonomerStore());
	}

	public static int getMonomerCountForRNA(String simpleNotation,
			MonomerStore monomerStore) throws NotationException {
		return getMonomerCount(simpleNotation,
				Monomer.NUCLIEC_ACID_POLYMER_TYPE, monomerStore);
	}

	/**
	 * This method returns the total number of monomers in a simple notation for
	 * a give polymer type
	 * 
	 * @param simpleNotation
	 * @param polymerType
	 * @return the total number of monomers in the notation
	 * @throws org.helm.notation.NotationException
	 */

	public static int getMonomerCount(String simpleNotation, String polymerType)
			throws NotationException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}

		return getMonomerCount(simpleNotation, polymerType,
				factory.getMonomerStore());

	}

	public static int getMonomerCount(String simpleNotation,
			String polymerType, MonomerStore monomerStore)
			throws NotationException {
		List<String> idList = getMonomerIDList(simpleNotation, polymerType,
				monomerStore);
		return idList.size();
	}

	/**
	 * This method replaces existing monomer with new monomer in simple polymer
	 * notation
	 * 
	 * @param simpleNotation
	 * @param polymerType
	 * @param existingMonomerID
	 * @param newMonomerID
	 * @return simple notation after replacement
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.NotationException
	 */
	public static String replaceMonomer(String simpleNotation,
			String polymerType, String existingMonomerID, String newMonomerID)
			throws MonomerException, IOException, JDOMException,
			NotationException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return replaceMonomer(simpleNotation, polymerType, existingMonomerID,
				newMonomerID, factory.getMonomerStore(), true);
	}

	public static String replaceMonomer(String simpleNotation,
			String polymerType, String existingMonomerID, String newMonomerID,
			MonomerStore monomerStore, boolean validate)
			throws MonomerException, IOException, JDOMException,
			NotationException {
		if (validate) {
			boolean valid = validateMonomerReplacement(polymerType,
					existingMonomerID, newMonomerID, monomerStore);
		}
		List<String> monomerIDs = getMonomerIDList(simpleNotation, polymerType,
				monomerStore);
		for (int i = 0; i < monomerIDs.size(); i++) {
			String tmp = monomerIDs.get(i);
			if (tmp.equals(existingMonomerID)) {
				monomerIDs.set(i, newMonomerID);
			}
		}
		return getSimpleNotation(monomerIDs, polymerType, monomerStore);
	}

	private static String getSimpleNotation(List<String> monomerIDList,
			String polymerType) throws NotationException, MonomerException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleNotation(monomerIDList, polymerType,
				factory.getMonomerStore());
	}

	private static String getSimpleNotation(List<String> monomerIDList,
			String polymerType, MonomerStore monomerStore)
			throws NotationException, MonomerException, JDOMException,
			IOException {
		StringBuffer sb = new StringBuffer();
		if (polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
			for (String id : monomerIDList) {
				Monomer monomer = getMonomer(id, polymerType, monomerStore);
				String monomerType = monomer.getMonomerType();
				String naturalAnalog = monomer.getNaturalAnalog();

				if (monomerType.equals(Monomer.BACKBONE_MOMONER_TYPE)) {
					if (naturalAnalog.equals(Monomer.ID_R)) {
						if (sb.length() > 0) {
							sb.append(SimpleNotationParser.GROUP_LEVEL_DELIMITER);
						}
					}
					if (id.length() > 1) {
						sb.append(SimpleNotationParser.MODIFICATION_START_SYMBOL);
						sb.append(id);
						sb.append(SimpleNotationParser.MODIFICATION_END_SYMBOL);
					} else {
						sb.append(id);
					}

				} else {
					sb.append(SimpleNotationParser.BRANCH_START_SYMBOL);
					if (id.length() > 1) {
						sb.append(SimpleNotationParser.MODIFICATION_START_SYMBOL);
						sb.append(id);
						sb.append(SimpleNotationParser.MODIFICATION_END_SYMBOL);
					} else {
						sb.append(id);
					}
					sb.append(SimpleNotationParser.BRANCH_END_SYMBOL);
				}
			}

		} else if (polymerType.equals(Monomer.PEPTIDE_POLYMER_TYPE)) {
			for (String id : monomerIDList) {
				if (sb.length() > 0) {
					sb.append(SimpleNotationParser.GROUP_LEVEL_DELIMITER);
				}
				if (id.length() > 1) {
					sb.append(SimpleNotationParser.MODIFICATION_START_SYMBOL);
					sb.append(id);
					sb.append(SimpleNotationParser.MODIFICATION_END_SYMBOL);
				} else {
					sb.append(id);
				}
			}
		} else if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			sb.append(monomerIDList.get(0));
		}
		return sb.toString();
	}

	public static String getSimpleNotationFromNucleotideList(
			List<Nucleotide> nucleotideList) {
		StringBuffer sb = new StringBuffer();
		for (Nucleotide nuc : nucleotideList) {

			if (sb.length() > 0) {
				sb.append(SimpleNotationParser.GROUP_LEVEL_DELIMITER);
			}
			String notion = nuc.getNotation();
			sb.append(notion);
		}
		return sb.toString();
	}

	protected static boolean validateMonomerReplacement(String polymerType,
			String existingMonomerID, String newMonomerID)
			throws MonomerException, IOException, JDOMException,
			NotationException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return validateMonomerReplacement(polymerType, existingMonomerID,
				newMonomerID, factory.getMonomerStore());
	}

	protected static boolean validateMonomerReplacement(String polymerType,
			String existingMonomerID, String newMonomerID,
			MonomerStore monomerStore) throws MonomerException, IOException,
			JDOMException, NotationException {

		if (null == polymerType || polymerType.length() == 0) {
			throw new NotationException(
					"Polymer type is required for monomer replacement");
		}

		if (null == existingMonomerID || existingMonomerID.length() == 0) {
			throw new NotationException(
					"Existing monomer ID is required for monomer replacement");
		}

		if (null == newMonomerID || newMonomerID.length() == 0) {
			throw new NotationException(
					"New monomer ID is required for monomer replacement");
		}
		Map<String, Monomer> monomers = monomerStore.getMonomers(polymerType);
		// Map<String, Monomer> monomers =
		// MonomerFactory.getInstance().getMonomerDB().get(polymerType);
		if (null == monomers || monomers.size() == 0) {
			throw new NotationException("Unknown polymer type [" + polymerType
					+ "] found");
		}

		if (!monomers.containsKey(existingMonomerID)) {
			throw new NotationException("Existing monomer ID ["
					+ existingMonomerID + "] is invalid in polymer type "
					+ polymerType);
		}

		if (!monomers.containsKey(newMonomerID)) {
			throw new NotationException("New monomer ID [" + newMonomerID
					+ "] is invalid in polymer type " + polymerType);
		}

		Monomer existingMonomer = monomers.get(existingMonomerID);
		Monomer newMonomer = monomers.get(newMonomerID);

		if (polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
			if (existingMonomer.getMonomerType().equals(
					newMonomer.getMonomerType())) {
				if (existingMonomer.getMonomerType().equals(
						Monomer.BACKBONE_MOMONER_TYPE)) {
					if (!existingMonomer.getNaturalAnalog().equals(
							newMonomer.getNaturalAnalog())) {
						throw new NotationException(
								"Existing monomer natural analog ["
										+ existingMonomer.getNaturalAnalog()
										+ "] and new monomer natural analog ["
										+ newMonomer.getNaturalAnalog()
										+ "] are different");
					}
				}
			} else {
				throw new NotationException("Existing monomer type ["
						+ existingMonomer.getMonomerType()
						+ "] and new monomer type ["
						+ newMonomer.getMonomerType() + "] are different");
			}
		}

		if (!newMonomer.attachmentEquals(existingMonomer)) {
			throw new NotationException("Existing monomer attachment ["
					+ existingMonomer.getAttachmentListString()
					+ "] and new monomer attachement ["
					+ newMonomer.getAttachmentListString() + "] are different");
		}

		return true;
	}

	/**
	 * This method returns the single letter natural analog peptide sequence
	 * 
	 * @param polymerNotation
	 *            -- simple peptide notation
	 * @return single letter peptide sequence
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.StructureException
	 */
	public static String getPeptideSequence(String polymerNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getPeptideSequence(polymerNotation, factory.getMonomerStore());
	}

	public static String getPeptideSequence(String polymerNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		List<Monomer> list = getMonomerList(polymerNotation,
				Monomer.PEPTIDE_POLYMER_TYPE, monomerStore);
		StringBuffer sb = new StringBuffer();
		for (Monomer monomer : list) {
			String single = monomer.getNaturalAnalog();
			if (null == single) {
				single = "X";
			}
			sb.append(single);
		}
		return sb.toString();
	}

	/**
	 * This method returns the modified peptide sequence
	 * 
	 * @param polymerNotation
	 *            -- simple peptide notation
	 * @param delimiter
	 *            used to separate amino acid
	 * @return modified amino acid sequenc with specified delimiter
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.StructureException
	 */
	public static String getModifiedPeptideSequence(String polymerNotation,
			String delimiter) throws NotationException, MonomerException,
			IOException, JDOMException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getModifiedPeptideSequence(polymerNotation, delimiter,
				factory.getMonomerStore());
	}

	public static String getModifiedPeptideSequence(String polymerNotation,
			String delimiter, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		List<Monomer> list = getMonomerList(polymerNotation,
				Monomer.PEPTIDE_POLYMER_TYPE, monomerStore);
		StringBuffer sb = new StringBuffer();
		String separator = "";
		if (null == delimiter) {
			separator = "";
		} else {
			separator = delimiter;
		}

		for (Monomer monomer : list) {
			if (sb.length() > 0) {
				sb.append(separator);
			}
			sb.append(monomer.getAlternateId());
		}
		return sb.toString();
	}

	public static String getSimpleCanonicalNotationForRNA(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleCanonicalNotationForRNA(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getSimpleCanonicalNotationForRNA(
			String simpleNotation, MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		validateSimpleNotation(simpleNotation,
				Monomer.NUCLIEC_ACID_POLYMER_TYPE, monomerStore);

		// check for cycle info
		// if first nucleotide starts with sugar monomer, last nucloetide must
		// end with phosphate monomer
		// if first nucleotide starts with phosphate monomer, last nculeotide
		// must not contain phosphate monomer
		List<Nucleotide> nList = getNucleotideList(simpleNotation, false,
				monomerStore);
		Nucleotide firstNuc = nList.get(0);
		Nucleotide lastNuc = nList.get(nList.size() - 1);
		String standardNotation = simpleNotation;
		if (firstNuc.getSugarMonomer() != null) {
			if (lastNuc.getPhosphateMonomer() == null) {
				throw new NotationException(
						"Non-cyclizable nucleic acid: last nucleotide contains no phosphate linker");
			}
		} else {
			if (lastNuc.getPhosphateMonomer() != null) {
				throw new NotationException(
						"Non-cyclizable nucleic acid: duplicate phosphate linkers at both ends");
			} else {
				int firstDelimiter = simpleNotation.indexOf(
						GROUP_LEVEL_DELIMITER, 0);
				String backPart = simpleNotation.substring(firstDelimiter + 1);
				String frontPart = simpleNotation.substring(0, firstDelimiter);
				standardNotation = backPart + frontPart;
			}
		}

		return getSimpleCanonicalNotation(standardNotation);
	}

	public static String getSimpleCanonicalNotationForPeptide(
			String simpleNotation) throws NotationException, MonomerException,
			StructureException, JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleCanonicalNotationForPeptide(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getSimpleCanonicalNotationForPeptide(
			String simpleNotation, MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		validateSimpleNotation(simpleNotation, Monomer.PEPTIDE_POLYMER_TYPE,
				monomerStore);
		return getSimpleCanonicalNotation(simpleNotation);
	}

	public static String getSimpleCanonicalNotationForChem(String simpleNotation)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleCanonicalNotationForChem(simpleNotation,
				factory.getMonomerStore());
	}

	public static String getSimpleCanonicalNotationForChem(
			String simpleNotation, MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		validateSimpleNotation(simpleNotation, Monomer.CHEMICAL_POLYMER_TYPE,
				monomerStore);
		return getSimpleCanonicalNotation(simpleNotation);
	}

	public static String getSimpleCanonicalNotation(String simpleNotation,
			String polymerType) throws NotationException, MonomerException,
			StructureException, JDOMException, IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleCanonicalNotation(simpleNotation, polymerType,
				factory.getMonomerStore());
	}

	public static String getSimpleCanonicalNotation(String simpleNotation,
			String polymerType, MonomerStore monomerStore)
			throws NotationException, MonomerException, StructureException,
			JDOMException, IOException {
		if (Monomer.NUCLIEC_ACID_POLYMER_TYPE.equals(polymerType)) {
			return getSimpleCanonicalNotationForRNA(simpleNotation,
					monomerStore);
		} else if (Monomer.PEPTIDE_POLYMER_TYPE.equals(polymerType)) {
			return getSimpleCanonicalNotationForPeptide(simpleNotation,
					monomerStore);
		} else if (Monomer.CHEMICAL_POLYMER_TYPE.equals(polymerType)) {
			return getSimpleCanonicalNotationForChem(simpleNotation,
					monomerStore);
		}
		throw new NotationException("Unknown Polymer Type: " + polymerType);
	}

	private static String getSimpleCanonicalNotation(String simpleNotation) {
		List<String> list = new ArrayList<String>();
		int start = 0;
		list.add(simpleNotation);
		int current = simpleNotation.indexOf(GROUP_LEVEL_DELIMITER, start);

		while (current > 0) {
			String backPart = simpleNotation.substring(current + 1);
			String frontPart = simpleNotation.substring(0, current);
			String newNotation = backPart + GROUP_LEVEL_DELIMITER + frontPart;
			list.add(newNotation);
			start = current + 1;
			current = simpleNotation.indexOf(GROUP_LEVEL_DELIMITER, start);
		}

		Collections.sort(list);
		return list.get(0);
	}

	public static Map.Entry<Integer, String> getSimpleCanonicalNotationMapEntry(
			String simpleNotation, String polymerType)
			throws NotationException, MonomerException, JDOMException,
			IOException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getSimpleCanonicalNotationMapEntry(simpleNotation, polymerType,
				factory.getMonomerStore());
	}

	public static Map.Entry<Integer, String> getSimpleCanonicalNotationMapEntry(
			String simpleNotation, String polymerType, MonomerStore monomerStore)
			throws NotationException, MonomerException, JDOMException,
			IOException, StructureException {
		List<String> list = new ArrayList<String>();
		if (Monomer.NUCLIEC_ACID_POLYMER_TYPE.equals(polymerType)) {
			List<Nucleotide> nucList = getNucleotideList(simpleNotation, false,
					monomerStore);
			Nucleotide firstNuc = nucList.get(0);
			Nucleotide lastNuc = nucList.get(nucList.size() - 1);
			int offset = 0;
			if (firstNuc.getSugarMonomer() != null) {
				if (lastNuc.getPhosphateMonomer() == null) {
					throw new NotationException(
							"Non-cyclizable nucleic acid: last nucleotide contains no phosphate linker");
				}
			} else {
				if (lastNuc.getPhosphateMonomer() != null) {
					throw new NotationException(
							"Non-cyclizable nucleic acid: duplicate phosphate linkers at both ends");
				} else {
					String firstNucNotation = firstNuc.getNotation();
					String lastNucNotation = lastNuc.getNotation();
					String newNotation = lastNucNotation + firstNucNotation;
					lastNuc.setNotation(newNotation);
					nucList.remove(0);
					offset = 1;
				}
			}

			for (int i = 0; i < nucList.size(); i++) {
				String notation = getSimpleNotationFromNucleotideList(nucList);
				list.add(notation + " " + offset);
				Nucleotide nuc = nucList.get(0);
				if (null != nuc.getSugarMonomer()) {
					offset++;
				}
				if (null != nuc.getBaseMonomer(monomerStore)) {
					offset++;
				}
				if (null != nuc.getPhosphateMonomer()) {
					offset++;
				}
				nucList.remove(0);
				nucList.add(nuc);
			}

		} else if (Monomer.PEPTIDE_POLYMER_TYPE.equals(polymerType)
				|| Monomer.CHEMICAL_POLYMER_TYPE.equals(polymerType)) {
			List<String> monomerIDList = getMonomerIDList(simpleNotation,
					polymerType, monomerStore);
			List<String> tmpList = new ArrayList<String>();
			for (String id : monomerIDList) {
				tmpList.add(id);
			}

			for (int i = 0; i < monomerIDList.size(); i++) {
				String monomerId = monomerIDList.get(i);
				String notation = getSimpleNotation(tmpList, polymerType,
						monomerStore);
				list.add(notation + " " + i);
				tmpList.remove(0);
				tmpList.add(monomerId);
			}
		}

		if (list.size() > 0) {
			Collections.sort(list);
			String entryString = list.get(0);
			String[] entries = entryString.split("\\s");
			String canNotation = entries[0];
			Integer position = new Integer(entries[1]);
			Map<Integer, String> map = new HashMap<Integer, String>();
			map.put(position, canNotation);
			return map.entrySet().iterator().next();
		}

		throw new NotationException("Unknown Polymer Type: " + polymerType);
	}

	/**
	 * This method returns the MoleculeInfo of simple polymer using a divide and
	 * conquer approach
	 * 
	 * @param notation
	 *            - simple notation
	 * @param polymerType
	 *            - RNA, PEPTIDE or CHEM
	 * @return MoleculeInfo of this simple polymer, all R groups are capped.
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws PluginException
	 * @throws StructureException
	 */
	public static MoleculeInfo getMoleculeInfo(String notation,
			String polymerType) throws NotationException, MonomerException,
			IOException, JDOMException, PluginException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return getMoleculeInfo(notation, polymerType, factory.getMonomerStore());
	}

	public static MoleculeInfo getMoleculeInfo(String notation,
			String polymerType, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, PluginException, StructureException {

		if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			notation = processNode(notation, Monomer.CHEMICAL_POLYMER_TYPE,
					monomerStore);
			// notation = preprocessChemNode(notation,monomerStore);
		}

		List<String> chunks = new ArrayList<String>();

		StringBuilder tmpNotationBuilder = new StringBuilder();

		SimpleNotationGroupIterator groupIterator = new SimpleNotationGroupIterator(
				notation);
		List<String> groups = new ArrayList<String>();
		// String[] groups = notation.split(GROUP_LEVEL_DELIMITER_REGEX);

		while (groupIterator.hasNextGroup()) {
			groups.add(groupIterator.nextGroup());

		}

		for (int i = 1; i <= groups.size(); i++) {
			if (tmpNotationBuilder.length() > 0) {
				tmpNotationBuilder.append(GROUP_LEVEL_DELIMITER);
			}
			tmpNotationBuilder.append(groups.get(i - 1));

			if ((i % NotationConstant.MONOMER_GROUP_COUNT_INTERVAL) == 0) {
				chunks.add(tmpNotationBuilder.toString());
				tmpNotationBuilder = new StringBuilder();
			}
		}

		if (tmpNotationBuilder.length() > 0) {
			chunks.add(tmpNotationBuilder.toString());
		}

		List<MoleculeInfo> capMiList = new ArrayList<MoleculeInfo>();
		// R2 caps
		for (int i = 0; i < chunks.size() - 1; i++) {
			String tmpNotation = chunks.get(i);
			List<String> monomerIDList = getMonomerIDList(tmpNotation,
					polymerType, monomerStore);
			Monomer monomer = getMonomer(
					monomerIDList.get(monomerIDList.size() - 1), polymerType,
					monomerStore);
			capMiList.add(monomer.getCapMoleculeInfo("R2"));
		}

		// R1 caps
		for (int i = 1; i < chunks.size(); i++) {
			String tmpNotation = chunks.get(i);
			List<String> monomerIDList = getMonomerIDList(tmpNotation,
					polymerType, monomerStore);
			Monomer monomer = getMonomer(monomerIDList.get(0), polymerType,
					monomerStore);
			capMiList.add(monomer.getCapMoleculeInfo("R1"));
		}

		if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			String tmpNotation = chunks.get(0);
			List<String> monomerIDList = getMonomerIDList(tmpNotation,
					polymerType, monomerStore);
			Monomer monomer = getMonomer(monomerIDList.get(0), polymerType,
					monomerStore);
			List<Attachment> attachments = monomer.getAttachmentList();
			for (Attachment att : attachments) {
				String label = att.getLabel();
				if (label.equals("R1") || label.equals("R2")) {
					// ignore
				} else {
					capMiList.add(monomer.getCapMoleculeInfo(label));
				}
			}
		}

		List<MoleculeInfo> chunkMiList = new ArrayList<MoleculeInfo>();
		for (int i = 0; i < chunks.size(); i++) {
			String tmpNotation = chunks.get(i);
			String complexNotation = getComplexNotation(tmpNotation,
					polymerType, monomerStore);
			String smiles = ComplexNotationParser.getComplexPolymerSMILES(
					complexNotation, monomerStore);
			MoleculeInfo tmpmi = StructureParser.getMoleculeInfo(smiles);
			chunkMiList.add(tmpmi);
		}

		return StructureParser.processMoleculeInfo(chunkMiList, capMiList);
	}

	/**
	 * This method preprocess the CHEMCHEM monomer node. Monomer could be
	 * predefined or added ad hoc which the system will generate an ID, and add
	 * to monomer factory
	 * 
	 * @param nodeDesc
	 *            - in the format of ID or extended smiles
	 * @return monomer alternateID
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 */
//	protected static String preprocessChemNode(String nodeDesc)
//			throws NotationException {
//		MonomerFactory factory;
//		try {
//			factory = MonomerFactory.getInstance();
//		} catch (Exception ex) {
//			throw new NotationException("Unable to initialize monomer factory",
//					ex);
//		}
//		return preprocessChemNode(nodeDesc, factory.getMonomerStore());
//	}
//
//	protected static String preprocessChemNode(String nodeDesc,
//			MonomerStore monomerStore) throws NotationException {
//
//		Map<String, Monomer> chemMonomers = monomerStore
//				.getMonomers(Monomer.CHEMICAL_POLYMER_TYPE);
//		// Map<String, Monomer> externalChemMonomers =
//		// factory.getExternalMonomerDB().get(Monomer.CHEMICAL_POLYMER_TYPE);
//		Map<String, Monomer> smilesDB = monomerStore.getSmilesMonomerDB();
//
//		String alternateId = null;
//		/*
//		 * if (externalChemMonomers!=null &&
//		 * externalChemMonomers.containsKey(nodeDesc)) { alternateId = nodeDesc;
//		 * } else
//		 */
//		if (chemMonomers.containsKey(nodeDesc)) {
//			alternateId = nodeDesc;
//		} else {
//			if (smilesDB.containsKey(nodeDesc)) {
//				Monomer tempM = smilesDB.get(nodeDesc);
//				if (tempM.getPolymerType()
//						.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
//					alternateId = tempM.getAlternateId();
//				} else {
//					throw new NotationException(
//							"Ad Hoc chemical monomer structure belongs to "
//									+ tempM.getPolymerType());
//				}
//			} else {
//				MonomerFactory factory;
//				try {
//					factory = MonomerFactory.getInstance();
//				} catch (Exception ex) {
//					throw new NotationException(
//							"Unable to initialize monomer factory", ex);
//				}
//
//				alternateId = generateNextMonomerID(
//						Monomer.CHEMICAL_POLYMER_TYPE, monomerStore);
//
//				Map<String, Attachment> ids = factory.getAttachmentDB();
//				Attachment R1HAtt = ids.get("R1-H");
//				Monomer m = new Monomer(Monomer.CHEMICAL_POLYMER_TYPE,
//						Monomer.UNDEFINED_MOMONER_TYPE, null, alternateId);
//				m.setCanSMILES(nodeDesc);
//				List<Attachment> al = new ArrayList<Attachment>();
//				int start = 0;
//				int pos = nodeDesc.indexOf("R", start);
//				String number = "";
//				while (pos >= 0) {
//					pos++;
//					String letter = nodeDesc.substring(pos, pos + 1);
//					while (letter.matches("\\d")) {
//						number = number + letter;
//						pos++;
//						letter = nodeDesc.substring(pos, pos + 1);
//					}
//
//					try {
//						Attachment tmpAtt = DeepCopy.copy(R1HAtt);
//						tmpAtt.setLabel("R" + number);
//						tmpAtt.setAlternateId("R" + number + "-H");
//						String oldSmi = tmpAtt.getCapGroupSMILES();
//						String newSmi = oldSmi.replace("R1", "R" + number);
//						tmpAtt.setCapGroupSMILES(newSmi);
//						al.add(tmpAtt);
//					} catch (Exception ex) {
//						throw new NotationException(
//								"Unable to create attachment by copying from attachment database",
//								ex);
//					}
//
//					start = pos;
//					pos = nodeDesc.indexOf("R", start);
//					number = "";
//				}
//
//				m.setAttachmentList(al);
//				try {
//					monomerStore.addNewMonomer(m);
//					MonomerFactory.setDBChanged(true);
//				} catch (Exception ex) {
//					throw new NotationException(
//							"Unable to add adhoc new monomer into monomer databse",
//							ex);
//				}
//			}
//		}
//		return alternateId;
//	}	
	private static Map<String, Integer> seedMap = new HashMap<String, Integer>();
	
	public static String getAdHocMonomerIDPrefix(String polymerType) {
		if ( polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
			return "CM#";
		}
		else if ( polymerType.equals(Monomer.PEPTIDE_POLYMER_TYPE)) {
			return "PM#";
		}
		else if ( polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
			return "NM#";
		}
		else {
			return "AM#";
		}
	}

	public static String generateNextAdHocMonomerID(String polymerType,
			MonomerStore store) {
		Map<String, Monomer> internalMonomers = null;
		try {
			internalMonomers = MonomerFactory.getInstance().getMonomerDB().get(polymerType);
		} catch (Exception e) {
			// TODO: handle exception
		}
		
		Map<String, Monomer> monomers = store.getMonomers(polymerType);
		
		Integer seed = seedMap.get(polymerType);
		if ( seed == null) {
			seed = 0;
		}
		seed++;
		seedMap.put( polymerType, seed);
		
		String result = getAdHocMonomerIDPrefix(polymerType) + seed;
		
		
		if ((monomers!=null && monomers.containsKey(result))
				|| (internalMonomers != null && internalMonomers.containsKey(result))) {
			return generateNextAdHocMonomerID(polymerType, store);
		} else {
			return result;
		}
	}


	public static void resetSeed() {
		seedMap = new HashMap<String, Integer>();
	}

	public static int getMatchingBracketPosition(char[] characters,
			int position, char openingBracket, char closingBracket) {
		if (position < (characters.length - 1)
				&& characters[position] == openingBracket) {
			int currentPosition = position;
			int openingBracketCount = 1;

			do {
				char currentCharacter = characters[++currentPosition];
				if (currentCharacter == openingBracket) {
					openingBracketCount++;
				} else if (currentCharacter == closingBracket) {
					openingBracketCount--;
				}
			} while (openingBracketCount > 0
					&& currentPosition < (characters.length - 1));

			if (characters[currentPosition] == closingBracket) {
				return currentPosition;
			} else {
				return -1;
			}
		} else {
			return -1;
		}
	}
	
	
	/*
	 * This function replaces smiles in simple notation with temporary ids
	 * 
	 */
	public static String getNotationByReplacingSmiles(String simpleNotation,String polymerType,MonomerStore store) throws NotationException, MonomerException, JDOMException, IOException{
		
		List<String> monomerIDs = getMonomerIDList(simpleNotation, polymerType,
				store);
		
		return getSimpleNotation(monomerIDs, polymerType, store);
		
	}

}
