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

import chemaxon.formats.MolImporter;
import chemaxon.marvin.calculations.ElementalAnalyserPlugin;
import chemaxon.marvin.plugin.PluginException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;

import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.MoleculeInfo;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

/**
 * This class provides methods that handle chemical structures
 * 
 * @author zhangtianhong
 */
public class StructureParser {

	public static final String SMILES_EXTENSION_SEPARATOR_REGEX = "\\|";
	public static final String SMILES_EXTENSION_SEPARATOR_SYMBOL = "|";
	public static final String EXTENSION_COMPONENT_SEPARATOR_REGEX = "\\$";
	public static final String EXTENSION_COMPONENT_SEPARATOR_SYMBOL = "$";
	public static final String SMILES_MIXTURE_DELIMITER_REGEX = "\\.";
	public static final String ATOM_POSITION_DELIMITER_REGEX = "\\$|;";
	public static final String ATOM_POSITION_DELIMITER_SYMBOL = ";";
	public static final String WHITE_SPACE_REGEX = "\\s";
	// 2014-05-13 -e :do not include enhanced stereochemistry features
	public static final String CHEMAXON_EXTENDEND_SMILES_FORMAT = "cxsmiles:u-e";
	public static final String UNIQUE_SMILES_FORMAT = "smiles:u";
	public static final String SMILES_FORMAT = "smiles";
	private static final int CANONICALIZATION_ROUND = 3;

	public static MoleculeInfo getMoleculeInfo(String smiles)
			throws IOException, PluginException {
		Molecule mol = getMolecule(smiles);
		ElementalAnalyserPlugin plugin = new ElementalAnalyserPlugin();
		plugin.setDoublePrecision(2);
		plugin.setMolecule(mol);
		plugin.run();
		MoleculeInfo mi = new MoleculeInfo();
		mi.setMolecularFormula(plugin.getFormula());
		mi.setMolecularWeight(plugin.getMass());
		mi.setExactMass(plugin.getExactMass());
		return mi;
	}

	/**
	 * This methods validates if the provided SMILES string is valid or not Can
	 * be called by any program interested in SMILES validation
	 * 
	 * @param smiles
	 *            SMILES string to be validated
	 * @return true or false
	 * @throws java.io.IOException
	 */
	public static boolean validateSmiles(String smiles) throws IOException {
		Molecule mol = getMolecule(smiles);
		for (int i = 0; i < mol.getAtomCount(); i++) {
			MolAtom a = mol.getAtom(i);
			a.valenceCheck();
			if (a.hasValenceError()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * This method removes everything in ChemAxon Extended SMILES extension, but
	 * Atom Position Input: [*][H].O[*].OC1=CC=C(C[C@H](N[*])C([*])=O)C=C1
	 * |r,$_R1;;;_R2;;;;;;;;;_R1;;_R2;;;$,c:14,t:3,5| Output:
	 * [*][H].O[*].OC1=CC=C(C[C@H](N[*])C([*])=O)C=C1
	 * |$_R1;;;_R2;;;;;;;;;_R1;;_R2;;;$|
	 * 
	 * @param extendedSMILES
	 * @return simpleExtendedSMILES
	 * @throws org.helm.notation.StructureException
	 */
	public static String getSimpleExtendedSMILES(String extendedSMILES)
			throws StructureException {
		String[] components = extendedSMILES
				.split(SMILES_EXTENSION_SEPARATOR_REGEX);
		String smi = null;
		String ext = null;

		if (components.length == 2) {
			smi = components[0];
			String extension = components[1];

			String[] extComponents = extension.split(
					EXTENSION_COMPONENT_SEPARATOR_REGEX, -1);

			if (extComponents.length == 3) {
				ext = extComponents[1];
			} else {
				throw new StructureException("Invalid extended SMILES format");
			}
		} else {
			throw new StructureException("Invalid extended SMILES format");
		}

		return smi + SMILES_EXTENSION_SEPARATOR_SYMBOL
				+ EXTENSION_COMPONENT_SEPARATOR_SYMBOL + ext
				+ EXTENSION_COMPONENT_SEPARATOR_SYMBOL
				+ SMILES_EXTENSION_SEPARATOR_SYMBOL;
	}

	/**
	 * This method returns unique cxsmiles i.e. cxsmiles:u format, to be used
	 * for monomer registration
	 * 
	 * @param extendedSMILES
	 * @return chemaxon extended unique smiles
	 * @throws StructureException
	 * @throws IOException
	 */
	public static String getUniqueExtendedSMILES(String extendedSMILES)
			throws StructureException, IOException {
		Molecule mol = getMolecule(extendedSMILES);
		return getUniqueExtendedSMILES(mol);
	}

	public static String getUniqueExtendedSMILES(Molecule molecule)
			throws StructureException, IOException {
		return getUniqueSmiles(molecule, CHEMAXON_EXTENDEND_SMILES_FORMAT);
	}

	/**
	 * This method returns unique smiles i.e. smiles:u format, to be used for
	 * compound registration
	 * 
	 * @param smiles
	 * @return unique smiles
	 * @throws IOException
	 */
	public static String getUniqueSmiles(String smiles) throws IOException {
		Molecule mol = StructureParser.getMolecule(smiles);
		return getUniqueSmiles(mol);
	}

	public static String getUniqueSmiles(Molecule molecule) throws IOException {
		return getUniqueSmiles(molecule, UNIQUE_SMILES_FORMAT);
	}

	public static String getUniqueSmiles(String smiles, String format)
			throws IOException {
		Molecule mol = StructureParser.getMolecule(smiles);
		return getUniqueSmiles(mol, format);
	}

	public static String getUniqueSmiles(Molecule molecule, String format)
			throws IOException {
		String oldSmi = "";
		molecule.implicitizeHydrogens(MolAtom.ALL_H);
		String newSmi = molecule.toFormat(format);
		int count = 0;
		while (!oldSmi.equals(newSmi) && count < CANONICALIZATION_ROUND) {
			count++;
			oldSmi = newSmi;

			Molecule mol = StructureParser.getMolecule(oldSmi);
			mol.implicitizeHydrogens(MolAtom.ALL_H);
			newSmi = mol.toFormat(format);
		}

		return newSmi;
	}

	public static String SMILES2ChemAxonPDB(String smiles) throws IOException {
		Molecule mol = getMolecule(smiles);
		mol.addExplicitHydrogens(0);
		mol.clean(3, null);
		return mol.toFormat("pdb");
	}

	/**
	 * convert SMILES or MOLFILE to Molecule
	 * 
	 * @param structureInput
	 *            input SMILES or MOLFILE string
	 * @return Molecule object
	 * @throws java.io.IOException
	 */
	public static Molecule getMolecule(String structureInput) throws IOException {
		Molecule molecule = null;
		if (null != structureInput) {
			InputStream is = new ByteArrayInputStream(structureInput.getBytes());
			MolImporter importer = new MolImporter(is);
			molecule = importer.read();
		}
		return molecule;
	}

	/**
	 * mol is changed after this step, make sure to clone it if need to preserve
	 * state the first R group with the specifid ID is removed
	 * 
	 * @param mol
	 * @param rGroupId
	 * @return MolAtom this is the atom in the new mol which will be used for
	 *         connection later, could be null
	 * @throws org.helm.notation.StructureException
	 */
	public static MolAtom removeRgroup(Molecule mol, int rGroupId)
			throws StructureException {
		MolAtom atom = null;
		MolAtom[] atoms = mol.getAtomArray();
		for (int i = 0; i < atoms.length; i++) {
			int rId = atoms[i].getRgroup();

			// found first R group
			if (rId == rGroupId) {
				atom = removeRgroup(mol, atoms[i]);
				break;
			}
		}
		if (null != atom) {
			return atom;
		} else {
			throw new StructureException(
					"Rgroup does not exist in the structure");
		}
	}

	/**
	 * MolAtom associates with the first R group with the specifid ID is
	 * returned
	 * 
	 * @param mol
	 * @param rGroupId
	 * @return MolAtom the MolAtom object that represents the first R group with
	 *         rGroupId
	 * @throws org.helm.notation.StructureException
	 */
	public static MolAtom getRgroupAtom(Molecule mol, int rGroupId)
			throws StructureException {
		MolAtom atom = null;
		MolAtom[] atoms = mol.getAtomArray();
		for (int i = 0; i < atoms.length; i++) {
			int rId = atoms[i].getRgroup();

			// found first R group
			if (rId == rGroupId) {
				atom = atoms[i];
				break;
			}
		}
		if (null != atom) {
			return atom;
		} else {
			throw new StructureException(
					"Rgroup does not exist in the structure");
		}
	}

	/**
	 * mol is changed after this step, make sure to clone it if need to preserve
	 * state It is possible that there are more than one R groups that have the
	 * same ID as rgroup
	 * 
	 * @param mol
	 * @param rgroup
	 * @return MolAtom which is the single connection MolAtom to rgroup
	 * @throws org.helm.notation.StructureException
	 */
	public static MolAtom removeRgroup(Molecule mol, MolAtom rgroup)
			throws StructureException {
		MolAtom atom = null;

		if (rgroup.getBondCount() == 1) {
			MolBond bond = rgroup.getBond(0);
			if (bond.getAtom1().equals(rgroup)) {
				atom = bond.getAtom2();
			} else {
				atom = bond.getAtom1();
			}
			mol.removeNode(rgroup);
		} else {
			throw new StructureException("Only one bond is allowed from Rgroup");
		}

		return atom;
	}

	/**
	 * convert extendedSMILES for mixture into a list of extendedSMILES for each
	 * component
	 * 
	 * @param mixtureExtendedSMILES
	 * @return List<String> list of individual SMILES string
	 * @throws org.helm.notation.StructureException
	 * @throws java.io.IOException
	 */
	public static List<String> decomposeMixture(String mixtureExtendedSMILES)
			throws StructureException, IOException {

		String mixtureSimpleExtendedSMILES = getSimpleExtendedSMILES(mixtureExtendedSMILES);

		String[] mainAndExtension = mixtureSimpleExtendedSMILES
				.split(SMILES_EXTENSION_SEPARATOR_REGEX);
		String main = mainAndExtension[0].replaceAll(WHITE_SPACE_REGEX, ""); // remove
																				// white
																				// space
		String extension = mainAndExtension[1];

		String[] allSmiles = main.split(SMILES_MIXTURE_DELIMITER_REGEX);
		String[] tokenAtoms = extension
				.split(ATOM_POSITION_DELIMITER_REGEX, -1); // need to remove
															// first and last
															// empty string

		List<String> allAtomList = new ArrayList<String>();
		for (int i = 0; i < tokenAtoms.length; i++) {
			if (i == 0 || i == tokenAtoms.length - 1) {
				if (tokenAtoms[i].length() != 0) {
					allAtomList.add(tokenAtoms[i]);
				}
			} else {
				allAtomList.add(tokenAtoms[i]);
			}

		}
		String[] allAtoms = new String[allAtomList.size()];
		allAtomList.toArray(allAtoms);

		int pos = 0;

		List<String> smilesExtList = new ArrayList<String>();

		for (int i = 0; i < allSmiles.length; i++) {
			String smiles = allSmiles[i];
			Molecule mol = getMolecule(smiles);
			int atomCount = mol.getAtomCount();

			StringBuffer sb = new StringBuffer();
			sb.append(SMILES_EXTENSION_SEPARATOR_SYMBOL
					+ EXTENSION_COMPONENT_SEPARATOR_SYMBOL); // |$
			for (int j = 0; j < atomCount; j++) {
				sb.append(allAtoms[pos]);
				sb.append(ATOM_POSITION_DELIMITER_SYMBOL);
				pos++;
			}
			sb.replace(sb.length() - 1, sb.length(),
					EXTENSION_COMPONENT_SEPARATOR_SYMBOL
							+ SMILES_EXTENSION_SEPARATOR_SYMBOL); // $|
			String ext = sb.toString();
			smilesExtList.add(smiles + " " + ext);
		}

		return smilesExtList;
	}

	/**
	 * This method merges molecule2 into molecule1,
	 * 
	 * @param molecule1
	 * @param molAtom1
	 *            atom to be removed, the connected atom is used for merging
	 * @param molecule2
	 * @param molAtom2
	 *            atom to be removed from molecule2, the connected atom is used
	 *            for merging
	 * @throws org.helm.notation.StructureException
	 */
	public static void mergeIgnoreStereo(Molecule molecule1, MolAtom molAtom1,
			Molecule molecule2, MolAtom molAtom2) throws StructureException {

		// cyclization is allowed
		if (molecule1 == molecule2) {
			molecule1.dearomatize();
			MolAtom atom1 = removeRgroup(molecule1, molAtom1);
			MolAtom atom2 = removeRgroup(molecule1, molAtom2);
			MolBond bond = new MolBond(atom1, atom2);
			molecule1.add(bond);
		} else {
			molecule1.dearomatize();
			molecule2.dearomatize();

			MolAtom atom1 = removeRgroup(molecule1, molAtom1);

			MolAtom atom2 = removeRgroup(molecule2, molAtom2);

			MolAtom[] atoms = molecule2.getAtomArray();
			for (int i = 0; i < atoms.length; i++) {
				molecule1.add(atoms[i]);
			}

			MolBond[] bonds = molecule2.getBondArray();
			for (int i = 0; i < bonds.length; i++) {
				molecule1.add(bonds[i]);
			}

			MolBond bond = new MolBond(atom1, atom2);
			molecule1.add(bond);
		}
	}
        
     /**
     * This method merges molecule2 to molecule1, and keeps the single stereo bond.
     * Old merge method is renamed to mergeIgnoreStereo, which does not preserve single stereo bond connected to R atom
     * @param molecule1 base molecule to merge to
     * @param molAtom1 R atom to be removed, connecting atom will be used for merging
     * @param molecule2 molecule to be merged
     * @param molAtom2 R atom to be removed
     * @throws StructureException 
     */
    public static void merge(Molecule molecule1, MolAtom molAtom1, Molecule molecule2, MolAtom molAtom2) throws StructureException {
        if (isSingleStereo(molAtom1) && isSingleStereo(molAtom2)) {
            throw new StructureException("Both R atoms are connected to chiral centers");
        }

        molecule1.dearomatize();
        molecule2.dearomatize();
        MolBond mergeBond = null;

        if (isSingleStereo(molAtom2)) {
            //keep chiral bond on molAtom2
            MolBond chiralBond = molAtom2.getBond(0);
            MolAtom atom2 = chiralBond.getOtherAtom(molAtom2);
            MolAtom atom1 = removeRgroup(molecule1, molAtom1);
            mergeBond = chiralBond.cloneBond(atom1, atom2);
            molecule1.add(mergeBond);
            removeRgroup(molecule2, molAtom2);

        } else {
            //keep bond on molAtom1, regardless stereo type
            MolBond keepBond = molAtom1.getBond(0);
            MolAtom atom1 = keepBond.getOtherAtom(molAtom1);
            MolAtom atom2 = removeRgroup(molecule2, molAtom2);
            mergeBond = keepBond.cloneBond(atom1, atom2);
            molecule1.add(mergeBond);
            removeRgroup(molecule1, molAtom1);
        }

        if (molecule1 != molecule2) {
            MolAtom[] atoms = molecule2.getAtomArray();
            for (MolAtom atom : atoms) {
                molecule1.add(atom);
            }

            MolBond[] bonds = molecule2.getBondArray();
            for (MolBond bond : bonds) {
                molecule1.add(bond);
            }
        }
    }
    
    /**
     * Check if the bond that R atom is part of is single stereo bond or not
     *
     * @param rAtom the R atom to be checked
     * @return true or false
     * @throws StructureException
     */
    protected static boolean isSingleStereo(MolAtom rAtom) throws StructureException {
        int bondCount = rAtom.getBondCount();
        if (bondCount != 1) {
            throw new StructureException("RGroup is allowed to have single connection to other atom");
        }

        MolBond bond = rAtom.getBond(0);
        int bondType = bond.getFlags() & MolBond.STEREO1_MASK;

        return bondType == MolBond.UP || bondType == MolBond.DOWN || bondType == MolBond.WAVY;
    }



	/**
	 * This method should be used by SimpleNotationParser and
	 * ComplexNotationParser, not exposed to public
	 * 
	 * @param chunks
	 *            - list of main MoleculeInfo
	 * @param caps
	 *            - list of cap MoleculeInfo
	 * @return MoleculeInfo after processing
	 * @throws NotationException
	 */
	protected static MoleculeInfo processMoleculeInfo(
			List<MoleculeInfo> chunks, List<MoleculeInfo> caps)
			throws NotationException {
		MoleculeInfo result = new MoleculeInfo();

		double mw = 0.0;
		double exact = 0.0;
		TreeMap<String, Integer> atomNumberMap = new TreeMap<String, Integer>();

		for (MoleculeInfo mi : chunks) {
			mw = mw + mi.getMolecularWeight();
			exact = exact + mi.getExactMass();
			String formula = mi.getMolecularFormula();
			processMoleculeFormula(atomNumberMap, formula, ADDITION);
		}

		for (MoleculeInfo mi : caps) {
			mw = mw - mi.getMolecularWeight();
			exact = exact - mi.getExactMass();
			String formula = mi.getMolecularFormula();
			processMoleculeFormula(atomNumberMap, formula, SUBTRACTION);
		}

		StringBuilder sb = new StringBuilder();
		Set<String> atoms = atomNumberMap.keySet();
		for (Iterator<String> i = atoms.iterator(); i.hasNext();) {
			String atom = i.next();
			Integer num = atomNumberMap.get(atom);
			sb.append(atom);
			sb.append(num.toString());
		}

		result.setMolecularWeight(mw);
		result.setExactMass(exact);
		result.setMolecularFormula(sb.toString());

		return result;
	}

	private static final int ADDITION = 1;
	private static final int SUBTRACTION = 2;

	private static void processMoleculeFormula(
			TreeMap<String, Integer> atomNumberMap, String formula,
			int operationType) throws NotationException {
		char[] chars = formula.toCharArray();

		String atom = "";
		String number = "";

		for (int i = 0; i < chars.length; i++) {
			String oneChar = String.valueOf(chars[i]);
			if (oneChar.matches("[A-Z]")) {
				if (atom.length() == 0) {
					atom = oneChar;
				} else {
					updateAtomNumberMap(atomNumberMap, atom, number,
							operationType);

					atom = oneChar;
					number = "";
				}
			} else if (oneChar.matches("[a-z]")) {
				if (atom.length() > 0) {
					atom = atom + oneChar;
				}
			} else {
				if (number.length() == 0) {
					number = oneChar;
				} else {
					number = number + oneChar;
				}
			}
		}

		updateAtomNumberMap(atomNumberMap, atom, number, operationType);
	}

	private static void updateAtomNumberMap(
			TreeMap<String, Integer> atomNumberMap, String atom, String number,
			int operationType) throws NotationException {
		if (null == number || number.length() == 0) {
			number = "1";
		}

		if (operationType == ADDITION) {
			if (atomNumberMap.containsKey(atom)) {
				Integer oldI = atomNumberMap.get(atom);
				Integer newI = new Integer(oldI.intValue()
						+ Integer.parseInt(number));
				atomNumberMap.put(atom, newI);
			} else {
				atomNumberMap.put(atom, new Integer(number));
			}
		} else {
			if (atomNumberMap.containsKey(atom)) {
				Integer oldI = atomNumberMap.get(atom);
				Integer newI = new Integer(oldI.intValue()
						- Integer.parseInt(number));
				atomNumberMap.put(atom, newI);
			} else {
				throw new NotationException("Atom " + atom
						+ " exists in cap group but not in main structure");
			}

		}

	}

	public static List<String> getRGroupsFromExtendedSmiles(
			String extendedSmiles) {
		List<String> list = new ArrayList<String>();
		String[] tokens = extendedSmiles.split("R", -1);
		if (tokens.length > 1) {
			for (int i = 1; i < tokens.length; i++) {
				String token = tokens[i];
				char[] chars = token.toCharArray();
				String numbers = "";
				for (int j = 0; j < chars.length; j++) {
					String letter = String.valueOf(chars[j]);
					if (letter.matches("[0-9]")) {
						numbers += letter;
					} else {
						break;
					}
				}

				if (numbers.length() > 0) {
					numbers = "R" + numbers;
					list.add(numbers);
				}
			}
		}

		return list;
	}

}
