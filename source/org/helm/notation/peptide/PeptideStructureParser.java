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
package org.helm.notation.peptide;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.StructureException;
import org.helm.notation.tools.*;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.jdom.JDOMException;

/**
 * This class provides methods that converter peptide chemical structures into
 * polymer notation Handles amide bond and disulfide bond only
 * 
 * @author zhangtianhong
 */
public class PeptideStructureParser {

	private List<String> IDList = new ArrayList<String>();
	// extended chemaxon smiles cxsmiles, contains R Group position
	private List<String> nonCappedSmilesList = new ArrayList<String>();
	private List<String> R1CappedSmilesList = new ArrayList<String>();
	private List<String> R2CappedSmilesList = new ArrayList<String>();
	private List<String> R3CappedSmilesList = new ArrayList<String>();
	private List<String> R1R2CappedSmilesList = new ArrayList<String>();
	private List<String> R1R3CappedSmilesList = new ArrayList<String>();
	private List<String> R2R3CappedSmilesList = new ArrayList<String>();
	private List<String> allCappedSmilesList = new ArrayList<String>();
	// unique smiles, smiles:u
	private List<String> nonCappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R1CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R2CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R3CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R1R2CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R1R3CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> R2R3CappedUniqueSmilesList = new ArrayList<String>();
	private List<String> allCappedUniqueSmilesList = new ArrayList<String>();
	private static PeptideStructureParser instance;
	private int seedID = 1;

	private PeptideStructureParser() {
	}

	public static PeptideStructureParser getInstance() {
		if (null == instance) {
			instance = new PeptideStructureParser();
		}
		return instance;
	}

	public int getSeedID() {
		return seedID;
	}

	public void setSeedID(int seedID) {
		this.seedID = seedID;
	}

	public void initAminoAcidLists() throws MonomerException, IOException,
			JDOMException, StructureException {
		Map<String, Monomer> idMonomerMap = MonomerFactory.getInstance()
				.getMonomerDB().get(Monomer.PEPTIDE_POLYMER_TYPE);
		Set<String> idSet = idMonomerMap.keySet();
		for (Iterator<String> i = idSet.iterator(); i.hasNext();) {
			String id = i.next();
			IDList.add(id);

			Monomer m = idMonomerMap.get(id);
			// String smiles = m.getCanSMILES();
			// Molecule mol = StructureParser.getMolecule(smiles);
			// mol.setAbsStereo(true);

			String nonCappedSmiles = getCappedSmiles(m, new int[0]);
			String nonCappedUniqueSmiles = getCappedUniqueSmiles(m, new int[0]);
			nonCappedSmilesList.add(nonCappedSmiles);
			nonCappedUniqueSmilesList.add(nonCappedUniqueSmiles);

			String r1CappedSmiles = getCappedSmiles(m, new int[] { 1 });
			String r1CappedUniqueSmiles = getCappedUniqueSmiles(m,
					new int[] { 1 });
			R1CappedSmilesList.add(r1CappedSmiles);
			R1CappedUniqueSmilesList.add(r1CappedUniqueSmiles);

			String r2CappedSmiles = getCappedSmiles(m, new int[] { 2 });
			String r2CappedUniqueSmiles = getCappedUniqueSmiles(m,
					new int[] { 2 });
			R2CappedSmilesList.add(r2CappedSmiles);
			R2CappedUniqueSmilesList.add(r2CappedUniqueSmiles);

			String r3CappedSmiles = getCappedSmiles(m, new int[] { 3 });
			String r3CappedUniqueSmiles = getCappedUniqueSmiles(m,
					new int[] { 3 });
			R3CappedSmilesList.add(r3CappedSmiles);
			R3CappedUniqueSmilesList.add(r3CappedUniqueSmiles);

			String r1r2CappedSmiles = getCappedSmiles(m, new int[] { 1, 2 });
			String r1r2CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[] {
					1, 2 });
			R1R2CappedSmilesList.add(r1r2CappedSmiles);
			R1R2CappedUniqueSmilesList.add(r1r2CappedUniqueSmiles);

			String r1r3CappedSmiles = getCappedSmiles(m, new int[] { 1, 3 });
			String r1r3CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[] {
					1, 3 });
			R1R3CappedSmilesList.add(r1r3CappedSmiles);
			R1R3CappedUniqueSmilesList.add(r1r3CappedUniqueSmiles);

			String r2r3CappedSmiles = getCappedSmiles(m, new int[] { 2, 3 });
			String r2r3CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[] {
					2, 3 });
			R2R3CappedSmilesList.add(r2r3CappedSmiles);
			R2R3CappedUniqueSmilesList.add(r2r3CappedUniqueSmiles);

			String allCappedSmiles = getCappedSmiles(m, new int[] { 1, 2, 3 });
			String allCapppedUniqueSmiles = getCappedUniqueSmiles(m, new int[] {
					1, 2, 3 });
			allCappedSmilesList.add(allCappedSmiles);
			allCappedUniqueSmilesList.add(allCapppedUniqueSmiles);
		}
	}

	public String molfile2notation(String molfile) throws StructureException,
			IOException {
		Molecule mol = StructureParser.getMolecule(molfile);
		return molecule2notation(mol);
	}

	public String smiles2notation(String smiles) throws StructureException,
			IOException {
		Molecule mol = StructureParser.getMolecule(smiles);
		return molecule2notation(mol);
	}

	public String molecule2notation(Molecule molecule)
			throws StructureException, IOException {
		Molecule tmp = molecule.cloneMolecule();
		Molecule[] frags = tmp.convertToFrags();
		if (frags.length > 1) {
			throw new StructureException("Input structure can't be a mixture");
		}

		setSeedID(1);

		PeptideFragment pf = new PeptideFragment();
		pf.setMolecule(molecule);
		pf.setLevel(1);
		Tree<PeptideFragment> fragTree = new Tree<PeptideFragment>(pf);

		breakPeptideBonds(fragTree);

		rollUp(fragTree);

		List<List<String>> monomerIDsList = fragTree.getHead()
				.getMonomerIDsList();
		List<List<Map<String, String>>> connectionMapsList = fragTree.getHead()
				.getConnectionMapsList();

		String notation = generateNotation(monomerIDsList, connectionMapsList);
		return notation;
	}

	private String generateNotation(List<List<String>> monomerIDsList,
			List<List<Map<String, String>>> connectionMapsList)
			throws StructureException, IOException {

		List<String> simpleNotations = new ArrayList<String>();
		for (List<String> monomerIDs : monomerIDsList) {
			String simpleNotation = getSimpleNotation(monomerIDs);
			simpleNotations.add(simpleNotation);
		}

		Map<String, List<String>> connectionMap = new TreeMap<String, List<String>>();
		for (int i = 0; i < connectionMapsList.size(); i++) {
			List<Map<String, String>> connectionMaps = connectionMapsList
					.get(i);
			for (int j = 0; j < connectionMaps.size(); j++) {
				Map<String, String> map = connectionMaps.get(j);
				if (!map.isEmpty()) {
					Set<String> keyset = map.keySet();
					for (Iterator<String> it = keyset.iterator(); it.hasNext();) {
						String id = it.next();
						String rgroup = map.get(id);
						String conString = (i + 1) + ":" + (j + 1) + ":"
								+ rgroup;
						if (connectionMap.containsKey(id)) {
							List<String> l = connectionMap.get(id);
							l.add(conString);
						} else {
							List<String> l = new ArrayList<String>();
							l.add(conString);
							connectionMap.put(id, l);
						}
					}
				}
			}
		}

		List<String> connections = new ArrayList<String>();
		Set<String> keyset = connectionMap.keySet();
		for (Iterator<String> it = keyset.iterator(); it.hasNext();) {
			String id = it.next();
			List<String> l = connectionMap.get(id);
			if (l.size() != 2) {
				throw new StructureException("Connection must be paired");
			}

			String s1 = l.get(0);
			String[] tokens1 = s1.split(":");
			String s2 = l.get(1);
			String[] tokens2 = s2.split(":");

			int rgroupNum1 = Integer.parseInt(tokens1[2].replace("R", ""));
			int rgroupNum2 = Integer.parseInt(tokens2[2].replace("R", ""));

			String connection = null;
			if (rgroupNum1 > rgroupNum2) {
				connection = "PEPTIDE" + tokens1[0] + ",PEPTIDE" + tokens2[0]
						+ "," + tokens1[1] + ":" + tokens1[2] + "-"
						+ tokens2[1] + ":" + tokens2[2];
			} else {
				connection = "PEPTIDE" + tokens2[0] + ",PEPTIDE" + tokens1[0]
						+ "," + tokens2[1] + ":" + tokens2[2] + "-"
						+ tokens1[1] + ":" + tokens1[2];
			}

			connections.add(connection);
		}

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < simpleNotations.size(); i++) {
			if (i > 0) {
				sb.append("|");
			}
			String notation = simpleNotations.get(i);
			sb.append("PEPTIDE");
			sb.append((i + 1));
			sb.append("{");
			sb.append(notation);
			sb.append("}");
		}

		sb.append("$");
		if (!connections.isEmpty()) {
			for (int i = 0; i < connections.size(); i++) {
				if (i > 0) {
					sb.append("|");
				}
				sb.append(connections.get(i));
			}
		}
		sb.append("$$$");

		return sb.toString();
	}

	private String getCappedSmiles(Monomer monomer, int[] rgroupIDs)
			throws IOException, StructureException {
		String molSmi = monomer.getCanSMILES();
		Molecule mol = StructureParser.getMolecule(molSmi);
		for (int rid : rgroupIDs) {
			MolAtom rAtom = null;
			try {
				rAtom = StructureParser.getRgroupAtom(mol, rid);
			} catch (StructureException se) {
			}
			if (null != rAtom) {
				Attachment att = monomer.getAttachment("R" + rid);
				String attSmi = att.getCapGroupSMILES();
				Molecule attMol = StructureParser.getMolecule(attSmi);
				MolAtom attAtom = StructureParser.getRgroupAtom(attMol, rid);
				StructureParser.merge(mol, rAtom, attMol, attAtom);
			}
		}
		mol.setAbsStereo(true);
		return StructureParser.getUniqueExtendedSMILES(mol);
	}

	private String getCappedUniqueSmiles(Monomer monomer, int[] rgroupIDs)
			throws IOException, StructureException {
		String molSmi = monomer.getCanSMILES();
		Molecule mol = StructureParser.getMolecule(molSmi);
		for (int rid : rgroupIDs) {
			MolAtom rAtom = null;
			try {
				rAtom = StructureParser.getRgroupAtom(mol, rid);
			} catch (StructureException se) {
			}
			if (null != rAtom) {
				Attachment att = monomer.getAttachment("R" + rid);
				String attSmi = att.getCapGroupSMILES();
				Molecule attMol = StructureParser.getMolecule(attSmi);
				MolAtom attAtom = StructureParser.getRgroupAtom(attMol, rid);
				StructureParser.merge(mol, rAtom, attMol, attAtom);
			}
		}
		mol.setAbsStereo(true);
		return StructureParser.getUniqueSmiles(mol);
	}

	private AminoAcidInfo getAminoAcidInfo(Molecule molecule)
			throws IOException, StructureException {
		String canSmi = StructureParser.getUniqueSmiles(molecule);

		String cxSmi = StructureParser.getUniqueExtendedSMILES(molecule);
		String monomerID = null;
		Map<String, String> rgroupMap = null;

		if (nonCappedUniqueSmilesList.contains(canSmi)) {
			int index = nonCappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = nonCappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);

		} else if (R1CappedUniqueSmilesList.contains(canSmi)) {
			int index = R1CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R1CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (R2CappedUniqueSmilesList.contains(canSmi)) {
			int index = R2CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R2CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (R3CappedUniqueSmilesList.contains(canSmi)) {
			int index = R3CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R3CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (R1R2CappedUniqueSmilesList.contains(canSmi)) {
			int index = R1R2CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R1R2CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (R1R3CappedUniqueSmilesList.contains(canSmi)) {
			int index = R1R3CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R1R3CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (R2R3CappedUniqueSmilesList.contains(canSmi)) {
			int index = R2R3CappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = R2R3CappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else if (allCappedUniqueSmilesList.contains(canSmi)) {
			int index = allCappedUniqueSmilesList.indexOf(canSmi);
			monomerID = IDList.get(index);
			String monomerCxsmi = allCappedSmilesList.get(index);
			rgroupMap = getRgroupMap(cxSmi, monomerCxsmi);
		} else {
			throw new StructureException("Unknown amino acid structure: "
					+ canSmi);
		}
		return new AminoAcidInfo(monomerID, rgroupMap);
	}

	private Map<String, String> getRgroupMap(String cxsmiles,
			String monomerCxsmiles) throws StructureException, IOException {
		String smi = getCoreSmiles(cxsmiles);
		String mSmi = getCoreSmiles(monomerCxsmiles);
		List<String> rgroups = getRgroups(cxsmiles);
		List<String> mRgroups = getRgroups(monomerCxsmiles);

		Map<String, String> map = new HashMap<String, String>();
		if (smi.equals(mSmi)) {
			if (rgroups.size() != mRgroups.size()) {
				throw new StructureException("Incomparable cxsmiles: "
						+ cxsmiles + " vs. " + monomerCxsmiles);
			}

			for (int i = 0; i < rgroups.size(); i++) {
				String key = rgroups.get(i);
				String value = mRgroups.get(i);
				map.put(key, value);
			}
		} else {
			HashMap<String, String> hashMap = new HashMap<String, String>();
			HashMap<String, String> reverseMap = new HashMap<String, String>();
			TreeMap<Integer, String> treeMap = new TreeMap<Integer, String>();
			for (String r : rgroups) {
				Integer i = new Integer(r.replace("R", ""));
				treeMap.put(i, r);
			}

			int largest = treeMap.lastKey().intValue();
			Set<Integer> keyset = treeMap.keySet();
			Integer[] keys = keyset.toArray(new Integer[0]);
			for (int i = keys.length - 1; i >= 0; i--) {
				largest++;
				Integer key = keys[i];
				String r = treeMap.get(key);
				String tmpR = "R" + largest;
				hashMap.put(r, tmpR);
				reverseMap.put(tmpR, r);
			}

			String newcxsmiles = cxsmiles;
			String rString = getRgroupPostionString(newcxsmiles);
			String[] atomTokens = rString.split(";", -1);
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < atomTokens.length; i++) {
				if (sb.length() > 0) {
					sb.append(";");
				}

				if (atomTokens[i].contains("_R")) {
					String r = atomTokens[i].replace("_", "");
					String tmpR = hashMap.get(r);
					sb.append("_");
					sb.append(tmpR);
				}
			}
			String newRstring = sb.toString();
			newcxsmiles = newcxsmiles.replace("$" + rString + "$", "$"
					+ newRstring + "$");
			newcxsmiles = StructureParser.getUniqueExtendedSMILES(newcxsmiles);

			Map<String, String> newRgroups = getRgroupMap(newcxsmiles,
					monomerCxsmiles);

			Set<String> nset = newRgroups.keySet();
			for (Iterator<String> i = nset.iterator(); i.hasNext();) {
				String newR = i.next();
				String r = reverseMap.get(newR);
				String monomerR = newRgroups.get(newR);
				map.put(r, monomerR);
			}
		}
		return map;
	}

	private String getCoreSmiles(String cxsmiles) {
		String[] tokens = cxsmiles.split(" ", -1);
		if (tokens.length == 2) {
			return tokens[0];
		} else {
			return cxsmiles;
		}
	}

	private List<String> getRgroups(String cxsmiles) {
		List<String> rgroups = new ArrayList<String>();
		String rString = getRgroupPostionString(cxsmiles);

		String[] atomTokens = rString.split(";", -1);

		for (int i = 0; i < atomTokens.length; i++) {
			if (atomTokens[i].contains("_R")) {
				String r = atomTokens[i].replace("_", "");
				rgroups.add(r);
			}
		}
		return rgroups;
	}

	private String getRgroupPostionString(String cxsmiles) {
		String[] tokens = cxsmiles.split("\\$", -1);
		String rString = "";
		if (tokens.length == 3) {
			rString = tokens[1];
		}
		return rString;
	}

	private String getSimpleNotation(List<String> ids) {
		StringBuilder sb = new StringBuilder();
		for (String id : ids) {
			if (sb.length() > 0) {
				sb.append(SimpleNotationParser.GROUP_LEVEL_DELIMITER);
			}

			if (id.length() > 1) {
				sb.append(SimpleNotationParser.MODIFICATION_START_SYMBOL);
			}

			sb.append(id);

			if (id.length() > 1) {
				sb.append(SimpleNotationParser.MODIFICATION_END_SYMBOL);
			}
		}
		return sb.toString();
	}

	public void breakPeptideBonds(Tree<PeptideFragment> tree) {
		PeptideFragment frag = tree.getHead();
		if (frag.getId() <= 0)
			frag.setId(seedID++);
		Molecule molecule = frag.getMolecule();
		int currentLevel = frag.getLevel();
		int currentId = frag.getId();
		// System.out.println("Level: " + currentLevel);
		// System.out.println("Id: " + currentId);
		// System.out.println("cxsmiles: " +
		// molecule.toFormat(CHEMAXON_EXTENDEND_SMILES_FORMAT));
		MolBond disulfideBond = findDisulfideBond(molecule);
		MolBond amideBond = findAmideBond(molecule);
		if (null != disulfideBond || null != amideBond) {
			List<Molecule> mols = null;
			if (null != disulfideBond) {
				// System.out.println("Breaking disulfide bond");
				mols = breakDisulfideBond(molecule, disulfideBond, currentId);
			} else {
				// System.out.println("Breaking amide bond");
				mols = breakAmideBond(molecule, amideBond, currentId);
			}

			int nextLevel = currentLevel + 1;

			if (mols.size() == 1) {
				PeptideFragment pf = new PeptideFragment();
				pf.setMolecule(mols.get(0));
				pf.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree = tree.addLeaf(pf);
				breakPeptideBonds(pfTree);
			} else {
				PeptideFragment pf1 = new PeptideFragment();
				pf1.setMolecule(mols.get(0));
				pf1.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree1 = tree.addLeaf(pf1);
				breakPeptideBonds(pfTree1);

				PeptideFragment pf2 = new PeptideFragment();
				pf2.setMolecule(mols.get(1));
				pf2.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree2 = tree.addLeaf(pf2);
				breakPeptideBonds(pfTree2);
			}
		}
	}

	public void breakDisulfideBonds(Tree<PeptideFragment> tree) {
		PeptideFragment frag = tree.getHead();
		if (frag.getId() <= 0)
			frag.setId(seedID++);
		Molecule molecule = frag.getMolecule();
		int currentLevel = frag.getLevel();
		int currentId = frag.getId();
		// System.out.println("Level: " + currentLevel);
		// System.out.println("Id: " + currentId);
		// System.out.println("cxsmiles: " +
		// molecule.toFormat(CHEMAXON_EXTENDEND_SMILES_FORMAT));
		MolBond disulfideBond = findDisulfideBond(molecule);
		if (null != disulfideBond) {

			List<Molecule> mols = breakDisulfideBond(molecule, disulfideBond,
					currentId);
			int nextLevel = currentLevel + 1;

			if (mols.size() == 1) {
				PeptideFragment pf = new PeptideFragment();
				pf.setMolecule(mols.get(0));
				pf.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree = tree.addLeaf(pf);
				breakDisulfideBonds(pfTree);
			} else {
				PeptideFragment pf1 = new PeptideFragment();
				pf1.setMolecule(mols.get(0));
				pf1.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree1 = tree.addLeaf(pf1);
				breakDisulfideBonds(pfTree1);

				PeptideFragment pf2 = new PeptideFragment();
				pf2.setMolecule(mols.get(1));
				pf2.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree2 = tree.addLeaf(pf2);
				breakDisulfideBonds(pfTree2);
			}
		}
	}

	public void breakAmideBonds(Tree<PeptideFragment> tree) {
		PeptideFragment frag = tree.getHead();
		if (frag.getId() <= 0)
			frag.setId(seedID++);
		Molecule molecule = frag.getMolecule();
		int currentLevel = frag.getLevel();
		int currentId = frag.getId();
		// System.out.println("Level: " + currentLevel);
		// System.out.println("Id: " + currentId);
		// System.out.println("cxsmiles: " +
		// molecule.toFormat(CHEMAXON_EXTENDEND_SMILES_FORMAT));
		MolBond amideBond = findAmideBond(molecule);
		if (null != amideBond) {

			List<Molecule> mols = breakAmideBond(molecule, amideBond, currentId);
			int nextLevel = currentLevel + 1;

			if (mols.size() == 1) {
				PeptideFragment pf = new PeptideFragment();
				pf.setMolecule(mols.get(0));
				pf.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree = tree.addLeaf(pf);
				breakAmideBonds(pfTree);
			} else {
				PeptideFragment pf1 = new PeptideFragment();
				pf1.setMolecule(mols.get(0));
				pf1.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree1 = tree.addLeaf(pf1);
				breakAmideBonds(pfTree1);

				PeptideFragment pf2 = new PeptideFragment();
				pf2.setMolecule(mols.get(1));
				pf2.setLevel(nextLevel);
				Tree<PeptideFragment> pfTree2 = tree.addLeaf(pf2);
				breakAmideBonds(pfTree2);
			}
		}
	}

	public void rollUp(Tree<PeptideFragment> tree) throws IOException,
			StructureException {
		int depth = tree.getDepth();
		int level = 0;
		while (level < depth) {
			singleRollUp(tree);
			level++;
		}
	}

	public void singleRollUp(Tree<PeptideFragment> tree) throws IOException,
			StructureException {
		PeptideFragment frag = tree.getHead();
		Molecule mol = frag.getMolecule();
		List<List<String>> monomerIDsList = frag.getMonomerIDsList();

		if (null == monomerIDsList || monomerIDsList.isEmpty()) {
			if (null != mol && !mol.isEmpty()) {
				generateAminoAcidIDInfo(frag);
			} else {
				Collection<Tree<PeptideFragment>> subtrees = tree.getSubTrees();
				if (subtrees.size() == 1) {
					Tree<PeptideFragment> subtree1 = subtrees.iterator().next();
					PeptideFragment frag1 = subtree1.getHead();
					List<List<String>> monomerIDsList1 = frag1
							.getMonomerIDsList();
					List<List<Map<String, String>>> connectionMapList1 = frag1
							.getConnectionMapsList();

					if (null != monomerIDsList1 && !monomerIDsList1.isEmpty()) {
						frag.setMonomerIDsList(monomerIDsList1);
						frag.setConnectionMapList(connectionMapList1);

						// System.out.println("ID: " + frag.getId());
						// System.out.println("Level: " + frag.getLevel());
						// System.out.println("Monomer IDs: " +
						// frag.getMonomerIDsList().size());
						// System.out.println("Connection Maps: " +
						// frag.getConnectionMapsList().size());
					} else {
						singleRollUp(subtree1);
					}

				} else {
					int parentFragmentID = frag.getId();
					Iterator<Tree<PeptideFragment>> iterator = subtrees
							.iterator();
					Tree<PeptideFragment> subtree1 = iterator.next();
					PeptideFragment frag1 = subtree1.getHead();
					List<List<String>> monomerIDsList1 = frag1
							.getMonomerIDsList();
					List<List<Map<String, String>>> connectionMapList1 = frag1
							.getConnectionMapsList();

					Tree<PeptideFragment> subtree2 = iterator.next();
					PeptideFragment frag2 = subtree2.getHead();
					List<List<String>> monomerIDsList2 = frag2
							.getMonomerIDsList();
					List<List<Map<String, String>>> connectionMapList2 = frag2
							.getConnectionMapsList();

					if (null != monomerIDsList1 && !monomerIDsList1.isEmpty()) {
						if (null != monomerIDsList2
								&& !monomerIDsList2.isEmpty()) {
							ConnectionInfo ci1 = getConnectionInfo(
									parentFragmentID, monomerIDsList1,
									connectionMapList1);
							ConnectionInfo ci2 = getConnectionInfo(
									parentFragmentID, monomerIDsList2,
									connectionMapList2);
							List<List<String>> parentMonomerIDsList = new ArrayList<List<String>>();
							List<List<Map<String, String>>> parentConnectionMapsList = new ArrayList<List<Map<String, String>>>();

							if (("R2".equals(ci1.getRgroup()) && "R1"
									.equals(ci2.getRgroup()))
									|| ("R2".equals(ci2.getRgroup()) && "R1"
											.equals(ci1.getRgroup()))) {
								boolean keepOrder = true;
								if ("R2".equals(ci2.getRgroup())
										&& "R1".equals(ci1.getRgroup())) {
									keepOrder = false;
								}

								int listIndex1 = ci1.getListIndex();
								List<String> monomerIDs1 = null;
								List<Map<String, String>> connectionMaps1 = null;
								for (int i = 0; i < monomerIDsList1.size(); i++) {
									if (i == listIndex1) {
										monomerIDs1 = monomerIDsList1.get(i);
										connectionMaps1 = connectionMapList1
												.get(i);
									} else {
										parentMonomerIDsList
												.add(monomerIDsList1.get(i));
										parentConnectionMapsList
												.add(connectionMapList1.get(i));
									}
								}
								connectionMaps1.get(ci1.getMonomerIndex())
										.remove("R" + parentFragmentID);

								int listIndex2 = ci2.getListIndex();
								List<String> monomerIDs2 = null;
								List<Map<String, String>> connectionMaps2 = null;
								for (int i = 0; i < monomerIDsList2.size(); i++) {
									if (i == listIndex2) {
										monomerIDs2 = monomerIDsList2.get(i);
										connectionMaps2 = connectionMapList2
												.get(i);
									} else {
										parentMonomerIDsList
												.add(monomerIDsList2.get(i));
										parentConnectionMapsList
												.add(connectionMapList2.get(i));
									}
								}
								connectionMaps2.get(ci2.getMonomerIndex())
										.remove("R" + parentFragmentID);

								List<String> newMonomerIDs = new ArrayList<String>();
								if (keepOrder) {
									newMonomerIDs.addAll(monomerIDs1);
									newMonomerIDs.addAll(monomerIDs2);
								} else {
									newMonomerIDs.addAll(monomerIDs2);
									newMonomerIDs.addAll(monomerIDs1);
								}
								parentMonomerIDsList.add(newMonomerIDs);

								List<Map<String, String>> newConnectionMaps = new ArrayList<Map<String, String>>();
								if (keepOrder) {
									newConnectionMaps.addAll(connectionMaps1);
									newConnectionMaps.addAll(connectionMaps2);
								} else {
									newConnectionMaps.addAll(connectionMaps2);
									newConnectionMaps.addAll(connectionMaps1);
								}
								parentConnectionMapsList.add(newConnectionMaps);

							} else {
								parentMonomerIDsList.addAll(monomerIDsList1);
								parentMonomerIDsList.addAll(monomerIDsList2);

								parentConnectionMapsList
										.addAll(connectionMapList1);
								parentConnectionMapsList
										.addAll(connectionMapList2);
							}
							frag.setMonomerIDsList(parentMonomerIDsList);
							frag.setConnectionMapList(parentConnectionMapsList);
							// System.out.println("ID: " + frag.getId());
							// System.out.println("Level: " + frag.getLevel());
							// System.out.println("Monomer IDs: " +
							// frag.getMonomerIDsList().size());
							// System.out.println("Connection Maps: " +
							// frag.getConnectionMapsList().size());

						} else {
							singleRollUp(subtree2);
						}
					} else {
						singleRollUp(subtree1);
						if (null == monomerIDsList2
								|| monomerIDsList2.isEmpty()) {
							singleRollUp(subtree2);
						}
					}
				}
			}
		}
	}

	private ConnectionInfo getConnectionInfo(int parentID,
			List<List<String>> monomerIDsList,
			List<List<Map<String, String>>> connectionMapList)
			throws StructureException {
		String connectionKey = "R" + parentID;
		int listIndex = -1;
		int monomerIndex = -1;
		String rgroup = "";
		String monomerID = "";

		for (int i = 0; i < connectionMapList.size(); i++) {
			List<Map<String, String>> lom = connectionMapList.get(i);
			for (int j = 0; j < lom.size(); j++) {
				Map<String, String> map = lom.get(j);
				if (map.containsKey(connectionKey)) {
					listIndex = i;
					monomerIndex = j;
					rgroup = map.get(connectionKey);
					break;
				}
			}
		}

		if (listIndex < 0 || monomerIndex < 0 || rgroup.length() == 0) {
			throw new StructureException("Unable to find connection info");
		} else {
			List<String> ids = monomerIDsList.get(listIndex);
			monomerID = ids.get(monomerIndex);
		}

		ConnectionInfo ci = new ConnectionInfo(listIndex, monomerIndex,
				monomerID, rgroup);

		return ci;
	}

	private void generateAminoAcidIDInfo(PeptideFragment frag)
			throws IOException, StructureException {
		Molecule mol = frag.getMolecule();
		List<List<String>> monomerIDsList = frag.getMonomerIDsList();
		List<List<Map<String, String>>> connectionMapsList = frag
				.getConnectionMapsList();
		if (!mol.isEmpty() && null == monomerIDsList
				&& null == connectionMapsList) {
			// is leaf, should be able to lookup ID
			AminoAcidInfo mi = getAminoAcidInfo(mol);
			String monomerID = mi.getMonomerID();
			Map<String, String> rgroupMap = mi.getRgroupMap();

			List<String> aaIDs = new ArrayList<String>();
			aaIDs.add(monomerID);
			List<List<String>> newMonomerIDsList = new ArrayList<List<String>>();
			newMonomerIDsList.add(aaIDs);
			frag.setMonomerIDsList(newMonomerIDsList);

			List<Map<String, String>> rgroupMapList = new ArrayList<Map<String, String>>();
			rgroupMapList.add(rgroupMap);
			List<List<Map<String, String>>> newConnectionMapList = new ArrayList<List<Map<String, String>>>();
			newConnectionMapList.add(rgroupMapList);
			frag.setConnectionMapList(newConnectionMapList);

			// System.out.println("ID: " + frag.getId());
			// System.out.println("Level: " + frag.getLevel());
			// System.out.println("Molecule: " +
			// frag.getMolecule().toFormat(CHEMAXON_EXTENDEND_SMILES_FORMAT));
			// System.out.println("Monomer IDs: " +
			// frag.getMonomerIDsList().size());
			// System.out.println("Connection Maps: " +
			// frag.getConnectionMapsList().size());
		}
	}

	private List<Molecule> breakDisulfideBond(Molecule mol, MolBond bond,
			int rGroupId) {
		List<Molecule> l = new ArrayList<Molecule>();

		MolAtom atom1 = bond.getAtom1();
		MolAtom atom2 = bond.getAtom2();

		MolAtom r1Atom = new MolAtom(MolAtom.RGROUP);
		r1Atom.setRgroup(rGroupId);
		MolAtom r2Atom = new MolAtom(MolAtom.RGROUP);
		r2Atom.setRgroup(rGroupId);

		MolBond bond1 = new MolBond(atom1, r1Atom);
		MolBond bond2 = new MolBond(atom2, r2Atom);

		mol.removeEdge(bond);
		mol.add(r1Atom);
		mol.add(bond1);
		mol.add(r2Atom);
		mol.add(bond2);

		Molecule[] fragments = mol.convertToFrags();

		if (fragments.length == 1) {
			l.add(fragments[0]);
		} else {
			l.add(fragments[0]);
			l.add(fragments[1]);
		}
		return l;
	}

	private List<Molecule> breakAmideBond(Molecule mol, MolBond bond,
			int rGroupId) {
		List<Molecule> l = new ArrayList<Molecule>();

		MolAtom atom1 = bond.getAtom1();
		MolAtom atom2 = bond.getAtom2();

		MolAtom r1Atom = new MolAtom(MolAtom.RGROUP);
		r1Atom.setRgroup(rGroupId);
		MolAtom r2Atom = new MolAtom(MolAtom.RGROUP);
		r2Atom.setRgroup(rGroupId);

		MolBond carbonylR2Bond = null;
		MolBond aminoR1Bond = null;

		if (isCarbonylCarbonAtom(atom1)) {
			carbonylR2Bond = new MolBond(atom1, r2Atom);
			aminoR1Bond = new MolBond(atom2, r1Atom);
		} else {
			carbonylR2Bond = new MolBond(atom2, r2Atom);
			aminoR1Bond = new MolBond(atom1, r1Atom);
		}
		mol.removeEdge(bond);
		mol.add(r2Atom);
		mol.add(carbonylR2Bond);

		mol.add(r1Atom);
		mol.add(aminoR1Bond);

		Molecule[] fragments = mol.convertToFrags();

		if (fragments.length == 1) {
			l.add(fragments[0]);
		} else {
			if (fragments[0].contains(carbonylR2Bond)) {
				l.add(fragments[0]);
				l.add(fragments[1]);
			} else {
				l.add(fragments[1]);
				l.add(fragments[0]);
			}
		}
		return l;
	}

	private MolBond findAmideBond(Molecule mol) {
		MolBond[] bonds = mol.getBondArray();

		for (MolBond bond : bonds) {
			if (isAmideBond(bond)) {
				return bond;
			}
		}
		return null;
	}

	private boolean isAmideBond(MolBond bond) {
		MolAtom atom1 = bond.getAtom1();
		MolAtom atom2 = bond.getAtom2();

		if (bond.getType() == 1) {
			if (isCarbonylCarbonAtom(atom1) && isAminoNitrogenAtom(atom2)) {
				return true;
			}

			if (isCarbonylCarbonAtom(atom2) && isAminoNitrogenAtom(atom1)) {
				return true;
			}
		}

		return false;

	}

	private boolean isCarbonylCarbonAtom(MolAtom atom) {
		if (atom.getAtno() == 6) {
			int bondCount = atom.getBondCount();
			int hcount = atom.getImplicitHcount();

			if (bondCount + hcount == 3) {
				for (int i = 0; i < bondCount; i++) {
					MolBond bond = atom.getBond(i);
					MolAtom otherAtom = bond.getOtherAtom(atom);
					if (bond.getType() == 2 && otherAtom.getAtno() == 8) {
						return true;
					}
				}
			}
		}
		return false;
	}

	// ignore terminal amine, which has two implicit hydrogens
	private boolean isAminoNitrogenAtom(MolAtom atom) {
		if (atom.getAtno() == 7) {
			int bondCount = atom.getBondCount();
			int hcount = atom.getImplicitHcount();
			if (bondCount + hcount == 3 && hcount < 2) {
				return true;
			}
		}
		return false;
	}

	private MolBond findDisulfideBond(List<Molecule> molecules) {
		for (Molecule mol : molecules) {
			MolBond bond = findDisulfideBond(mol);
			if (null != bond) {
				return bond;
			}
		}
		return null;

	}

	private MolBond findDisulfideBond(Molecule mol) {
		MolBond[] bonds = mol.getBondArray();

		for (MolBond bond : bonds) {
			if (isDisulfideBond(bond)) {
				return bond;
			}
		}
		return null;
	}

	private boolean isDisulfideBond(MolBond bond) {
		MolAtom atom1 = bond.getAtom1();
		MolAtom atom2 = bond.getAtom2();
		if (bond.getType() == 1) {
			if (isDisulfideSulfurAtom(atom1) && isDisulfideSulfurAtom(atom2)) {
				return true;
			}
		}
		return false;
	}

	// no -SH, has to be R-S-S-R
	private boolean isDisulfideSulfurAtom(MolAtom atom) {
		if (atom.getAtno() == 16) {
			int bondCount = atom.getBondCount();
			if (bondCount == 2) {
				return true;
			}
		}
		return false;
	}
}
