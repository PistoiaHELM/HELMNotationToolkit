/**
 * *****************************************************************************
 * Copyright C 2012, The Pistoia Alliance
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ****************************************************************************
 */
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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.jdom.JDOMException;

/**
 * This class provides methods that converter peptide chemical structures into
 * polymer notation handles amide bond and disulfide bond only. 
 * Here are the steps taken:
 * 1. find all amide (C(=O)-N) bond and disulfide (S-S) bond in the structure
 * 2. exclude preserved amide bonds such as urea, urethane and imide
 * 3. fragment the molecule by breaking these bonds and keeping history in a tree model
 * 4. close opened lactam in any fragment
 * 5. match structure of each fragment with monomer core structure 
 * 6. generate extended smiles for fragment without match, and use for in line HELM notation
 * 7. combine broken side chain amine from Asparagine (N) and Glutamine (Q) with Aspartic Acid (D) or Glutamic Acid (E)
 * 8. generate HELM notation based on monomer ID and tree model
 * 
 * Limitations:
 * 1. Only supports R1, R2 and R3 attachments
 * 2. R1, R2, R3 assignment to each fragment could be off (based on alpha amino acid configuration)
 * 
 * @author zhangtianhong
 */
public class PeptideStructureParser {

    private final List<String> IDList = new ArrayList<String>();
    // unique smiles, smiles:u
    private final List<String> nonCappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R1CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R2CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R3CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R1R2CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R1R3CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> R2R3CappedUniqueSmilesList = new ArrayList<String>();
    private final List<String> allCappedUniqueSmilesList = new ArrayList<String>();
    private static PeptideStructureParser instance;
    private int seedID = 1;

    private static final String AMINE_ATTACHEMENT_LABEL = "R1-R3";
    private static final String CARBONYL_ATTACHEMENT_LABEL = "R2-R3";
    private static final String SULFIDE_ATTACCHMENT_LABEL = "R3";

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

    public void initAminoAcidLists() throws MonomerException, IOException, JDOMException, StructureException {
        Map<String, Monomer> idMonomerMap = MonomerFactory.getInstance().getMonomerDB().get(Monomer.PEPTIDE_POLYMER_TYPE);
        Set<String> idSet = idMonomerMap.keySet();
        for (Iterator<String> i = idSet.iterator(); i.hasNext();) {
            String id = i.next();
            IDList.add(id);

            Monomer m = idMonomerMap.get(id);

            String nonCappedUniqueSmiles = getCappedUniqueSmiles(m, new int[0]);
            nonCappedUniqueSmilesList.add(nonCappedUniqueSmiles);

            String r1CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{1});
            R1CappedUniqueSmilesList.add(r1CappedUniqueSmiles);

            String r2CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{2});
            R2CappedUniqueSmilesList.add(r2CappedUniqueSmiles);

            String r3CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{3});
            R3CappedUniqueSmilesList.add(r3CappedUniqueSmiles);

            String r1r2CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{1, 2});
            R1R2CappedUniqueSmilesList.add(r1r2CappedUniqueSmiles);

            String r1r3CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{1, 3});
            R1R3CappedUniqueSmilesList.add(r1r3CappedUniqueSmiles);

            String r2r3CappedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{2, 3});
            R2R3CappedUniqueSmilesList.add(r2r3CappedUniqueSmiles);

            String allCapppedUniqueSmiles = getCappedUniqueSmiles(m, new int[]{1, 2, 3});
            allCappedUniqueSmilesList.add(allCapppedUniqueSmiles);
        }
    }

    public String molfile2notation(String molfile) throws StructureException, IOException {
        Molecule mol = StructureParser.getMolecule(molfile);
        return molecule2notation(mol);
    }

    public String smiles2notation(String smiles) throws StructureException, IOException {
        Molecule mol = StructureParser.getMolecule(smiles);
        return molecule2notation(mol);
    }

    public String molecule2notation(Molecule molecule) throws StructureException, IOException {
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

        closeOpenedLactam(fragTree);

        rollUp(fragTree);

        List<List<String>> monomerIDsList = fragTree.getHead().getMonomerIDsList();
        List<List<Map<String, String>>> connectionMapsList = fragTree.getHead().getConnectionMapsList();

        combineOpendSideChainAmide(monomerIDsList, connectionMapsList);

        String notation = generateNotation(monomerIDsList, connectionMapsList);
        return notation;
    }

    private void combineOpendSideChainAmide(List<List<String>> monomerIDsList, List<List<Map<String, String>>> connectionMapsList) {
        //combine "Asparatic Acid" and side chain "amide" to "Asparagine" which is am-R1 and D-R3 connection at the same level => N
        //combine "Glutamic Acid" and side chain "amide" to "Glutamine" which is am-R1 and E-R3 connection at the same level => Q

        List<Integer> aminePolymerIndices = new ArrayList<Integer>();
        List<Integer> otherPolymerIndices = new ArrayList<Integer>();
        for (int index = 0; index < monomerIDsList.size(); index++) {
            List<String> monomerIDs = monomerIDsList.get(index);
            if (monomerIDs.size() == 1 && monomerIDs.get(0).equals("am")) {
                aminePolymerIndices.add(index);
            } else {
                otherPolymerIndices.add(index);
            }
        }

        List<Integer> removalIndices = new ArrayList<Integer>();
        if (!aminePolymerIndices.isEmpty()) {
            for (Integer index : aminePolymerIndices) {
                List<Map<String, String>> connectionMaps = connectionMapsList.get(index);
                Map<String, String> amineConnections = connectionMaps.get(0);
                String level = amineConnections.keySet().iterator().next();

                for (Integer otherIndex : otherPolymerIndices) {
                    List<String> monomerIDs = monomerIDsList.get(otherIndex);
                    boolean foundMatch = false;
                    for (int monomerIndex = 0; monomerIndex < monomerIDs.size(); monomerIndex++) {
                        String monomerID = monomerIDs.get(monomerIndex);
                        if (monomerID.equals("D") || monomerID.equals("E")) {
                            List<Map<String, String>> otherConnectionMaps = connectionMapsList.get(otherIndex);
                            Map<String, String> otherConnections = otherConnectionMaps.get(monomerIndex);
                            String connectionR = otherConnections.get(level);

                            if (null != connectionR && connectionR.equals("R3")) {
                                //mark amine polymer for removal
                                removalIndices.add(index);

                                //modify other polymer and its connection
                                String newMonomerID = monomerID;
                                if (monomerID.equals("D")) {
                                    newMonomerID = "N";
                                } else if (monomerID.equals("E")) {
                                    newMonomerID = "Q";
                                }
                                monomerIDs.set(monomerIndex, newMonomerID);
                                otherConnections.remove(level);
                                foundMatch = true;
                                break;
                            }
                        }
                    }
                    if (foundMatch) {
                        break;
                    }
                }
            }
        }

        //remove last one first so that the index won't be messed up
        Collections.sort(removalIndices);
        for (int i = removalIndices.size() - 1; i >= 0; i--) {
            int removalIndex = removalIndices.get(i);
            monomerIDsList.remove(removalIndex);
            connectionMapsList.remove(removalIndex);
        }
    }

    private String generateNotation(List<List<String>> monomerIDsList, List<List<Map<String, String>>> connectionMapsList) throws StructureException, IOException {
        List<String> simpleNotations = new ArrayList<String>();
        for (List<String> monomerIDs : monomerIDsList) {
            String simpleNotation = getSimpleNotation(monomerIDs);
            simpleNotations.add(simpleNotation);
        }

        Map<String, List<String>> connectionMap = new TreeMap<String, List<String>>();
        for (int i = 0; i < connectionMapsList.size(); i++) {
            List<Map<String, String>> connectionMaps = connectionMapsList.get(i);
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

    private String getCappedUniqueSmiles(Monomer monomer, int[] rgroupIDs) throws IOException, StructureException {
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

    private void closeOpenedLactam(Molecule molecule) {
        MolAtom[] atoms = molecule.getAtomArray();
        Map<Integer, List<MolAtom>> levelMap = new HashMap<Integer, List<MolAtom>>();
        for (MolAtom atom : atoms) {
            if (atom.getAtno() == MolAtom.RGROUP) {
                int level = atom.getRgroup();
                if (levelMap.containsKey(level)) {
                    levelMap.get(level).add(atom);
                } else {
                    List<MolAtom> rAtoms = new ArrayList<MolAtom>();
                    rAtoms.add(atom);
                    levelMap.put(level, rAtoms);
                }
            }
        }

        for (Integer level : levelMap.keySet()) {
            List<MolAtom> rAtoms = levelMap.get(level);
            if (rAtoms.size() == 2) {
                MolAtom rAtom1 = rAtoms.get(0);
                MolBond bond1 = rAtom1.getBond(0);
                MolAtom atom1 = bond1.getOtherAtom(rAtom1);

                MolAtom rAtom2 = rAtoms.get(1);
                MolBond bond2 = rAtom2.getBond(0);
                MolAtom atom2 = bond2.getOtherAtom(rAtom2);

                MolBond newBond = new MolBond(atom2, atom1);

                molecule.removeNode(rAtom1);
                molecule.removeEdge(bond1);
                molecule.removeNode(rAtom2);
                molecule.removeEdge(bond2);

                molecule.add(newBond);
            }
        }
    }

    private AminoAcidInfo molecule2AminoAcidInfo(Molecule molecule) throws IOException, StructureException {
        //perform matching    
        String canSmi = StructureParser.getUniqueSmiles(molecule);
        String monomerID = null;
        Map<String, String> rgroupMap = new HashMap<String, String>();

        if (nonCappedUniqueSmilesList.contains(canSmi)) {
            int index = nonCappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R1CappedUniqueSmilesList.contains(canSmi)) {
            int index = R1CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R2CappedUniqueSmilesList.contains(canSmi)) {
            int index = R2CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R3CappedUniqueSmilesList.contains(canSmi)) {
            int index = R3CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R1R2CappedUniqueSmilesList.contains(canSmi)) {
            int index = R1R2CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R1R3CappedUniqueSmilesList.contains(canSmi)) {
            int index = R1R3CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (R2R3CappedUniqueSmilesList.contains(canSmi)) {
            int index = R2R3CappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        } else if (allCappedUniqueSmilesList.contains(canSmi)) {
            int index = allCappedUniqueSmilesList.indexOf(canSmi);
            monomerID = IDList.get(index);
        }

        //map level with rgroup
        List<MolAtom> rAtoms = new ArrayList<MolAtom>();
        for (MolAtom atom : molecule.getAtomArray()) {
            if (atom.getAtno() == MolAtom.RGROUP) {
                rAtoms.add(atom);
            }
        }

        for (MolAtom rAtom : rAtoms) {
            String label = rAtom.getExtraLabel();
            if (label.startsWith(SULFIDE_ATTACCHMENT_LABEL)) {
                rgroupMap.put("R" + rAtom.getRgroup(), "R3");
            } else if (label.startsWith(AMINE_ATTACHEMENT_LABEL)) {
                //could be R1 or R3
                MolAtom nAtom = rAtom.getBond(0).getOtherAtom(rAtom);
                if (isAlphaAminoAcidNitrogenAtom(nAtom)) {
                    rgroupMap.put("R" + rAtom.getRgroup(), "R1");
                } else {
                    if (rAtoms.size() == 1) {
                        rgroupMap.put("R" + rAtom.getRgroup(), "R1");
                    } else {
                        rgroupMap.put("R" + rAtom.getRgroup(), "R3");
                    }
                }
            } else {
                //could be R2 or R3
                MolAtom cAtom = rAtom.getBond(0).getOtherAtom(rAtom);
                if (isAlphaAminoAcidCarbonylCarbonAtom(cAtom)) {
                    rgroupMap.put("R" + rAtom.getRgroup(), "R2");
                } else {
                    if (rAtoms.size() == 1) {
                        rgroupMap.put("R" + rAtom.getRgroup(), "R2");
                    } else {
                        rgroupMap.put("R" + rAtom.getRgroup(), "R3");
                    }
                }
            }
        }

        if (null == monomerID) {
            monomerID = StructureParser.getUniqueExtendedSMILES(molecule);
            List<Integer> rNumbers = new ArrayList<Integer>();
            for (String oldR : rgroupMap.keySet()) {
                Integer no = Integer.parseInt(oldR.substring(1));
                rNumbers.add(no);
            }
            Collections.shuffle(rNumbers);

            for (int i = rNumbers.size() - 1; i >= 0; i--) {
                String oldR = "R" + rNumbers.get(i);
                String newR = rgroupMap.get(oldR);
                String newX = newR.replace("R", "X");
                monomerID = monomerID.replace(oldR, newX);
            }
            monomerID = monomerID.replaceAll("X", "R");
        } else {
            if (monomerID.equals("Ggu")) {
                //gama glutamic acid has the same structure as glutamic acid, attachment points R2 and R3 flipped
                //Since R group assignment is based on alpha amino acid, so change Ggu to E
                monomerID = "E";
            }
        }

        return new AminoAcidInfo(monomerID, rgroupMap);
    }

    private boolean isAlphaAminoAcidNitrogenAtom(MolAtom atom) {
        if (atom.getAtno() == 7) {
            int bondCount = atom.getBondCount();

            for (int i = 0; i < bondCount; i++) {
                MolBond bond = atom.getBond(i);
                MolAtom otherAtom = bond.getOtherAtom(atom);

                if (bond.getType() == 1 && otherAtom.getAtno() == 6) {
                    int carbonBondCount = otherAtom.getBondCount();
                    for (int j = 0; j < carbonBondCount; j++) {
                        MolBond carbonBond = otherAtom.getBond(j);
                        MolAtom neighborAtom = carbonBond.getOtherAtom(otherAtom);
                        if (carbonBond.getType() == 1 && neighborAtom.getAtno() == 6 && isCarbonylCarbonAtom(neighborAtom)) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    private boolean isAlphaAminoAcidCarbonylCarbonAtom(MolAtom atom) {
        if (atom.getAtno() == 6) {
            int bondCount = atom.getBondCount();
            boolean foundOxygen = false;
            boolean foundAlphaNitrogen = false;
            for (int i = 0; i < bondCount; i++) {
                MolBond bond = atom.getBond(i);
                MolAtom otherAtom = bond.getOtherAtom(atom);
                if (bond.getType() == 2 && otherAtom.getAtno() == 8) {
                    foundOxygen = true;
                } else if (bond.getType() == 1 && otherAtom.getAtno() == 6) {
                    int carbonBondCount = otherAtom.getBondCount();

                    for (int j = 0; j < carbonBondCount; j++) {
                        MolBond carbonBond = otherAtom.getBond(j);
                        MolAtom neighborAtom = carbonBond.getOtherAtom(otherAtom);
                        if (carbonBond.getType() == 1 && neighborAtom.getAtno() == 7) {
                            foundAlphaNitrogen = true;
                            break;
                        }
                    }
                }
            }
            if (foundOxygen && foundAlphaNitrogen) {
                return true;
            }
        }
        return false;
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
        if (frag.getId() <= 0) {
            frag.setId(seedID++);
        }
        Molecule molecule = frag.getMolecule();
        int currentLevel = frag.getLevel();
        int currentId = frag.getId();
        MolBond disulfideBond = findDisulfideBond(molecule);
        MolBond amideBond = findAmideBond(molecule);
        if (null != disulfideBond || null != amideBond) {
            List<Molecule> mols = null;
            if (null != disulfideBond) {
                mols = breakDisulfideBond(molecule, disulfideBond, currentId);
            } else {
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
        if (frag.getId() <= 0) {
            frag.setId(seedID++);
        }
        Molecule molecule = frag.getMolecule();
        int currentLevel = frag.getLevel();
        int currentId = frag.getId();
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
        if (frag.getId() <= 0) {
            frag.setId(seedID++);
        }
        Molecule molecule = frag.getMolecule();
        int currentLevel = frag.getLevel();
        int currentId = frag.getId();
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

    private void closeOpenedLactam(Tree<PeptideFragment> tree) throws IOException, StructureException {
        //handle lactam

        Collection<Tree<PeptideFragment>> subtrees = tree.getSubTrees();
        if (subtrees.isEmpty()) {
            Molecule mol = tree.getHead().getMolecule();
            closeOpenedLactam(mol);
        } else {
            for (Iterator it = subtrees.iterator(); it.hasNext();) {
                Tree<PeptideFragment> subtree = (Tree<PeptideFragment>) it.next();
                closeOpenedLactam(subtree);
            }
        }
    }

    private void rollUp(Tree<PeptideFragment> tree) throws IOException, StructureException {
        int depth = tree.getDepth();
        int level = 0;
        while (level < depth) {
            singleRollUp(tree);
            level++;
        }
    }

    private void singleRollUp(Tree<PeptideFragment> tree) throws IOException, StructureException {
        PeptideFragment frag = tree.getHead();
        Molecule mol = frag.getMolecule();
        List<List<String>> monomerIDsList = frag.getMonomerIDsList();

        if (null == monomerIDsList || monomerIDsList.isEmpty()) {
            if (null != mol && !mol.isEmpty()) {
                // is leaf, should be able to lookup ID
                AminoAcidInfo mi = molecule2AminoAcidInfo(mol);
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
            } else {
                Collection<Tree<PeptideFragment>> subtrees = tree.getSubTrees();
                if (subtrees.size() == 1) {
                    Tree<PeptideFragment> subtree1 = subtrees.iterator().next();
                    PeptideFragment frag1 = subtree1.getHead();
                    List<List<String>> monomerIDsList1 = frag1.getMonomerIDsList();
                    List<List<Map<String, String>>> connectionMapList1 = frag1.getConnectionMapsList();

                    if (null != monomerIDsList1 && !monomerIDsList1.isEmpty()) {
                        frag.setMonomerIDsList(monomerIDsList1);
                        frag.setConnectionMapList(connectionMapList1);
                    } else {
                        singleRollUp(subtree1);
                    }

                } else {
                    int parentFragmentID = frag.getId();
                    Iterator<Tree<PeptideFragment>> iterator = subtrees.iterator();
                    Tree<PeptideFragment> subtree1 = iterator.next();
                    PeptideFragment frag1 = subtree1.getHead();
                    List<List<String>> monomerIDsList1 = frag1.getMonomerIDsList();
                    List<List<Map<String, String>>> connectionMapList1 = frag1.getConnectionMapsList();

                    Tree<PeptideFragment> subtree2 = iterator.next();
                    PeptideFragment frag2 = subtree2.getHead();
                    List<List<String>> monomerIDsList2 = frag2.getMonomerIDsList();
                    List<List<Map<String, String>>> connectionMapList2 = frag2.getConnectionMapsList();

                    if (null != monomerIDsList1 && !monomerIDsList1.isEmpty()) {
                        if (null != monomerIDsList2 && !monomerIDsList2.isEmpty()) {
                            ConnectionInfo ci1 = getConnectionInfo(parentFragmentID, monomerIDsList1, connectionMapList1);
                            ConnectionInfo ci2 = getConnectionInfo(parentFragmentID, monomerIDsList2, connectionMapList2);
                            List<List<String>> parentMonomerIDsList = new ArrayList<List<String>>();
                            List<List<Map<String, String>>> parentConnectionMapsList = new ArrayList<List<Map<String, String>>>();

                            if (("R2".equals(ci1.getRgroup()) && "R1".equals(ci2.getRgroup()))
                                    || ("R2".equals(ci2.getRgroup()) && "R1".equals(ci1.getRgroup()))) {
                                boolean keepOrder = true;
                                if ("R2".equals(ci2.getRgroup()) && "R1".equals(ci1.getRgroup())) {
                                    keepOrder = false;
                                }

                                int listIndex1 = ci1.getListIndex();
                                List<String> monomerIDs1 = null;
                                List<Map<String, String>> connectionMaps1 = null;
                                for (int i = 0; i < monomerIDsList1.size(); i++) {
                                    if (i == listIndex1) {
                                        monomerIDs1 = monomerIDsList1.get(i);
                                        connectionMaps1 = connectionMapList1.get(i);
                                    } else {
                                        parentMonomerIDsList.add(monomerIDsList1.get(i));
                                        parentConnectionMapsList.add(connectionMapList1.get(i));
                                    }
                                }
                                connectionMaps1.get(ci1.getMonomerIndex()).remove("R" + parentFragmentID);

                                int listIndex2 = ci2.getListIndex();
                                List<String> monomerIDs2 = null;
                                List<Map<String, String>> connectionMaps2 = null;
                                for (int i = 0; i < monomerIDsList2.size(); i++) {
                                    if (i == listIndex2) {
                                        monomerIDs2 = monomerIDsList2.get(i);
                                        connectionMaps2 = connectionMapList2.get(i);
                                    } else {
                                        parentMonomerIDsList.add(monomerIDsList2.get(i));
                                        parentConnectionMapsList.add(connectionMapList2.get(i));
                                    }
                                }
                                connectionMaps2.get(ci2.getMonomerIndex()).remove("R" + parentFragmentID);

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

                                parentConnectionMapsList.addAll(connectionMapList1);
                                parentConnectionMapsList.addAll(connectionMapList2);
                            }
                            frag.setMonomerIDsList(parentMonomerIDsList);
                            frag.setConnectionMapList(parentConnectionMapsList);

                        } else {
                            singleRollUp(subtree2);
                        }
                    } else {
                        singleRollUp(subtree1);
                        if (null == monomerIDsList2 || monomerIDsList2.isEmpty()) {
                            singleRollUp(subtree2);
                        }
                    }
                }
            }
        }
    }

    private ConnectionInfo getConnectionInfo(int parentID, List<List<String>> monomerIDsList, List<List<Map<String, String>>> connectionMapList) throws StructureException {
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

        ConnectionInfo ci = new ConnectionInfo(listIndex, monomerIndex, monomerID, rgroup);

        return ci;
    }

    private List<Molecule> breakDisulfideBond(Molecule mol, MolBond bond, int rGroupId) {
        List<Molecule> l = new ArrayList<Molecule>();

        MolAtom atom1 = bond.getAtom1();
        MolAtom atom2 = bond.getAtom2();

        MolAtom r1Atom = new MolAtom(MolAtom.RGROUP);
        r1Atom.setRgroup(rGroupId);
        r1Atom.setExtraLabel(SULFIDE_ATTACCHMENT_LABEL);
        MolAtom r2Atom = new MolAtom(MolAtom.RGROUP);
        r2Atom.setRgroup(rGroupId);
        r2Atom.setExtraLabel(SULFIDE_ATTACCHMENT_LABEL);

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

    private List<Molecule> breakAmideBond(Molecule mol, MolBond bond, int rGroupId) {
        List<Molecule> l = new ArrayList<Molecule>();

        MolAtom atom1 = bond.getAtom1();
        MolAtom atom2 = bond.getAtom2();

        MolAtom r1Atom = new MolAtom(MolAtom.RGROUP);
        r1Atom.setRgroup(rGroupId);
        r1Atom.setExtraLabel(AMINE_ATTACHEMENT_LABEL);
        MolAtom r2Atom = new MolAtom(MolAtom.RGROUP);
        r2Atom.setRgroup(rGroupId);
        r2Atom.setExtraLabel(CARBONYL_ATTACHEMENT_LABEL);

        MolBond carbonylR2Bond;
        MolBond aminoR1Bond;

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
            boolean found = false;
            if (isCarbonylCarbonAtom(atom1) && isAminoNitrogenAtom(atom2)) {
                found = true;
            } else if  (isCarbonylCarbonAtom(atom2) && isAminoNitrogenAtom(atom1)) {
                found = true;
            }
            
            if (found && !isPreservedAmideBond(bond)) {
                return true;
            }
        }

        return false;
    }

    private boolean isPreservedAmideBond(MolBond amideBond) {
        MolAtom cAtom;
        MolAtom nAtom;
        if (amideBond.getAtom1().getAtno() == 6) {
            cAtom = amideBond.getAtom1();
            nAtom = amideBond.getAtom2();
        } else {
            cAtom = amideBond.getAtom2();
            nAtom = amideBond.getAtom1();
        }

        //urea N-C(=O)-N
        int match = 0;
        for (int i = 0; i < cAtom.getBondCount(); i++) {
            MolBond bond = cAtom.getBond(i);
            MolAtom atom = bond.getOtherAtom(cAtom);
            if (bond.getType() == 2 && atom.getAtno() == 8) {
                match++;
            } else if (bond.getType() == 1 && atom.getAtno() == 7) {
                match++;
            }
        }
        if (match == 3) {
            return true;
        }

        //urethane O-C(=O)-N
        match = 0;
        for (int i = 0; i < cAtom.getBondCount(); i++) {
            MolBond bond = cAtom.getBond(i);
            MolAtom atom = bond.getOtherAtom(cAtom);
            if (bond.getType() == 2 && atom.getAtno() == 8) {
                match++;
            } else if (bond.getType() == 1 && atom.getAtno() == 7) {
                match++;
            } else if (bond.getType() == 1 && atom.getAtno() == 8) {
                match++;
            }
        }
        if (match == 3) {
            return true;
        }

        //imide C(=O)-N-C(=O)
        match = 0;
        for (int i = 0; i < nAtom.getBondCount(); i++) {
            MolBond bond = nAtom.getBond(i);
            MolAtom atom = bond.getOtherAtom(nAtom);
            if (bond.getType() == 1 && atom.getAtno() == 6 && isCarbonylCarbonAtom(atom)) {
                match++;
            }
        }
        if (match == 2) {
            return true;
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

    private boolean isAminoNitrogenAtom(MolAtom atom) {
        if (atom.getAtno() == 7) {
            int bondCount = atom.getBondCount();
            int hcount = atom.getImplicitHcount();

            if (bondCount + hcount == 3) {
                return true;
            }
        }
        return false;
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
