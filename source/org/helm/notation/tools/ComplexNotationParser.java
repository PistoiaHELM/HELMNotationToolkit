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
import org.helm.notation.StructureException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.ComplexPolymer;
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import org.helm.notation.model.RgroupStructure;
import org.helm.notation.model.PolymerEdge;
import org.helm.notation.model.PolymerNode;
import org.helm.notation.model.RNAPolymerNode;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Set;
import java.util.TreeMap;

import org.jdom.JDOMException;

/**
 * 10/4/2011
 * @author zhangtianhong
 * Added support for generic edges, included in valiateComplexNotation(), getCanonicalNotation(), getMoleculeInfo() and getComplexPolymerSMILES()
 * 
 */
/**
 * This class provides methods that handle complex polymer extendendNotation.
 * Complex Polymer Notation is based on graph theory, and uses simple polymer
 * notation and chemical structure as building blocks. The following is a
 * extendendNotation for complex polymer RNA1{R(A)P.[mR](U)P.R([BiA])P.R(C)P}
 * |RNA2{R([BiA])P.R(U)P.R(G)P.R(C)P} |CHEM1{bLink} $RNA1,CHEM1,8:R2-1:R1
 * |CHEM1,RNA2,1:R2-2:R2 Complex notation will follow the following format
 * starting 7/14/2008 NodeList$EdgeList$BasePairList$NodeLabelList$OtherTBD,
 * where EdgeList, BasePairList, NodeLabelList and OtherTBD could be empty, but
 * still keep the positional $ 1. Node List 2. Edge List 3. Base Pair List 4.
 * Node Label List 5. OtherTBD Be careful that monomer XML may contain $, so we
 * have to walk the notation from left to right
 * 
 * @author ZHANGTIANHONG
 */
public class ComplexNotationParser {

	public static final String TOP_LEVEL_STOPPER = "}$";
	public static final String TOP_LEVEL_DELIMITER = "$";
	public static final String LIST_LEVEL_DELIMITER = "|";
	public static final String LIST_LEVEL_DELIMITER_REGEX = "\\|";
	public static final String NODE_LABEL_START_SYMBOL = "{";
	public static final String NODE_LABEL_END_SYMBOL = "}";
	public static final String INVALID_POLYMER_NODE = "Polymer node notation must be in the format of PolymerTypeNumber{Notation}";
	public static final String INVALID_NODE_ID = "Polymer node ID must be in the format of PolymerTypeNumber";
	public static final String DEFAULT_PADDING_CHAR = " ";
	public static final String DEFAULT_BASE_PAIR_CHAR = "|";

	/**
	 * This function checks the given MonomerStore. It is returned unchanged,
	 * when it is neither null nor empty. Else a fresh MonomerStore is fetched
	 * from the MonomerFactory. The function is used to ensure a default value.
	 * 
	 * @param monomerStore
	 * @return monomerStore or default store
	 */
	private static MonomerStore checkForMonomerStore(MonomerStore monomerStore) {
		// check for default value
		MonomerStore combinedMonomerStore = monomerStore;
		// SM 2014-03-03 unit test xHelmNotationParserTest.testXHelmValidation
		// fails
		// do we need to check for empty monomer store?
		// if an invalid xhelm format (empty monomer section) is parsed, the
		// validation is not working anymore because the local monomerStore is
		// used
		if (combinedMonomerStore == null) {// ||
											// combinedMonomerStore.isMonomerStoreEmpty())
											// {
			MonomerFactory factory = null;
			try {
				factory = MonomerFactory.getInstance();
			} catch (Exception ex) {
				Logger.getLogger(ComplexNotationParser.class.getName()).log(
						Level.SEVERE, null, ex);
				return null;
			}
			combinedMonomerStore = factory.getMonomerStore();
		}

		return combinedMonomerStore;
	}

	/**
	 * This methods returns the unique SMILES string for complex polymer
	 * extendendNotation, if all monomers have specific structures
	 * 
	 * @param extendedNotation
	 *            text string for complex polymer extendendNotation
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static String getComplexPolymerSMILES(String extendedNotation)
			throws NotationException, IOException, MonomerException,
			StructureException, JDOMException {

		return getComplexPolymerSMILES(extendedNotation, null);
	}

	public static String getComplexPolymerSMILES(String extendedNotation,
			MonomerStore monomerStore) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		List<Molecule> list = getComplexPolymerStructure(extendedNotation,
				monomerStoreToUse);
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < list.size(); i++) {
			Molecule m = list.get(i);
			String smi = m.toFormat("smiles");
			if (sb.length() > 0) {
				sb.append(".");
			}
			sb.append(smi);
		}

		String mixtureSmiles = sb.toString();
		Molecule mol = StructureParser.getMolecule(mixtureSmiles);

		return mol.toFormat("smiles:u");
	}

	public static String getComplexPolymerCanonicalSmiles(
			String extendedNotation) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {

		return getComplexPolymerCanonicalSmiles(extendedNotation, null);
	}

	public static String getComplexPolymerCanonicalSmiles(
			String extendedNotation, MonomerStore monomerStore)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		List<Molecule> list = getComplexPolymerStructure(extendedNotation,
				monomerStoreToUse);
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < list.size(); i++) {
			Molecule m = list.get(i);
			String smi = StructureParser.getUniqueSmiles(m);
			if (sb.length() > 0) {
				sb.append(".");
			}
			sb.append(smi);
		}

		String mixtureSmiles = sb.toString();
		return StructureParser.getUniqueSmiles(mixtureSmiles);
	}

	/**
	 * This method determines if the extendendNotation contains generic
	 * structure (cannot generate SMILES)
	 * 
	 * @param extendendNotation
	 *            the complete structure extendendNotation
	 * @return true or false
	 */
	public static boolean containsGenericStructure(String extendendNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException {
		String allNodeString = getAllNodeString(extendendNotation);
		List<PolymerNode> nodeList = getPolymerNodeList(allNodeString);
		for (int i = 0; i < nodeList.size(); i++) {
			PolymerNode node = nodeList.get(i);
			if (node.getId().startsWith(Monomer.CHEMICAL_POLYMER_TYPE)) {
				String monomerID = node.getLabel();
				Monomer monomer = SimpleNotationParser.getMonomer(monomerID,
						Monomer.CHEMICAL_POLYMER_TYPE);

				if (null == monomer.getCanSMILES()
						|| monomer.getCanSMILES().length() == 0) {
					return true;
				}

				if (monomer.containAnyAtom()) {
					return true;
				}
			}
		}

		return false;
	}

	/**
	 * This method returns list of Molecule with specific structure, all R
	 * groups are filled
	 * 
	 * @param extendedNotation
	 * @return list of RgroupStructure that are not connected
	 * @throws java.io.IOException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static List<Molecule> getComplexPolymerStructure(
			String extendedNotation) throws IOException, NotationException,
			MonomerException, StructureException, JDOMException {

		return getComplexPolymerStructure(extendedNotation, null);
	}

	public static List<Molecule> getComplexPolymerStructure(
			String extendedNotation, MonomerStore monomerStore)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		int totalMonomerCount = getTotalMonomerCount(extendedNotation,
				monomerStoreToUse);
		if (totalMonomerCount > NotationConstant.MONOMER_COUNT_THRESHOLD) {
			throw new NotationException("Total monomer count ["
					+ totalMonomerCount + "] is above support threshold ["
					+ NotationConstant.MONOMER_COUNT_THRESHOLD + "]");
		}

		ComplexPolymer complexPolymer = parse(extendedNotation, monomerStoreToUse);
		validateComplexPolymer(complexPolymer, monomerStoreToUse);

		List<PolymerNode> nodeList = complexPolymer.getPolymerNodeList();
		List<PolymerEdge> edgeList = complexPolymer.getPolymerEdgeList();

		if (null == edgeList) {
			edgeList = new ArrayList<PolymerEdge>();
		} else {
			for (int i = 0; i < edgeList.size(); i++) {
				PolymerEdge edge = edgeList.get(i);
				if (edge.getEdgeType() != PolymerEdge.STANDARD_EDGE) {
					throw new NotationException(
							"Polymer notation contains non-standard connection");
				}
			}
		}

		// convert each polymer node into RgroupStructure
		// add original poymer ID into rgroup Atom key
		Map<String, RgroupStructure> nodeStrucMap = getPolymerNodeStructureMap(
				nodeList, monomerStoreToUse);
		for (String nodeID : nodeStrucMap.keySet()) {
			RgroupStructure struc = nodeStrucMap.get(nodeID);
			Map<String, MolAtom> rgMap = struc.getRgroupMap();
			Map<String, MolAtom> newRgMap = new HashMap<String, MolAtom>();
			Set keySet = rgMap.keySet();
			for (Iterator itr = keySet.iterator(); itr.hasNext();) {
				String key = (String) itr.next();
				newRgMap.put(nodeID + ":" + key, rgMap.get(key));
			}
			struc.setRgroupMap(newRgMap);
		}

		// list string could be duplicating for self-connection and multiple
		// connections
		// single element list indicates standalone polymer node
		Map<Integer, List<String>> groupNodeListMap = getGroupNodeListMap(
				nodeList, edgeList);

		// convert each connected polymer node group into RgroupStructure
		Map<Integer, RgroupStructure> groupStructureMap = new HashMap<Integer, RgroupStructure>();
		Set groupSet = groupNodeListMap.keySet();
		for (Iterator it = groupSet.iterator(); it.hasNext();) {
			Integer group = (Integer) it.next();
			List<String> nList = groupNodeListMap.get(group);
			if (nList.size() == 1) {
				RgroupStructure struc = nodeStrucMap.get(nList.get(0));
				groupStructureMap.put(group, struc);
			} else {
				List<PolymerEdge> connEdgeList = getConnectedEdgeList(nList,
						edgeList);

				Map<List<String>, RgroupStructure> mergedNodeStrucMap = new HashMap<List<String>, RgroupStructure>();

				// reduce the number of edges to 0 by merging the two nodes into
				// a new node, and updating mergedNodeStrucMap
				// final merged node structure map size is one since all nodes
				// are connected togather
				while (connEdgeList.size() > 0) {
					PolymerEdge edge = connEdgeList.get(0);
					String sourceId = edge.getSourceNode();
					String targetId = edge.getTargetNode();
					List<String> sourceIDList = new ArrayList<String>();
					RgroupStructure sourceStruc;
					Entry<List<String>, RgroupStructure> sourceEntry = getStructureEntry(
							sourceId, mergedNodeStrucMap);
					if (null == sourceEntry) {
						sourceStruc = nodeStrucMap.get(sourceId);
					} else {
						sourceIDList = sourceEntry.getKey();
						sourceStruc = sourceEntry.getValue();
					}

					List<String> targetIDList = new ArrayList<String>();
					RgroupStructure targetStruc;
					Entry<List<String>, RgroupStructure> targetEntry = getStructureEntry(
							targetId, mergedNodeStrucMap);
					if (null == targetEntry) {
						targetStruc = nodeStrucMap.get(targetId);
					} else {
						targetIDList = targetEntry.getKey();
						targetStruc = targetEntry.getValue();
					}

					String sourceR = edge.getSourceUID(); // e.g. RNA1:5:R1
					String targetR = edge.getTargetUID();
					Map<String, MolAtom> sourceAtomMap = sourceStruc
							.getRgroupMap();
					Map<String, MolAtom> targetAtomMap = targetStruc
							.getRgroupMap();

					Molecule tempMol = sourceStruc.getMolecule();
					StructureParser.merge(tempMol,
							(MolAtom) sourceAtomMap.get(sourceR),
							targetStruc.getMolecule(),
							(MolAtom) targetAtomMap.get(targetR));

					// new structure after connection
					RgroupStructure dirtyStruc = new RgroupStructure();
					dirtyStruc.setMolecule(tempMol);
					Map<String, MolAtom> rgroupMap = new HashMap<String, MolAtom>();
					for (String key : sourceAtomMap.keySet()) {
						MolAtom atom = sourceAtomMap.get(key);
						rgroupMap.put(key, atom);
					}
					for (String key : targetAtomMap.keySet()) {
						MolAtom atom = targetAtomMap.get(key);
						rgroupMap.put(key, atom);
					}
					rgroupMap.remove(sourceR);
					rgroupMap.remove(targetR);
					dirtyStruc.setRgroupMap(rgroupMap);

					// new ID list
					List<String> newIDList = new ArrayList<String>();
					if (sourceIDList.isEmpty() && targetIDList.isEmpty()) {
						newIDList.add(sourceId);
						newIDList.add(targetId);
					} else if (sourceIDList.isEmpty()
							&& !targetIDList.isEmpty()) {
						newIDList.addAll(targetIDList);
						if (!newIDList.contains(sourceId)) {
							newIDList.add(sourceId);
						}
						mergedNodeStrucMap.remove(targetIDList);

					} else if (!sourceIDList.isEmpty()
							&& targetIDList.isEmpty()) {
						newIDList.addAll(sourceIDList);
						if (!newIDList.contains(targetId)) {
							newIDList.add(targetId);
						}
						mergedNodeStrucMap.remove(sourceIDList);
					} else {
						newIDList.addAll(sourceIDList);
						for (String id : targetIDList) {
							if (!newIDList.contains(id)) {
								newIDList.add(id);
							}
						}
						mergedNodeStrucMap.remove(sourceIDList);
						mergedNodeStrucMap.remove(targetIDList);
					}
					mergedNodeStrucMap.put(newIDList, dirtyStruc);
					connEdgeList.remove(edge);
				}

				// only one entry finally
				RgroupStructure mergerdStru = mergedNodeStrucMap.values()
						.iterator().next();
				groupStructureMap.put(group, mergerdStru);
			}
		}

		// remove all remaining R groups not used in inter polymer connection
		List<Molecule> l = getMoleculeList(nodeList, groupStructureMap,
				monomerStoreToUse);
		return l;
	}

	/*
	 * private static List<Molecule> getMoleculeList(List<PolymerNode> nodeList,
	 * Map<Integer, RgroupStructure> groupStructureMap) throws
	 * NotationException, MonomerException, IOException, JDOMException,
	 * StructureException { MonomerFactory factory =
	 * MonomerFactory.getInstance(); MonomerStore
	 * store=factory.getMonomerStore(); return
	 * getMoleculeList(nodeList,groupStructureMap,store); }
	 */

	private static List<Molecule> getMoleculeList(List<PolymerNode> nodeList,
			Map<Integer, RgroupStructure> groupStructureMap,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		List<Molecule> l = new ArrayList<Molecule>();
		Collection<RgroupStructure> c = groupStructureMap.values();
		for (Iterator i = c.iterator(); i.hasNext();) {
			RgroupStructure struc = (RgroupStructure) i.next();
			Molecule mol = struc.getMolecule();
			Map<String, MolAtom> rgroupMap = struc.getRgroupMap();

			Set keyset = rgroupMap.keySet();
			for (Iterator it = keyset.iterator(); it.hasNext();) {
				String key = (String) it.next(); // NodeID:MonomerNumber:R#
				String[] keyItems = key
						.split(PolymerEdge.MONOMER_ATTACHEMENT_SEPARATOR);
				String nodeId = keyItems[0];
				String polymerType = PolymerNode.getPolymerType(nodeId);
				String nodeNotation = getPolymerNotation(nodeId, nodeList);
				List<String> monomerIDList = SimpleNotationParser
						.getMonomerIDList(nodeNotation, polymerType,
								monomerStore);

				int monomerNumber = Integer.parseInt(keyItems[1]);
				String monomerId = monomerIDList.get(monomerNumber - 1);
				Monomer monomer = SimpleNotationParser.getMonomer(monomerId,
						polymerType, monomerStore);
				List<Attachment> attachments = monomer.getAttachmentList();

				String rlabel = keyItems[2];
				for (int j = 0; j < attachments.size(); j++) {
					Attachment att = attachments.get(j);
					if (att.getCapGroupSMILES() == null) {
						MonomerParser.fillAttachmentInfo(att);
					}
					if (att.getLabel().equals(rlabel)) {
						Molecule attMol = StructureParser.getMolecule(att
								.getCapGroupSMILES());

						int rgroupId = Integer.parseInt(rlabel.substring(1));
						MolAtom attAtom = StructureParser.getRgroupAtom(attMol,
								rgroupId);
						StructureParser.merge(mol, rgroupMap.get(key), attMol,
								attAtom);
					}

				}
			}
			l.add(mol);
		}
		return l;
	}

	private static Entry<List<String>, RgroupStructure> getStructureEntry(
			String nodeId, Map<List<String>, RgroupStructure> mergedNodeStrucMap) {
		Entry<List<String>, RgroupStructure> result = null;
		for (Entry<List<String>, RgroupStructure> entry : mergedNodeStrucMap
				.entrySet()) {
			List<String> ids = entry.getKey();
			if (ids.contains(nodeId)) {
				result = entry;
				break;
			}
		}
		return result;
	}

	/**
	 * This method breaks PolymerNode lists into connected groups based on edge
	 * info List<String> could have one element for standalone polymer node,
	 * duplicated elements for self connection or multiple connection
	 * 
	 * @param allNodeList
	 * @param allEdgeList
	 * @return map of groupId and list of polymer node
	 */
	private static Map<Integer, List<String>> getGroupNodeListMap(
			List<PolymerNode> allNodeList, List<PolymerEdge> allEdgeList) {
		Map<Integer, List<String>> groupNodeListMap = new HashMap<Integer, List<String>>();
		int groupID = 0;

		// put connected nodes into groups
		for (int i = 0; i < allEdgeList.size(); i++) {
			PolymerEdge edge = allEdgeList.get(i);
			String sourceId = edge.getSourceNode();
			String targetId = edge.getTargetNode();
			List<Integer> overlapGroups = new ArrayList<Integer>();
			Set<Integer> groupSet = groupNodeListMap.keySet();
			for (Integer in : groupSet) {
				List<String> li = groupNodeListMap.get(in);
				if (li.contains(sourceId) || li.contains(targetId)) {
					overlapGroups.add(in);
				}
			}

			if (overlapGroups.isEmpty()) {
				groupID++;
				List<String> l = new ArrayList<String>();
				l.add(sourceId);
				l.add(targetId);
				groupNodeListMap.put(new Integer(groupID), l);
			} else {
				if (overlapGroups.size() > 1) {
					for (int j = 0; j < overlapGroups.size(); j++) {
						Integer i1 = overlapGroups.get(j);
						List<String> l1 = groupNodeListMap.get(i1);
						for (int k = j + 1; k < overlapGroups.size(); k++) {
							Integer i2 = overlapGroups.get(k);
							List<String> l2 = groupNodeListMap.get(i2);
							l1.addAll(l2);
							groupNodeListMap.remove(i2);
						}
					}
				}

				Integer i3 = overlapGroups.get(0);
				List<String> l3 = groupNodeListMap.get(i3);
				l3.add(sourceId);
				l3.add(targetId);
			}
		}

		// put standalone node into separate group
		for (int i = 0; i < allNodeList.size(); i++) {
			PolymerNode node = allNodeList.get(i);
			String nodeId = node.getId();
			boolean foundParent = false;
			Set<Integer> groupSet = groupNodeListMap.keySet();
			for (Integer in : groupSet) {
				List<String> li = groupNodeListMap.get(in);
				if (li.contains(nodeId)) {
					foundParent = true;
					break;
				}
			}
			if (!foundParent) {
				groupID++;
				List<String> l = new ArrayList<String>();
				l.add(nodeId);
				groupNodeListMap.put(new Integer(groupID), l);
			}
		}

		return groupNodeListMap;
	}

	/**
	 * This methods generates a List of PolymerEdge for the connected Nodes
	 * 
	 * @param connectedNodeList
	 * @param allEdgeList
	 * @return list of PolymerEdge that are connected
	 */
	private static List<PolymerEdge> getConnectedEdgeList(
			List<String> connectedNodeList, List<PolymerEdge> allEdgeList) {
		List<PolymerEdge> l = new ArrayList<PolymerEdge>();
		for (int i = 0; i < allEdgeList.size(); i++) {
			PolymerEdge edge = allEdgeList.get(i);
			String sourceId = edge.getSourceNode();
			String targetId = edge.getTargetNode();
			if (connectedNodeList.contains(sourceId)
					&& connectedNodeList.contains(targetId)) {
				l.add(edge);
			}
		}
		return l;
	}

	public static boolean validateComplexNotation(String extendedNotation)
			throws NotationException, MonomerException, IOException,
			StructureException, JDOMException {
		return validateComplexNotation(extendedNotation, null);
	}

	/**
	 * This methods validates the complex polymer extendendNotation
	 * 
	 * @param extendedNotation
	 * @return true or false
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean validateComplexNotation(String extendedNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		// validate
		validateNotationFormat(extendedNotation);
		ComplexPolymer cp = parse(extendedNotation, monomerStoreToUse);
		validateComplexPolymer(cp, monomerStoreToUse);

		return true;
	}

	/**
	 * This methods checks Nodes and Edges matches, and checks validity of each
	 * node label
	 * 
	 * @param complexPolymer
	 * @return true or false
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.helm.notation.StructureException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean validateComplexPolymer(ComplexPolymer complexPolymer)
			throws NotationException, IOException, MonomerException,
			StructureException, JDOMException {

		return validateComplexPolymer(complexPolymer, null);
	}

	public static boolean validateComplexPolymer(ComplexPolymer complexPolymer,
			MonomerStore monomerStore) throws NotationException, IOException,
			MonomerException, StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		List<PolymerNode> nodeList = complexPolymer.getPolymerNodeList();
		List<PolymerEdge> edgeList = complexPolymer.getPolymerEdgeList();
		List<PolymerEdge> bpList = complexPolymer.getBasePairList();
		Map<String, String> annotationMap = complexPolymer
				.getPolymerNodeAnnotationMap();

		// must
		if (null == nodeList || nodeList.size() == 0) {
			throw new NotationException(
					"Complex notation must have at least one polymer node");
		}
		Map<String, String> nodeMap = new HashMap<String, String>();
		for (int i = 0; i < nodeList.size(); i++) {
			PolymerNode node = nodeList.get(i);
			nodeMap.put(node.getId(), node.getLabel());
		}
		if (nodeList.size() > nodeMap.size()) {
			throw new NotationException("Polymer node IDs are not unique");
		}

		// this call also validates simple polymers
		Map<String, RgroupStructure> polymerNodeStructureMap = getPolymerNodeStructureMap(
				nodeList, monomerStoreToUse);

		if (null != edgeList) {
			for (PolymerEdge edge : edgeList) {
				if (edge.getEdgeType() == PolymerEdge.STANDARD_EDGE) {
					validateStandardEdge(edge, nodeMap, polymerNodeStructureMap);
				} else if (edge.getEdgeType() == PolymerEdge.GENERIC_EDGE) {
					validateGenericEdge(edge, nodeMap, polymerNodeStructureMap);
				} else if (edge.getEdgeType() == PolymerEdge.PAIR_EDGE) {
					validatePairEdge(edge, nodeMap);
				} else {
					throw new NotationException("Invalid edge type: "
							+ edge.getEdgeType());
				}
			}
		}

		if (null != bpList) {
			for (int i = 0; i < bpList.size(); i++) {
				PolymerEdge edge = bpList.get(i);
				validatePairEdge(edge, nodeMap);
			}
		}

		if (null != annotationMap && annotationMap.size() > 0) {
			Set ids = annotationMap.keySet();
			for (Iterator i = ids.iterator(); i.hasNext();) {
				String id = (String) i.next();
				if (!nodeMap.containsKey(id)) {
					throw new NotationException(
							"Polymer annotation contains unknown polymer node ID");
				}
			}
		}

		return true;
	}

	private static void validateStandardEdge(PolymerEdge edge,
			Map<String, String> nodeMap,
			Map<String, RgroupStructure> polymerNodeStructureMap)
			throws NotationException {
		// source node check
		String node = edge.getSourceNode();
		String connection = edge.getSourceConnection();
		validateStandardAttachment(node, connection, nodeMap,
				polymerNodeStructureMap);

		// target node check
		node = edge.getTargetNode();
		connection = edge.getTargetConnection();
		validateStandardAttachment(node, connection, nodeMap,
				polymerNodeStructureMap);
	}

	private static void validatePairEdge(PolymerEdge edge,
			Map<String, String> nodeMap) throws NotationException {
		if (edge.getEdgeType() != PolymerEdge.PAIR_EDGE) {
			throw new NotationException("Invalid base pair edge: "
					+ edge.getEdgeNotation() + " => " + edge.toString());
		}

		if (!(nodeMap.containsKey(edge.getSourceNode()))) {
			throw new NotationException(
					"Polymer edge contains unknown polymer node ID");
		}
		if (!(nodeMap.containsKey(edge.getTargetNode()))) {
			throw new NotationException(
					"Polymer edge contains unknown polymer node ID");
		}
	}

	private static void validateGenericEdge(PolymerEdge edge,
			Map<String, String> nodeMap,
			Map<String, RgroupStructure> polymerNodeStructureMap)
			throws NotationException {
		// source node check
		int attType = edge.getSourceAttachmentType();
		if (attType == PolymerEdge.STANDARD_EDGE_ATTACHMENT) {
			String node = edge.getSourceNode();
			String connection = edge.getSourceConnection();
			validateStandardAttachment(node, connection, nodeMap,
					polymerNodeStructureMap);
		} else if (attType == PolymerEdge.GENERIC_EDGE_ATTACHMENT) {
			String[] nodes = edge.getSourceNodes();
			validateGenericAttachment(nodes, nodeMap);
		}

		// target node check
		attType = edge.getTargetAttachmentType();
		if (attType == PolymerEdge.STANDARD_EDGE_ATTACHMENT) {
			String node = edge.getTargetNode();
			String connection = edge.getTargetConnection();
			validateStandardAttachment(node, connection, nodeMap,
					polymerNodeStructureMap);
		} else if (attType == PolymerEdge.GENERIC_EDGE_ATTACHMENT) {
			String[] nodes = edge.getTargetNodes();
			validateGenericAttachment(nodes, nodeMap);
		}
	}

	private static void validateStandardAttachment(String node,
			String connection, Map<String, String> nodeMap,
			Map<String, RgroupStructure> polymerNodeStructureMap)
			throws NotationException {
		RgroupStructure rstructure = polymerNodeStructureMap.get(node);

		if (!(nodeMap.containsKey(node))) {
			throw new NotationException(
					"Polymer edge contains unknown polymer node ID");
		}

		if (null == rstructure) {
			throw new NotationException(
					"Polymer edge contains polymer node without connection point");
		}

		Map<String, MolAtom> rMap = rstructure.getRgroupMap();
		if (null == rMap || !rMap.containsKey(connection)) {
			throw new NotationException(
					"Polymer edge contains polymer node without connection point");
		}
	}

	private static void validateGenericAttachment(String[] nodes,
			Map<String, String> nodeMap) throws NotationException {
		for (String node : nodes) {
			if (!(nodeMap.containsKey(node))) {
				throw new NotationException(
						"Polymer edge contains unknown polymer node ID");
			}
		}
	}

	/**
	 * This method convert the complex polymer extendendNotation into an
	 * ComplexPolymer object Will throw NotationException if graph
	 * extendendNotation format is invalid,
	 * 
	 * @param extendendNotation
	 * @return ComplexPolymer
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 */
	public static ComplexPolymer parse(String extendendNotation)
			throws NotationException, MonomerException, JDOMException,
			IOException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return parse(extendendNotation, factory.getMonomerStore());
	}

	public static ComplexPolymer parse(String extendendNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, JDOMException, IOException {
		ComplexPolymer cp = new ComplexPolymer();

		String allNodeString = getAllNodeString(extendendNotation);
		List<PolymerNode> nodes = getPolymerNodeList(allNodeString,
				monomerStore);
		cp.setPolymerNodeList(nodes);

		String allEdgeString = getAllEdgeString(extendendNotation);
		List<PolymerEdge> edges = getPolymerEdgeList(allEdgeString);
		cp.setPolymerEdgeList(edges);

		String allBasePairString = getAllBasePairString(extendendNotation);
		List<PolymerEdge> basePairEdges = getPolymerEdgeList(allBasePairString);
		cp.setBasePairList(basePairEdges);

		String allNodeLabelString = getAllNodeLabelString(extendendNotation);
		Map<String, String> annotationMap = getPolymerNodeIDAnnotationMap(allNodeLabelString);
		cp.setPolymerNodeAnnotationMap(annotationMap);

		return cp;
	}

	public static String getAllNodeString(String extendedNotation)
			throws NotationException {
		return getComponentString(extendedNotation, 0, false);
	}

	public static String getAllEdgeString(String extendedNotation)
			throws NotationException {
		return getComponentString(extendedNotation, 1, false);
	}

	public static String getAllBasePairString(String extendedNotation)
			throws NotationException {
		return getComponentString(extendedNotation, 2, false);
	}

	public static String getAllNodeLabelString(String extendedNotation)
			throws NotationException {
		return getComponentString(extendedNotation, 3, false);
	}

	public static String getOtherString(String extendedNotation)
			throws NotationException {
		return getComponentString(extendedNotation, 4, true);
	}

	/**
	 * This method returns the string for the specified component
	 * 
	 * @param extendedNotation
	 * @param componentIndex
	 *            0 based component index, i.e, first component has an index of
	 *            0
	 * @return string of component
	 */
	private static String getComponentString(String extendedNotation,
			int componentIndex, boolean isLast) throws NotationException {
		validateNotationFormat(extendedNotation);
		int count = 0;
		int startPos = 0;
		int endPos = 0;

		int firstPos = extendedNotation.indexOf(TOP_LEVEL_STOPPER, 0) + 1;

		if (componentIndex == 0) {
			endPos = firstPos;
		} else {
			endPos = firstPos;
			while (count < componentIndex) {
				startPos = endPos + 1;
				endPos = extendedNotation
						.indexOf(TOP_LEVEL_DELIMITER, startPos);
				count++;
			}
		}

		if (isLast) {
			return extendedNotation.substring(startPos);
		} else {
			return extendedNotation.substring(startPos, endPos);
		}
	}

	public static boolean validateNotationFormat(String extendedNotation)
			throws NotationException {
		int count = 0;
		int startPos = 0;
		int pos = 0;

		pos = extendedNotation.indexOf(TOP_LEVEL_STOPPER, 0);
		if (pos < 0) {
			throw new NotationException(
					"Invalid complex notation format,missing first component stopper }$");
		} else {
			startPos = pos + TOP_LEVEL_STOPPER.length();
		}

		while (count < 3) {
			pos = extendedNotation.indexOf(TOP_LEVEL_DELIMITER, startPos);
			if (pos < 0) {
				throw new NotationException(
						"Invalid complex notation format, must have four positional delimieters $$$$");
			} else {
				startPos = pos + 1;
				count++;
			}
		}

		return true;
	}

	/**
	 * This methods converts the allNodeString to a List of PolymerNode Note on
	 * 2/25/2011 by Tianhong: CHEM nodes could contain adhoc monomers, which
	 * contains SMILES string in the format of ID?SMILES SMILES could have
	 * LIST-LEVEL_DELIMITER which is pipe | Need to find them and do a
	 * preprocess, add adhoc monomer into monomer factory
	 * 
	 * @param allNodeString
	 *            - allNodeString or complex notation string
	 * @return List<PolymerNode>
	 * @throws org.helm.notation.NotationException
	 */
	public static List<PolymerNode> getPolymerNodeList(String allNodeString)
			throws NotationException {

		return getPolymerNodeList(allNodeString, null);
	}

	public static List<PolymerNode> getPolymerNodeList(String allNodeString,
			MonomerStore monomerStore) throws NotationException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		List<PolymerNode> list = new ArrayList<PolymerNode>();
		if (null != allNodeString && allNodeString.length() > 0) {
			String nodes = allNodeString;
			if (allNodeString.indexOf(TOP_LEVEL_STOPPER) > 0) {
				nodes = getAllNodeString(allNodeString);
			}

			int nodeStartPos = 0;
			int delimiterStartPos = 0;
			int endPos = 0;
			String nodeString = null;

			delimiterStartPos = nodes.indexOf(NODE_LABEL_START_SYMBOL,
					nodeStartPos);
			endPos = nodes.indexOf(NODE_LABEL_END_SYMBOL, nodeStartPos) + 1;

			while (delimiterStartPos > 0 && endPos > 0) {
				nodeString = nodes.substring(nodeStartPos, endPos);
				PolymerNode node = getPolymerNode(nodeString, monomerStoreToUse);
				list.add(node);

				nodeStartPos = endPos + 1;
				delimiterStartPos = nodes.indexOf(NODE_LABEL_START_SYMBOL,
						nodeStartPos);
				endPos = nodes.indexOf(NODE_LABEL_END_SYMBOL, nodeStartPos) + 1;
			}
		}
		return list;
	}

	/**
	 * This method converts a nodeString into a PolymerNode object, throws
	 * validation error if invalid
	 * 
	 * @param nodeString
	 *            , the string for each polymer node in the complex
	 *            extendendNotation
	 * @return and object of PolymerNode
	 * @throws org.helm.notation.NotationException
	 */

	private static PolymerNode getPolymerNode(String nodeString)
			throws NotationException {

		return getPolymerNode(nodeString, null);
	}

	private static PolymerNode getPolymerNode(String nodeString,
			MonomerStore monomerStore) throws NotationException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		int startSymbolPos = nodeString.indexOf(NODE_LABEL_START_SYMBOL);
		int endSymbolPos = nodeString.indexOf(NODE_LABEL_END_SYMBOL);

		if (startSymbolPos <= 0 || endSymbolPos <= 0) {
			throw new NotationException(INVALID_POLYMER_NODE + ": "
					+ nodeString);
		}
		if (endSymbolPos != nodeString.length() - 1) {
			throw new NotationException(INVALID_POLYMER_NODE + ": "
					+ nodeString);
		}
		String id = nodeString.substring(0, startSymbolPos);
		validateNodeID(id);
		String label = nodeString.substring(startSymbolPos + 1, endSymbolPos);

		// add adhoc chem monomer into monomer database if adhoc
		if (id.startsWith(Monomer.CHEMICAL_POLYMER_TYPE)) {
//			label = SimpleNotationParser.preprocessChemNode(label,
//					monomerStoreToUse);
			
			label = SimpleNotationParser.processNode(label,Monomer.CHEMICAL_POLYMER_TYPE,monomerStoreToUse);
		}

		PolymerNode node = new PolymerNode();
		node.setId(id);
		node.setLabel(label);
		return node;
	}

	/**
	 * This method returns the polymer extendendNotation for each polymer node
	 * 
	 * @param nodeId
	 * @param nodeList
	 * @return polymer extendendNotation for the given node
	 */
	private static String getPolymerNotation(String nodeId,
			List<PolymerNode> nodeList) {
		String notation = null;

		for (int i = 0; i < nodeList.size(); i++) {
			PolymerNode node = nodeList.get(i);
			if (node.getId().equals(nodeId)) {
				notation = node.getLabel();
				break;
			}
		}

		return notation;
	}

	/**
	 * nodeID must be in the format of Letters followed by Numbers, throws
	 * NotationException if invalid
	 * 
	 * @param nodeID
	 * @throws org.helm.notation.NotationException
	 */
	private static boolean validateNodeID(String nodeID)
			throws NotationException {
		char[] chars = nodeID.toCharArray();
		if (!(String.valueOf(chars[0]).matches("[A-Za-z]"))) {
			throw new NotationException(INVALID_NODE_ID + ": " + nodeID);
		}
		if (!(String.valueOf(chars[chars.length - 1]).matches("[0-9]"))) {
			throw new NotationException(INVALID_NODE_ID + ": " + nodeID);
		}
		boolean foundNum = false;

		for (int i = 0; i < chars.length; i++) {
			if (String.valueOf(chars[i]).matches("[0-9]")) {
				foundNum = true;
			}
			if (foundNum) {
				if ((String.valueOf(chars[i]).matches("[A-Za-z]"))) {
					throw new NotationException(INVALID_NODE_ID + ": " + nodeID);
				}
			}
		}

		return true;
	}

	/**
	 * This methods converts the allEdgeString to a List of PolymerEdge, could
	 * contain empty or invalid edges
	 * 
	 * @param allEdgeString
	 * @return List<PolymerEdge>
	 * @throws org.helm.notation.NotationException
	 */
	private static List<PolymerEdge> getPolymerEdgeList(String allEdgeString)
			throws NotationException, MonomerException {

		List<PolymerEdge> list = new ArrayList<PolymerEdge>();
		if (null != allEdgeString && allEdgeString.length() > 0) {
			String[] allEdges = allEdgeString.split(LIST_LEVEL_DELIMITER_REGEX);
			for (int i = 0; i < allEdges.length; i++) {
				PolymerEdge edge = EdgeParser.parse(allEdges[i]);
				list.add(edge);
			}
		}
		return list;
	}

	public static List<RNAPolymerNode> getRNAPolymerNodeList(
			String complexNotation) throws NotationException, MonomerException,
			IOException, JDOMException, StructureException {
		List<PolymerNode> list = getPolymerNodeList(getAllNodeString(complexNotation));
		String allNodeAnnotationString = getAllNodeLabelString(complexNotation);
		Map<String, String> map = getPolymerNodeIDAnnotationMap(allNodeAnnotationString);
		List<RNAPolymerNode> l = new ArrayList<RNAPolymerNode>();
		for (int i = 0; i < list.size(); i++) {
			PolymerNode node = list.get(i);
			if (node.getType().equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
				RNAPolymerNode rnaNode = new RNAPolymerNode(node);
				List<Nucleotide> nucList = SimpleNotationParser
						.getNucleotideList(node.getLabel());
				String seq = SimpleNotationParser
						.getNucleotideSequence(nucList);
				rnaNode.setSequence(seq);
				String modSeq = SimpleNotationParser
						.getModifiedNucleotideSequence(nucList);
				rnaNode.setModifiedSequence(modSeq);
				if (map.containsKey(rnaNode.getId())) {
					rnaNode.setAnotation(map.get(rnaNode.getId()));
				}
				l.add(rnaNode);
			}
		}

		return l;
	}

	/**
	 * generate formated siRNA sequence with default padding char " " and
	 * base-pair char "|"
	 * 
	 * @param complexNotaton
	 * @return string array of formated nucleotide sequence
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.StructureException
	 */
	public static String[] getFormatedSirnaSequences(String complexNotaton)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		return getFormatedSirnaSequences(complexNotaton, DEFAULT_PADDING_CHAR,
				DEFAULT_BASE_PAIR_CHAR);
	}

	/**
	 * generates formated siRNA sequence with provided padding and base-pair
	 * character, such as AGCUUUGGTT |||||||| TTUCGAAACC
	 * 
	 * @param complexNotation
	 * @param paddingChar
	 * @param basePairChar
	 * @return string array of formated nucleotide sequence
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.StructureException
	 */
	public static String[] getFormatedSirnaSequences(String complexNotation,
			String paddingChar, String basePairChar) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		if (null == paddingChar || paddingChar.length() != 1) {
			throw new NotationException(
					"Padding string must be single character");
		}

		if (null == basePairChar || basePairChar.length() != 1) {
			throw new NotationException(
					"Base pair string must be single character");
		}

		List<RNAPolymerNode> rnaList = getRNAPolymerNodeList(complexNotation);
		int count = rnaList.size();
		if (count == 0) {
			return new String[0];
		} else if (count == 1) {
			return new String[] { rnaList.get(0).getSequence() };
		} else if (count == 2) {
			String rna1Seq = null;
			String rna2Seq = null;
			String rna1Annotation = null;
			String rna2Annotation = null;

			for (RNAPolymerNode node : rnaList) {
				if (node.getId().equals("RNA1")) {
					rna1Seq = node.getSequence();
					rna1Annotation = node.getAnotation();
				} else if (node.getId().equals("RNA2")) {
					rna2Seq = node.getSequence();
					rna2Annotation = node.getAnotation();
				}
			}
			String reverseRna2Seq = NucleotideSequenceParser
					.getReverseSequence(rna2Seq);

			String bpString = getAllBasePairString(complexNotation);
			if (null == bpString || bpString.length() == 0) {
				return new String[] { rna1Seq, rna2Seq };
			} else {
				Map<Integer, Integer> monomerPositionMap = getSirnaMonomerPositionMap(bpString);
				Map<Integer, Integer> seqPositionMap = new HashMap<Integer, Integer>();
				Set<Integer> monomerSet = monomerPositionMap.keySet();
				for (Integer key : monomerSet) {
					Integer value = monomerPositionMap.get(key);
					Integer seqKey = new Integer(key.intValue() / 3 + 1);
					Integer seqValue = new Integer(value.intValue() / 3 + 1);
					seqPositionMap.put(seqKey, seqValue);
				}

				Set<Integer> seqSet = seqPositionMap.keySet();
				List<Integer> seqList = new ArrayList<Integer>();
				for (Integer key : seqSet) {
					seqList.add(key);
				}
				Collections.sort(seqList);

				int rna1First = seqList.get(0).intValue();
				int rna2Last = seqPositionMap.get(seqList.get(0)).intValue();
				int rna1Last = seqList.get(seqList.size() - 1).intValue();
				int rna2First = seqPositionMap.get(
						seqList.get(seqList.size() - 1)).intValue();

				if ((rna1Last - rna1First) != (rna2Last - rna2First)) {
					throw new NotationException(
							"siRNA matching lengths are different");
				}

				int rna1LeftOverhang = rna1First - 1;
				int rna1RightOverhang = rna1Seq.length() - rna1Last;
				int rna2LeftOverhang = rna2Seq.length() - rna2Last;
				int rna2RightOverhang = rna2First - 1;
				StringBuffer[] sbs = new StringBuffer[3];
				for (int i = 0; i < sbs.length; i++) {
					sbs[i] = new StringBuffer();
				}

				if (rna1LeftOverhang >= rna2LeftOverhang) {
					sbs[0].append(rna1Seq);

					for (int i = 0; i < rna1LeftOverhang; i++) {
						sbs[1].append(paddingChar);
					}
					for (int i = rna1First; i < (rna1Last + 1); i++) {
						Integer in = new Integer(i);
						if (seqPositionMap.containsKey(in)) {
							sbs[1].append(basePairChar);
						} else {
							sbs[1].append(paddingChar);
						}
					}

					for (int i = 0; i < rna1LeftOverhang - rna2LeftOverhang; i++) {
						sbs[2].append(paddingChar);
					}
					sbs[2].append(reverseRna2Seq);
				} else {
					for (int i = 0; i < rna2LeftOverhang - rna1LeftOverhang; i++) {
						sbs[0].append(paddingChar);
					}
					sbs[0].append(rna1Seq);

					for (int i = 0; i < rna2LeftOverhang; i++) {
						sbs[1].append(paddingChar);
					}
					for (int i = rna1First; i < (rna1Last + 1); i++) {
						Integer in = new Integer(i);
						if (seqPositionMap.containsKey(in)) {
							sbs[1].append(basePairChar);
						} else {
							sbs[1].append(paddingChar);
						}
					}

					sbs[2].append(reverseRna2Seq);
				}

				if (rna1RightOverhang >= rna2RightOverhang) {
					for (int i = 0; i < rna1RightOverhang; i++) {
						sbs[1].append(paddingChar);
					}

					for (int i = 0; i < rna1RightOverhang - rna2RightOverhang; i++) {
						sbs[2].append(paddingChar);
					}
				} else {
					for (int i = 0; i < rna2RightOverhang - rna1RightOverhang; i++) {
						sbs[0].append(paddingChar);
					}

					for (int i = 0; i < rna2RightOverhang - rna1RightOverhang; i++) {
						sbs[1].append(paddingChar);
					}
				}

				if ((rna1Annotation != null && rna1Annotation
						.equalsIgnoreCase("AS"))
						|| (rna2Annotation != null && rna2Annotation
								.equalsIgnoreCase("SS"))) {
					return new String[] { reverseString(sbs[2].toString()),
							reverseString(sbs[1].toString()),
							reverseString(sbs[0].toString()) };
				} else {
					return new String[] { sbs[0].toString(), sbs[1].toString(),
							sbs[2].toString() };
				}
			}

		} else {
			throw new NotationException(
					"Structure contains more than two RNA sequences");
		}

	}

	private static Map<Integer, Integer> getSirnaMonomerPositionMap(
			String basePairString) throws NotationException {
		String[] pairs = basePairString.split(LIST_LEVEL_DELIMITER_REGEX);
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		for (String pair : pairs) {
			String[] edgeComponents = pair
					.split(PolymerEdge.EDGE_COMPONENT_DELIMITER);
			String pairInfo = edgeComponents[2];
			String[] edgeProperties = pairInfo
					.split(PolymerEdge.CONNECTION_SEPARATOR);

			String str1 = edgeProperties[0]
					.split(PolymerEdge.MONOMER_ATTACHEMENT_SEPARATOR)[0];
			Integer int1 = new Integer(str1);

			String str2 = edgeProperties[1]
					.split(PolymerEdge.MONOMER_ATTACHEMENT_SEPARATOR)[0];
			Integer int2 = new Integer(str2);

			if (pair.startsWith("RNA1")) {
				map.put(int1, int2);
			} else if (pair.startsWith("RNA2")) {
				map.put(int2, int1);
			} else {
				throw new NotationException(
						"Structure contains more than two RNA sequences");
			}
		}
		return map;
	}

	private static String reverseString(String source) {
		int i;
		int len = source.length();
		StringBuffer dest = new StringBuffer();

		for (i = (len - 1); i >= 0; i--) {
			dest.append(source.charAt(i));
		}
		return dest.toString();
	}

	private static Map<String, String> getPolymerNodeIDAnnotationMap(
			String allNodeAnnotationString) {
		Map<String, String> map = new HashMap<String, String>();
		if (null != allNodeAnnotationString
				&& allNodeAnnotationString.trim().length() > 0) {
			String[] items = allNodeAnnotationString
					.split(LIST_LEVEL_DELIMITER_REGEX);
			for (int i = 0; i < items.length; i++) {
				String item = items[i];
				int start = item.indexOf(NODE_LABEL_START_SYMBOL);
				int end = item.indexOf(NODE_LABEL_END_SYMBOL);
				String nodeId = item.substring(0, start);
				String annotation = item.substring(start + 1, end);
				map.put(nodeId, annotation);
			}
		}
		return map;
	}

	/**
	 * Generate canonical notation without validation
	 * 
	 * @param complexNotation
	 * @return canonical notation
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws ClassNotFoundException
	 * @throws StructureException
	 * @throws JDOMException
	 */
	public static String getCanonicalNotation(String complexNotation)
			throws NotationException, MonomerException, IOException,
			ClassNotFoundException, StructureException, JDOMException {

		return getCanonicalNotation(complexNotation, null);
	}

	public static String getCanonicalNotation(String complexNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, ClassNotFoundException,
			StructureException, JDOMException {

		return getCanonicalNotation(complexNotation, false, monomerStore);
	}

	/**
	 * Generate canonical notation for the the input complex notation
	 * 
	 * @param complexNotation
	 * @param includeValidation
	 * @return canonical notation
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws ClassNotFoundException
	 * @throws StructureException
	 * @throws JDOMException
	 */

	public static String getCanonicalNotation(String complexNotation,
			boolean includeValidation) throws NotationException,
			MonomerException, IOException, ClassNotFoundException,
			StructureException, JDOMException {

		return getCanonicalNotation(complexNotation, includeValidation, null);
	}

	public static String getCanonicalNotation(String complexNotation,
			boolean includeValidation, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			ClassNotFoundException, StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		validateNotationFormat(complexNotation);
		ComplexPolymer cp = parse(complexNotation, monomerStoreToUse);

		if (includeValidation) {
			validateComplexPolymer(cp, monomerStoreToUse);
		}

		List<PolymerNode> nodeList = cp.getPolymerNodeList();
		List<PolymerEdge> edgeList = cp.getPolymerEdgeList();
		// TODO deal with adhoc monomers in all polymer types
		// deal with ad hoc CHEM monomer here, use smiles instead of temp ID
		for (PolymerNode node : nodeList) {
			if (node.getType().equals(Monomer.CHEMICAL_POLYMER_TYPE)
					&& node.getLabel().startsWith(
							SimpleNotationParser.AD_HOC_MONOMER_ID_PREFIX)) {

				Monomer m = monomerStoreToUse.getMonomer(
						Monomer.CHEMICAL_POLYMER_TYPE, node.getLabel());

				String smiles = m.getCanSMILES();
				String uniSmi = StructureParser.getUniqueExtendedSMILES(smiles);

				List<String> rgroups = StructureParser
						.getRGroupsFromExtendedSmiles(uniSmi);
				Map<String, String> oldNewRgroupMap = new HashMap<String, String>();
				for (int i = 0; i < rgroups.size(); i++) {
					String oldId = rgroups.get(i);
					String newId = "Q" + (i + 1);
					oldNewRgroupMap.put(oldId, newId);
					uniSmi = uniSmi.replaceAll(oldId, newId);
				}
				uniSmi = uniSmi.replaceAll("Q", "R");
				node.setLabel(uniSmi);

				String nodeId = node.getId();

				// deal with edges
				for (PolymerEdge edge : edgeList) {
					// find the edges of interest, could be more than one
					if (edge.getSourceNode().equals(nodeId)) {
						String oldR = edge.getSourceR();
						String newR = oldNewRgroupMap.get(oldR);
						if (null != newR) {
							newR = newR.replaceAll("Q", "R");
							edge.setSourceR(newR);
						}
					}

					if (edge.getTargetNode().equals(nodeId)) {
						String oldR = edge.getTargetR();
						String newR = oldNewRgroupMap.get(oldR);
						if (null != newR) {
							newR = newR.replaceAll("Q", "R");
							edge.setTargetR(newR);
						}
					}
				}
			}
		}

		// all backbone cycles
		List<String> backboneCycleNodes = new ArrayList<String>();
		List<String> branchCycleNodes = new ArrayList<String>();
		for (PolymerEdge edge : edgeList) {
			if (edge.isBackboneSelfCycle()) {
				backboneCycleNodes.add(edge.getSourceNode());
			} else if (edge.isSelfCycle()) {
				branchCycleNodes.add(edge.getSourceNode());
			}
		}

		// deal with backbone cycles self cycles: both node and edge needs to be
		// modified if backbone cycles are involved
		for (String nodeId : backboneCycleNodes) {
			// deal with node
			int offset = 0;
			int monomerCount = 0;
			for (PolymerNode node : nodeList) {
				// find the node of interest, should be only one
				if (node.getId().equals(nodeId)) {
					Map.Entry<Integer, String> canEntry = SimpleNotationParser
							.getSimpleCanonicalNotationMapEntry(
									node.getLabel(), node.getType());
					node.setLabel(canEntry.getValue());
					offset = canEntry.getKey().intValue();
					monomerCount = SimpleNotationParser.getMonomerCount(
							node.getLabel(), node.getType());
					break;
				}
			}

			// deal with edges
			for (PolymerEdge edge : edgeList) {
				// find the edges of interest, could be more than one
				if (edge.getSourceNode().equals(nodeId)
						|| edge.getTargetNode().equals(nodeId)) {
					if (edge.isBackboneSelfCycle()) {
						// canonicalize connection
						List<String> connections = new ArrayList<String>();
						connections.add(edge.getConnection());
						connections.add(edge.getReverseConnection());
						Collections.sort(connections);
						edge.setConnection(connections.get(0));
					} else {
						// need to modify monomer positions in connection string
						if (edge.getSourceNode().equals(nodeId)) {
							int oldSrcMonomerNum = edge
									.getSourceMonomerNumber();
							int newSrcMonomerNumber;
							if (oldSrcMonomerNum > offset) {
								newSrcMonomerNumber = oldSrcMonomerNum - offset;
							} else {
								newSrcMonomerNumber = monomerCount
										+ oldSrcMonomerNum - offset;
							}
							edge.setSourceMonomerNumber(newSrcMonomerNumber);
						}

						if (edge.getTargetNode().equals(nodeId)) {
							int oldTarMonomerNum = edge
									.getTargetMonomerNumber();
							int newTarMonomerNumber;
							if (oldTarMonomerNum > offset) {
								newTarMonomerNumber = oldTarMonomerNum - offset;
							} else {
								newTarMonomerNumber = monomerCount
										+ oldTarMonomerNum - offset;
							}
							edge.setTargetMonomerNumber(newTarMonomerNumber);
						}
					}
				}
			}
		}

		// deal with branch self cycles: self-connection edge should be
		// canonicalized
		for (String nodeId : branchCycleNodes) {
			// deal with edges
			for (PolymerEdge edge : edgeList) {
				// find the edges of interest, could be more than one
				if (edge.isSelfCycle() && !edge.isBackboneSelfCycle()
						&& edge.getSourceNode().equals(nodeId)) {
					// canonicalize connection
					List<String> connections = new ArrayList<String>();
					connections.add(edge.getConnection());
					connections.add(edge.getReverseConnection());
					Collections.sort(connections);
					edge.setConnection(connections.get(0));
				}
			}
		}

		Map<String, String> idLabelMap = new HashMap<String, String>();
		Map<String, List<String>> labelIdMap = new TreeMap<String, List<String>>();

		for (PolymerNode node : nodeList) {
			idLabelMap.put(node.getId(), node.getLabel());
			if (labelIdMap.containsKey(node.getLabel())) {
				List<String> l = labelIdMap.get(node.getLabel());
				l.add(node.getId());
			} else {
				List<String> l = new ArrayList<String>();
				l.add(node.getId());
				labelIdMap.put(node.getLabel(), l);
			}
		}

		Set<String> sortedLabelSet = labelIdMap.keySet();
		List<List<String[]>> lol = new ArrayList<List<String[]>>();
		for (String key : sortedLabelSet) {
			List<String> value = labelIdMap.get(key);
			List<String[]> al = new ArrayList<String[]>();
			al.add(value.toArray(new String[0]));
			// PermutationAndExpansion.permutate(al, value);
			PermutationAndExpansion.expand(lol, al);
		}

		List<List<String>> nodeIdPermutations = PermutationAndExpansion
				.linearize(lol);
		List<String> notationList = new ArrayList<String>();
		for (List<String> sortedIdList : nodeIdPermutations) {
			String notation = generateNotationBasedNodeOrder(sortedIdList,
					nodeList, edgeList);
			notationList.add(notation);
		}
		Collections.sort(notationList);
		return notationList.get(0);
	}

	private static String generateNotationBasedNodeOrder(
			List<String> sortedNodeIdList, List<PolymerNode> nodeList,
			List<PolymerEdge> edgeList) throws NotationException {
		Map<String, String> idLabelMap = new HashMap<String, String>();
		for (PolymerNode node : nodeList) {
			idLabelMap.put(node.getId(), node.getLabel());
		}

		Map<String, Integer> typeCountMap = new HashMap<String, Integer>();
		Map<String, String> oldNewIdMap = new HashMap<String, String>();

		StringBuilder sb = new StringBuilder();
		for (String oldId : sortedNodeIdList) {
			String label = idLabelMap.get(oldId);
			String type = oldId.split("\\d")[0];
			if (typeCountMap.containsKey(type)) {
				int count = typeCountMap.get(type).intValue();
				typeCountMap.put(type, new Integer(count + 1));
			} else {
				typeCountMap.put(type, new Integer(1));
			}
			int counter = typeCountMap.get(type).intValue();
			if (sb.length() > 0) {
				sb.append("|");
			}
			String newId = type + counter;
			oldNewIdMap.put(oldId, newId);
			sb.append(newId + "{" + label + "}");
		}
		sb.append("$");

		if (null != edgeList) {
			// canonicalize each polymer edge based on polymer node sort
			for (PolymerEdge edge : edgeList) {
				EdgeParser.canonicalize(edge, sortedNodeIdList);
			}

			List<PolymerEdge> sortedEdgeList = EdgeParser.sort(edgeList,
					sortedNodeIdList);

			// swap new node id with old node id
			for (PolymerEdge edge : sortedEdgeList) {
				String oldSource = edge.getSourceNode();
				String newSource = replaceNodeId(oldSource, oldNewIdMap);
				edge.setSourceNode(newSource);

				String oldTarget = edge.getTargetNode();
				String newTarget = replaceNodeId(oldTarget, oldNewIdMap);
				edge.setTargetNode(newTarget);
			}

			StringBuilder edgeSb = new StringBuilder();
			for (PolymerEdge edge : sortedEdgeList) {
				if (edgeSb.length() > 0) {
					edgeSb.append("|");
				}
				edgeSb.append(edge.toString());
			}

			sb.append(edgeSb.toString());
		}
		sb.append("$$$");
		return sb.toString();
	}

	private static String replaceNodeId(String oldNodeId,
			Map<String, String> oldNewIdMap) {
		String[] oldNodes = oldNodeId.split(
				PolymerEdge.NODE_CONCATENATOR_REGEX, -1);
		String[] newNodes = new String[oldNodes.length];
		for (int i = 0; i < oldNodes.length; i++) {
			newNodes[i] = oldNewIdMap.get(oldNodes[i]);
		}
		StringBuilder sb = new StringBuilder();
		for (String newNode : newNodes) {
			if (sb.length() > 0) {
				sb.append(PolymerEdge.NODE_CONCATENATOR);
			}
			sb.append(newNode);
		}
		return sb.toString();
	}

	/**
	 * The method generates the complex the notation for two complex notations
	 * combined. There will be no chemcial conection and base pairing between
	 * the components represented by each complex notation. *
	 * 
	 * @param complexNotation1
	 *            -- complex notation for the first component to be combined
	 * @param complexNotation2
	 *            -- complex notation for the second component to be combined
	 * @return the complex notation for the combined structure
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 */
	public static String getCombinedComlexNotation(String complexNotation1,
			String complexNotation2) throws NotationException, MonomerException {
		if (null == complexNotation1 || complexNotation1.trim().length() == 0) {
			return complexNotation2;
		}
		if (null == complexNotation2 || complexNotation2.trim().length() == 0) {
			return complexNotation1;
		}
		StringBuffer sb = new StringBuffer();
		int rnaCount = 0;
		int peptideCount = 0;
		int chemCount = 0;
		String nodeString1 = getAllNodeString(complexNotation1);
		String edgeString1 = getAllEdgeString(complexNotation1);
		String basePairString1 = getAllBasePairString(complexNotation1);
		String nodeLabelString1 = getAllNodeLabelString(complexNotation1);
		List<PolymerNode> nodeList1 = getPolymerNodeList(nodeString1);
		for (PolymerNode node : nodeList1) {
			String type = node.getType();
			if (type.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
				rnaCount++;
			} else if (type.equals(Monomer.PEPTIDE_POLYMER_TYPE)) {
				peptideCount++;
			} else if (type.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
				chemCount++;
			} else {
				throw new NotationException("Unknown polymer type found ["
						+ type + "]");
			}
		}

		Map<String, String> nodeIdMap = new HashMap<String, String>();
		String nodeString2 = getAllNodeString(complexNotation2);
		String edgeString2 = getAllEdgeString(complexNotation2);
		String basePairString2 = getAllBasePairString(complexNotation2);
		String nodeLabelString2 = getAllNodeLabelString(complexNotation2);
		List<PolymerNode> nodeList2 = getPolymerNodeList(nodeString2);
		for (PolymerNode node : nodeList2) {
			String type = node.getType();
			String oldId = node.getId();
			String newId = null;
			if (type.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
				rnaCount++;
				newId = type + rnaCount;
			} else if (type.equals(Monomer.PEPTIDE_POLYMER_TYPE)) {
				peptideCount++;
				newId = type + peptideCount;
			} else if (type.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
				chemCount++;
				newId = type + chemCount;
			} else {
				throw new NotationException("Unknown polymer type found ["
						+ type + "]");
			}
			nodeIdMap.put(oldId, newId);
		}

		Set keyset = nodeIdMap.keySet();
		for (Iterator i = keyset.iterator(); i.hasNext();) {
			String oldId = (String) i.next();
			String newId = nodeIdMap.get(oldId);
			nodeString2 = nodeString2.replaceAll(oldId, newId);
			edgeString2 = edgeString2.replaceAll(oldId, newId);
			basePairString2 = basePairString2.replaceAll(oldId, newId);
			nodeLabelString2 = nodeLabelString2.replaceAll(oldId, newId);
		}

		// build all nodes first
		sb.append(nodeString1);
		sb.append("|");
		sb.append(nodeString2);
		sb.append("$");

		// build all edges
		sb.append(edgeString1);
		if (edgeString1.length() > 0) {
			sb.append("|");
		}
		sb.append(edgeString2);
		sb.append("$");

		// build all base pairs
		sb.append(basePairString1);
		if (basePairString1.length() > 0) {
			sb.append("|");
		}
		sb.append(basePairString2);
		sb.append("$");

		// build all node labels
		sb.append(nodeLabelString1);
		if (nodeLabelString1.length() > 0) {
			sb.append("|");
		}
		sb.append(nodeLabelString2);
		sb.append("$");

		return sb.toString();
	}

	/**
	 * This method will automatically add base pair info into notation only if
	 * it contains TWO RNA polymer nodes and there is no base pairing info
	 * 
	 * @param complexNotation
	 * @return complex notation with base pairing
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.StructureException
	 */
	public static String hybridize(String complexNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		String result = null;
		List<RNAPolymerNode> l = getRNAPolymerNodeList(complexNotation);
		String basePairString = getAllBasePairString(complexNotation);
		if (l.size() == 2 && basePairString.length() == 0) {
			String nodeString = getAllNodeString(complexNotation);
			String edgeString = getAllEdgeString(complexNotation);
			String labelString = getAllNodeLabelString(complexNotation);
			basePairString = getBasePairString(l.get(0), l.get(1));
			result = nodeString + "$" + edgeString + "$" + basePairString + "$"
					+ labelString + "$";
		} else {
			result = complexNotation;
		}
		return result;
	}

	private static String getBasePairString(RNAPolymerNode node1,
			RNAPolymerNode node2) throws NotationException, IOException,
			JDOMException, MonomerException, StructureException {
		String basePair = "";
		String seq1 = node1.getSequence().replaceAll("T", "U");
		String seq2 = node2.getSequence().replaceAll("T", "U");
		char[] chars = seq2.toCharArray();
		StringBuffer sb = new StringBuffer();
		for (int i = chars.length; i > 0; i--) {
			String symbol = String.valueOf(chars[i - 1]);
			String compSymbol = NucleotideSequenceParser.complementMap
					.get(symbol);
			sb.append(compSymbol);
		}
		String compSeq2 = sb.toString();
		String maxSeqMatch = NucleotideSequenceParser.getMaxMatchFragment(seq1,
				compSeq2);
		int seqMatchLength = maxSeqMatch.length();
		int seq1NucStart = -1;
		int seq1MonomerStart = 0;
		int seq2NucStart = -1;
		int seq2MonomerStart = 0;

		List<Nucleotide> seq1NucList = SimpleNotationParser.getNucleotideList(
				node1.getLabel(), false);
		List<Nucleotide> seq2NucList = SimpleNotationParser.getNucleotideList(
				node2.getLabel(), false);

		if (seqMatchLength > 0) {
			// get the starting monomer position for sequence 1
			seq1NucStart = seq1.indexOf(maxSeqMatch);
			for (int i = 0; i < seq1NucStart; i++) {
				Nucleotide nuc = seq1NucList.get(i);
				int monomerCount = SimpleNotationParser
						.getMonomerCountForRNA(nuc.getNotation());
				seq1MonomerStart = seq1MonomerStart + monomerCount;
			}

			// get the starting monomer position for sequence 2
			int compSeq2NucStart = compSeq2.indexOf(maxSeqMatch);
			seq2NucStart = seq2.length() - seqMatchLength - compSeq2NucStart;
			for (int i = 0; i < seq2NucStart; i++) {
				Nucleotide nuc = seq2NucList.get(i);
				int monomerCount = SimpleNotationParser
						.getMonomerCountForRNA(nuc.getNotation());
				seq2MonomerStart = seq2MonomerStart + monomerCount;
			}

			// build the matching monomer position
			for (int i = 0; i < seqMatchLength; i++) {
				Nucleotide nuc1 = seq1NucList.get(i + seq1NucStart);
				if (null == nuc1.getBaseMonomer()) {
					throw new NotationException(
							"Nucleotide without base cannot be hybridized with others");
				}
				if (i == 0) {
					seq1MonomerStart = seq1MonomerStart + 2;
				} else {
					seq1MonomerStart = seq1MonomerStart + 3;
				}

				Nucleotide nuc2 = seq2NucList.get(i + seq2NucStart);
				if (null == nuc2.getBaseMonomer()) {
					throw new NotationException(
							"Nucleotide without base cannot be hybridized with others");
				}
				if (i == 0) {
					seq2MonomerStart = seq2MonomerStart + 2;
				} else {
					seq2MonomerStart = seq2MonomerStart + 3;
				}
			}

			// what if there is an X with two monomers in the middle? exception
			// should have been thrown
			// now build the base pair string, offset by 3
			for (int i = seqMatchLength; i > 0; i--) {
				int seq1MonomerPos = seq1MonomerStart - (i - 1) * 3;
				int seq2MonomerPos = seq2MonomerStart - (seqMatchLength - i)
						* 3;

				if (basePair.length() > 0) {
					basePair = basePair + "|";
				}
				basePair = basePair + "RNA1,RNA2," + seq1MonomerPos + ":pair-"
						+ seq2MonomerPos + ":pair";
			}
		}
		return basePair;
	}

	/**
	 * This method decomposes complex notation into complex notations of
	 * covalently connected monomers, base pairing info is lost
	 * 
	 * @param complexNotation
	 * @return complex notation for connected polymers
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.MonomerException
	 * @throws org.jdom.JDOMException
	 * @throws java.io.IOException
	 */
	public static String[] decompose(String complexNotation)
			throws NotationException, MonomerException, JDOMException,
			IOException {
		ComplexPolymer cp = parse(complexNotation);
		Map<String, String> annotationMap = cp.getPolymerNodeAnnotationMap();

		Map<Integer, List<String>> groupedNodes = getGroupNodeListMap(
				cp.getPolymerNodeList(), cp.getPolymerEdgeList());
		if (groupedNodes.size() == 1) {
			return new String[] { complexNotation };
		}

		Map<Integer, String> groupNotationMap = new HashMap<Integer, String>();
		Set groupSet = groupedNodes.keySet();
		for (Iterator i = groupSet.iterator(); i.hasNext();) {
			Integer groupId = (Integer) i.next();
			List<String> allNodes = groupedNodes.get(groupId);
			List<String> connectedNodes = new ArrayList<String>();
			// remove duplicated nodes from list
			for (String node : allNodes) {
				if (!connectedNodes.contains(node)) {
					connectedNodes.add(node);
				}
			}

			List<PolymerEdge> connectedEdges = getConnectedEdgeList(
					connectedNodes, cp.getPolymerEdgeList());

			// Is this sorted based on order from complexNotation?
			StringBuffer nodeSB = new StringBuffer();
			for (String nodeId : connectedNodes) {
				String polymerNotation = getPolymerNotation(nodeId,
						cp.getPolymerNodeList());
				if (nodeSB.length() > 0) {
					nodeSB.append("|");
				}
				nodeSB.append(nodeId + "{" + polymerNotation + "}");
			}

			StringBuffer edgeSB = new StringBuffer();
			for (PolymerEdge edge : connectedEdges) {
				String edgeString = edge.toString();
				if (edgeSB.length() > 0) {
					edgeSB.append("|");
				}
				edgeSB.append(edgeString);
			}

			StringBuffer annotationSB = new StringBuffer();
			for (String nodeId : connectedNodes) {
				if (annotationMap.containsKey(nodeId)) {
					if (annotationSB.length() > 0) {
						annotationSB.append("|");
					}
					annotationSB.append(nodeId + "{"
							+ annotationMap.get(nodeId) + "}");
				}
			}

			String result = nodeSB.toString() + "$" + edgeSB.toString() + "$$"
					+ annotationSB.toString() + "$";
			groupNotationMap.put(groupId, result);
		}

		// sort based on polymerNode list
		List<String> results = new ArrayList<String>();
		List<String> usedNodes = new ArrayList<String>();
		Set<Entry<Integer, List<String>>> entrySet = groupedNodes.entrySet();
		for (PolymerNode polymerNode : cp.getPolymerNodeList()) {
			String nodeId = polymerNode.getId();
			if (!usedNodes.contains(nodeId)) {
				for (Entry<Integer, List<String>> entry : entrySet) {
					Integer groupId = entry.getKey();
					List<String> nodes = entry.getValue();
					if (nodes.contains(nodeId)) {
						String notation = groupNotationMap.get(groupId);
						results.add(notation);
						usedNodes.addAll(nodes);
						break;
					}
				}

			}

		}

		return results.toArray(new String[0]);
	}

	/**
	 * This method replace existing monomer with new monomer for a given polymer
	 * type in the complex notation
	 * 
	 * @param complexNotation
	 * @param polymerType
	 * @param existingMonomerID
	 * @param newMonomerID
	 * @return complex notation after replacement
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.NotationException
	 */

	public static String replaceMonomer(String complexNotation,
			String polymerType, String existingMonomerID, String newMonomerID)
			throws MonomerException, IOException, JDOMException,
			NotationException {

		return replaceMonomer(complexNotation, polymerType, existingMonomerID,
				newMonomerID, null, true);
	}

	public static String replaceMonomer(String complexNotation,
			String polymerType, String existingMonomerID, String newMonomerID,
			MonomerStore monomerStore, boolean validate)
			throws MonomerException, IOException, JDOMException,
			NotationException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		if (validate)
			SimpleNotationParser.validateMonomerReplacement(polymerType,
					existingMonomerID, newMonomerID, monomerStoreToUse);
		if (null == complexNotation || complexNotation.length() == 0) {
			return complexNotation;
		}

		if (existingMonomerID.equals(newMonomerID)) {
			return complexNotation;
		}

		String allNodeString = getAllNodeString(complexNotation);
		String restOfNotation = complexNotation.substring(allNodeString
				.length());
		List<PolymerNode> polymers = getPolymerNodeList(allNodeString);
		StringBuffer sb = new StringBuffer();
		for (PolymerNode polymer : polymers) {
			if (sb.length() > 0) {
				sb.append(ComplexNotationParser.LIST_LEVEL_DELIMITER);
			}
			sb.append(polymer.getId());
			sb.append(ComplexNotationParser.NODE_LABEL_START_SYMBOL);

			if (polymer.getType().equals(polymerType)) {
				String notation = polymer.getLabel();
				String result = SimpleNotationParser.replaceMonomer(notation,
						polymerType, existingMonomerID, newMonomerID,
						monomerStoreToUse, validate);
				sb.append(result);
			} else {
				sb.append(polymer.getLabel());
			}
			sb.append(ComplexNotationParser.NODE_LABEL_END_SYMBOL);
		}
		return sb.toString() + restOfNotation;
	}

	/**
	 * Standardize all forms of HELM notation into complex notation
	 * 
	 * @param notation
	 *            - in the format of Complex Notation or Simple Notation for
	 *            RNA, Peptide, and Chem
	 * @return complex notation
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 * @throws org.helm.notation.NotationException
	 * @throws org.helm.notation.StructureException
	 */
	public static String standardize(String notation) throws MonomerException,
			IOException, JDOMException, NotationException, StructureException {
		MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
		return standardize(notation, factory.getMonomerStore());

	}

	public static String standardize(String notation, MonomerStore monomerStore)
			throws MonomerException, IOException, JDOMException,
			NotationException, StructureException {
		if (null == notation || notation.length() == 0) {
			return notation;
		}

		// convert to complex notation from simple notation
		if (notation.indexOf(ComplexNotationParser.TOP_LEVEL_DELIMITER) >= 0) {
			ComplexNotationParser.validateComplexNotation(notation,
					monomerStore);
			return notation;
		} else {
			Set<String> polymerTypes = MonomerFactory.getInstance()
					.getMonomerDB().keySet();
			for (String polymerType : polymerTypes) {
				try {
					SimpleNotationParser.validateSimpleNotation(notation,
							polymerType);
					return SimpleNotationParser.getComplexNotation(notation,
							polymerType);
				} catch (Exception e) {
				}
			}
			throw new NotationException(
					"Invalid HELM notation for simple polymer: " + notation);
		}
	}

	/**
	 * Returns true if the any of the polymerNodes that are not nucleotides have
	 * has a modified nucleotide. The assertion is based on breaking the
	 * notation into single Nucleotides and using the Nucleotide.isModified()
	 * 
	 * @param polymerNodes
	 * @return true or false
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */

	public static boolean hasNucleotideModification(
			List<PolymerNode> polymerNodes) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {

		return hasNucleotideModification(polymerNodes, null);
	}

	public static boolean hasNucleotideModification(
			List<PolymerNode> polymerNodes, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		for (PolymerNode node : polymerNodes) {
			String oneNodeID = node.getId();
			String notation = node.getLabel();

			if (oneNodeID.startsWith(Monomer.CHEMICAL_POLYMER_TYPE)) {
			} else if (oneNodeID.startsWith(Monomer.PEPTIDE_POLYMER_TYPE)) {
			} else {
				List<Nucleotide> nucleotidelist = SimpleNotationParser
						.getStrictNucleotideList(notation, false,
								monomerStoreToUse);
				for (Nucleotide oneNucleotide : nucleotidelist) {
					if (oneNucleotide.isModified()) {
						return true;
					}
				}
			}
		}
		return false;
	}

	public static int getTotalMonomerCount(String notation)
			throws NotationException {

		return getTotalMonomerCount(notation, null);
	}

	public static int getTotalMonomerCount(String notation,
			MonomerStore monomerStore) throws NotationException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		int totalMonomerCount = 0;

		List<PolymerNode> nodes = ComplexNotationParser
				.getPolymerNodeList(notation, monomerStoreToUse);
		for (PolymerNode node : nodes) {
			String polymerType = node.getType();
			String label = node.getLabel();
			int monomerCount = SimpleNotationParser.getMonomerCount(label,
					polymerType, monomerStoreToUse);
			totalMonomerCount = totalMonomerCount + monomerCount;
		}
		return totalMonomerCount;
	}

	/**
	 * Perform calculation without validation
	 * 
	 * @param extendedNotation
	 *            - complex polymer notation
	 * @return MoleculeInfo of complex polymer
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws PluginException
	 * @throws StructureException
	 */

	public static MoleculeInfo getMoleculeInfo(String extendedNotation)
			throws NotationException, MonomerException, IOException,
			JDOMException, PluginException, StructureException {

		return getMoleculeInfo(extendedNotation, null);
	}

	public static MoleculeInfo getMoleculeInfo(String extendedNotation,
			MonomerStore monomerStore) throws NotationException,
			MonomerException, IOException, JDOMException, PluginException,
			StructureException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		return getMoleculeInfo(extendedNotation, false, monomerStoreToUse);
	}

	/**
	 * This method returns the MoleculeInfo of input polymer notation using the
	 * divide and conquer approach in SimpleNotationParser Should return the
	 * same result as StructureParser.getMoleculeInfo(String smiles) method, but
	 * runs faster for standard edges Ignore fuzzy edges in the calculation,
	 * good enough for large structures
	 * 
	 * @param extendedNotation
	 * @param includeValidation
	 * @return MoleculeInfo object
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws PluginException
	 * @throws StructureException
	 */
	public static MoleculeInfo getMoleculeInfo(String extendedNotation,
			boolean includeValidation) throws NotationException,
			MonomerException, IOException, JDOMException, PluginException,
			StructureException {
		return getMoleculeInfo(extendedNotation, includeValidation, null);
	}

	public static MoleculeInfo getMoleculeInfo(String extendedNotation,
			boolean includeValidation, MonomerStore monomerStore)
			throws NotationException, MonomerException, IOException,
			JDOMException, PluginException, StructureException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		ComplexPolymer cp = parse(extendedNotation);

		if (includeValidation) {
			validateComplexPolymer(cp, monomerStoreToUse);
		}

		List<PolymerNode> nodeList = cp.getPolymerNodeList();
		List<PolymerEdge> edgeList = cp.getPolymerEdgeList();

		// deal with all connections between polymers
		List<MoleculeInfo> caps = new ArrayList<MoleculeInfo>();
		for (PolymerEdge edge : edgeList) {
			if (edge.getEdgeType() == PolymerEdge.GENERIC_EDGE) {
				// ignore generic connections
				continue;
			}

			String sourceId = edge.getSourceNode();
			PolymerNode sourceNode = getPolymerNodeFromID(nodeList, sourceId);
			int sourceMonomerNumber = edge.getSourceMonomerNumber();
			String sourceR = edge.getSourceR();
			MoleculeInfo sourceCapMI = getCapMoleculeInfo(sourceNode,
					sourceMonomerNumber, sourceR);
			caps.add(sourceCapMI);

			String targetId = edge.getTargetNode();
			PolymerNode targetNode = getPolymerNodeFromID(nodeList, targetId);
			int targetMonomerNumber = edge.getTargetMonomerNumber();
			String targetR = edge.getTargetR();
			MoleculeInfo targetCapMI = getCapMoleculeInfo(targetNode,
					targetMonomerNumber, targetR);
			caps.add(targetCapMI);
		}

		// deal with all simple polymers
		List<MoleculeInfo> chunks = new ArrayList<MoleculeInfo>();
		for (PolymerNode node : nodeList) {
			String polymerType = node.getType();
			String label = node.getLabel();
			MoleculeInfo tmpMI = SimpleNotationParser.getMoleculeInfo(label,
					polymerType, monomerStoreToUse);
			chunks.add(tmpMI);
		}

		return StructureParser.processMoleculeInfo(chunks, caps);
	}

	private static MoleculeInfo getCapMoleculeInfo(PolymerNode node,
			int monomerNumber, String rgroup) throws NotationException,
			MonomerException, IOException, JDOMException, PluginException {
		String polymerType = node.getType();
		String notation = node.getLabel();

		List<String> monomerIdList = SimpleNotationParser.getMonomerIDList(
				notation, polymerType);
		Monomer monomer = SimpleNotationParser.getMonomer(
				monomerIdList.get(monomerNumber - 1), polymerType);
		return monomer.getCapMoleculeInfo(rgroup);
	}

	private static PolymerNode getPolymerNodeFromID(List<PolymerNode> nodes,
			String nodeID) {
		for (PolymerNode node : nodes) {
			if (node.getId().equals(nodeID)) {
				return node;
			}
		}
		return null;
	}

	/*
	 * private static Map<String, RgroupStructure>
	 * getPolymerNodeStructureMap(List<PolymerNode> nodeList) throws
	 * IOException, NotationException, MonomerException, StructureException,
	 * JDOMException { MonomerFactory factory = null; try { factory =
	 * MonomerFactory.getInstance(); } catch (Exception ex) { throw new
	 * NotationException("Unable to initialize monomer factory", ex); } return
	 * getPolymerNodeStructureMap(nodeList,factory.getMonomerStore()); }
	 */

	private static Map<String, RgroupStructure> getPolymerNodeStructureMap(
			List<PolymerNode> nodeList, MonomerStore monomerStore)
			throws IOException, NotationException, MonomerException,
			StructureException, JDOMException {

		MonomerStore monomerStoreToUse = checkForMonomerStore(monomerStore);

		Map<String, RgroupStructure> nodeStrucMap = new HashMap<String, RgroupStructure>();
		for (int i = 0; i < nodeList.size(); i++) {
			PolymerNode node = nodeList.get(i);
			String nodeId = node.getId();
			String nodeLabel = node.getLabel();
			String polymerType = PolymerNode.getPolymerType(nodeId);
			RgroupStructure struc = SimpleNotationParser
					.getSimplePolymerStructure(nodeLabel, polymerType,
							monomerStoreToUse);
			if (null == struc.getMolecule()) {
				throw new NotationException(
						"Polymer notation contains non-specific monomer structure");
			}
			nodeStrucMap.put(nodeId, struc);
		}
		return nodeStrucMap;
	}
	
	/*
	 * This function replaces smiles in complex notation with temporary ids
	 * 
	 */
	public static String getNotationByReplacingSmiles(String helmString,MonomerStore monomerStore) throws NotationException, MonomerException, JDOMException, IOException{
		
		String allNodeString = getAllNodeString(helmString);
		String restOfNotation = helmString.substring(allNodeString
				.length());
		List<PolymerNode> polymers = getPolymerNodeList(allNodeString);
		StringBuffer sb = new StringBuffer();
		for (PolymerNode polymer : polymers) {
			if (sb.length() > 0) {
				sb.append(ComplexNotationParser.LIST_LEVEL_DELIMITER);
			}
			sb.append(polymer.getId());
			sb.append(ComplexNotationParser.NODE_LABEL_START_SYMBOL);

			String notation = polymer.getLabel();
			String result = SimpleNotationParser.getNotationByReplacingSmiles(notation, polymer.getType(), monomerStore);
			sb.append(result);
			
			sb.append(ComplexNotationParser.NODE_LABEL_END_SYMBOL);
		}
		return sb.toString() + restOfNotation;
		
		
		
		
	}
}
