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

import org.helm.notation.NotationException;
import org.helm.notation.model.PolymerEdge;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Edge Notation in the format of:
 * SourcePolymerNode1[+SourcePolymerNode2+..],TargetPolymerNode1
 * [+TargetPolymerNode2+..],SourceConnection-TargetConnection where
 * SourceConnection and TargetConnection could be one of the three formats:
 * Standard Connection: MonomerPosition:RGroupID Hydrogen Bond:
 * MonomerPosition:pair Generic Connection: generic:
 * 
 * @author zhangtianhong
 */
public class EdgeParser {

	/**
	 * parses input edgeNotation into a PolymerEdge object
	 * 
	 * @param edgeNotation
	 * @return a PolymerEdge object
	 * @throws NotationException
	 */
	public static PolymerEdge parse(String edgeNotation)
			throws NotationException {
		PolymerEdge edge = new PolymerEdge();
		if (null == edgeNotation || edgeNotation.length() == 0) {
			throw new NotationException("Empty polymer edge notation");
		}

		String[] edgeComponents = edgeNotation.split(
				PolymerEdge.EDGE_COMPONENT_DELIMITER, -1);
		if (edgeComponents.length != 3) {
			throw new NotationException(
					"Invalid polymer edge notation: must contain three components");
		}

		String sourceNode = edgeComponents[0];
		if (checkNode(sourceNode)) {
			edge.setSourceNode(sourceNode);
		}

		String targetNode = edgeComponents[1];
		if (checkNode(targetNode)) {
			edge.setTargetNode(targetNode);
		}

		String connection = edgeComponents[2];
		PolymerEdge tmp = parseConnection(connection);
		edge.setConnection(connection);
		edge.setEdgeType(tmp.getEdgeType());
		edge.setSourceAttachmentType(tmp.getSourceAttachmentType());
		edge.setTargetAttachmentType(tmp.getTargetAttachmentType());

		edge.setEdgeNotation(edgeNotation);

		return edge;
	}

	/**
	 * canonicalize input polymerEdge based on sortedNodeIDs
	 * 
	 * @param polymerEdge
	 * @param sortedNodeIDs
	 */
	public static void canonicalize(PolymerEdge polymerEdge,
			List<String> sortedNodeIDs) throws NotationException {
		String[] sourceNodes = polymerEdge.getSourceNodes();
		if (sourceNodes.length > 1) {
			String newSourceNode = canonicalizeNode(sourceNodes, sortedNodeIDs);
			polymerEdge.setSourceNode(newSourceNode);
		}

		String[] targetNodes = polymerEdge.getTargetNodes();
		if (targetNodes.length > 1) {
			String newTargetNode = canonicalizeNode(targetNodes, sortedNodeIDs);
			polymerEdge.setTargetNode(newTargetNode);
		}

		int sourceIndex = getSourceNodeIndex(polymerEdge, sortedNodeIDs);
		int targetIndex = getTargetNodeIndex(polymerEdge, sortedNodeIDs);

		if (sourceIndex == targetIndex) {
			if (sourceNodes.length > 1 || targetNodes.length > 1) {
				sourceIndex = getSourceNodeMinIndex(polymerEdge, sortedNodeIDs);
				targetIndex = getTargetNodeMinIndex(polymerEdge, sortedNodeIDs);
			}
		}

		// reverse it if necessary
		if (sourceIndex > targetIndex) {
			String oldSourceNode = polymerEdge.getSourceNode();
			String oldTargetNode = polymerEdge.getTargetNode();
			polymerEdge.setSourceNode(oldTargetNode);
			polymerEdge.setTargetNode(oldSourceNode);
			polymerEdge.setConnection(polymerEdge.getReverseConnection());
		}
	}

	/**
	 * sort canonical polymer edges based on sorted nodes
	 * 
	 * @param canonicalPolymerEdges
	 * @param sortedNodeIDs
	 * @return sorted polymer edge list
	 * @throws NotationException
	 */
	public static List<PolymerEdge> sort(
			List<PolymerEdge> canonicalPolymerEdges, List<String> sortedNodeIDs)
			throws NotationException {
		List<String> sortedEdgeNodeList = getSortedEdgeNodeList(
				canonicalPolymerEdges, sortedNodeIDs);

		List<PolymerEdge> sortedEdgeList = new ArrayList<PolymerEdge>();
		List<String> sortedEdgeKeyList = new ArrayList<String>();
		Map<String, List<PolymerEdge>> groupEdgeMap = new HashMap<String, List<PolymerEdge>>();
		for (int i = 0; i < sortedEdgeNodeList.size(); i++) {
			String id = sortedEdgeNodeList.get(i);

			// allow self connection
			for (int j = i; j < sortedEdgeNodeList.size(); j++) {
				String nextId = sortedEdgeNodeList.get(j);
				String key = id + nextId;
				for (PolymerEdge edge : canonicalPolymerEdges) {

					if (edge.getSourceNode().equals(id)
							&& edge.getTargetNode().equals(nextId)) {
						PolymerEdge tmpEdge = new PolymerEdge();
						tmpEdge.setSourceNode(edge.getSourceNode());
						tmpEdge.setTargetNode(edge.getTargetNode());
						tmpEdge.setConnection(edge.getConnection());
						if (groupEdgeMap.containsKey(key)) {
							List<PolymerEdge> l = groupEdgeMap.get(key);
							l.add(tmpEdge);
						} else {
							List<PolymerEdge> l = new ArrayList<PolymerEdge>();
							l.add(tmpEdge);
							groupEdgeMap.put(key, l);
							sortedEdgeKeyList.add(key);
						}
					}
				}
			}
		}

		for (String key : sortedEdgeKeyList) {
			List<PolymerEdge> l = groupEdgeMap.get(key);
			Collections.sort(l, new PolymerEdgeComparator());
			sortedEdgeList.addAll(l);
		}
		return sortedEdgeList;
	}

	private static boolean checkNode(String node) throws NotationException {
		if (null == node || node.length() == 0) {
			throw new NotationException("Node is empty");
		}

		String[] nodeIDs = node.split(PolymerEdge.NODE_CONCATENATOR_REGEX, -1);
		if (nodeIDs.length > 1) {
			String type = getPolymerType(nodeIDs[0]);
			for (int i = 1; i < nodeIDs.length; i++) {
				String tmp = getPolymerType(nodeIDs[i]);
				if (!tmp.equals(type)) {
					throw new NotationException(
							"Group polymers must of the same type");
				}
			}
		}
		return true;
	}

	private static PolymerEdge parseConnection(String connection)
			throws NotationException {
		PolymerEdge pe = new PolymerEdge();
		pe.setConnection(connection);

		if (null == connection || connection.length() == 0) {
			throw new NotationException("Empty connection notation");
		}

		String[] attachments = connection.split(
				PolymerEdge.CONNECTION_SEPARATOR, -1);
		if (attachments.length != 2) {
			throw new NotationException(
					"Connection must contain two parts separated with "
							+ PolymerEdge.CONNECTION_SEPARATOR);
		}

		int srcAtt = parseAttachment(attachments[0]);
		pe.setSourceAttachmentType(srcAtt);

		int tgtAtt = parseAttachment(attachments[1]);
		pe.setTargetAttachmentType(tgtAtt);

		int edgeType = 0;
		if (srcAtt == PolymerEdge.STANDARD_EDGE_ATTACHMENT
				&& tgtAtt == PolymerEdge.STANDARD_EDGE_ATTACHMENT) {
			edgeType = PolymerEdge.STANDARD_EDGE;
		} else if (srcAtt == PolymerEdge.PAIR_EDGE_ATTACHMENT
				&& tgtAtt == PolymerEdge.PAIR_EDGE_ATTACHMENT) {
			edgeType = PolymerEdge.PAIR_EDGE;
		} else if (srcAtt == PolymerEdge.GENERIC_EDGE_ATTACHMENT
				|| tgtAtt == PolymerEdge.GENERIC_EDGE_ATTACHMENT) {
			edgeType = PolymerEdge.GENERIC_EDGE;
		} else {
			throw new NotationException("Invalid connection notation: "
					+ connection);
		}
		pe.setEdgeType(edgeType);

		return pe;
	}

	private static int parseAttachment(String attachment)
			throws NotationException {

		int attachmentType = 0;
		if (null == attachment && attachment.length() == 0) {
			throw new NotationException("Empty attachment notation");
		}

		String[] parts = attachment.split(
				PolymerEdge.MONOMER_ATTACHEMENT_SEPARATOR, -1);
		if (parts.length != 2) {
			throw new NotationException(
					"Attachment must contain two parts separated with "
							+ PolymerEdge.MONOMER_ATTACHEMENT_SEPARATOR);
		}

		if ((parts[0].matches("[0-9]++"))) {
			if (parts[1].startsWith("R")) {
				try {
					MonomerParser.validateAttachmentLabel(parts[1]);
					attachmentType = PolymerEdge.STANDARD_EDGE_ATTACHMENT;
				} catch (Exception e) {
				}
			} else if (parts[1].equalsIgnoreCase(PolymerEdge.PAIR_EDGE_KEY)) {
				attachmentType = PolymerEdge.PAIR_EDGE_ATTACHMENT;
			}
		} else if (parts[0].equalsIgnoreCase(PolymerEdge.GENERIC_EDGE_KEY)) {
			attachmentType = PolymerEdge.GENERIC_EDGE_ATTACHMENT;
		} else {
			throw new NotationException("Invalid attachment notation: "
					+ attachment);
		}

		return attachmentType;
	}

	private static String getPolymerType(String nodeID) {
		if (null == nodeID) {
			return null;
		}

		char[] chars = nodeID.toCharArray();
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < chars.length; i++) {
			if (String.valueOf(chars[i]).matches("[0-9]")) {
				break;
			} else {
				sb.append(chars[i]);
			}
		}
		return sb.toString();
	}

	private static List<String> getSortedEdgeNodeList(
			List<PolymerEdge> polymerEdges, List<String> sortedNodeIDs)
			throws NotationException {
		List<String> edgeNodeList = new ArrayList<String>();

		// sort on node index
		Map<Integer, List<String>> indexNodesMap = new TreeMap<Integer, List<String>>();
		for (PolymerEdge pe : polymerEdges) {
			String srcNode = pe.getSourceNode();
			int srcIndex = getSourceNodeIndex(pe, sortedNodeIDs);
			if (indexNodesMap.containsKey(srcIndex)) {
				List<String> nodes = indexNodesMap.get(srcIndex);
				if (!nodes.contains(srcNode)) {
					nodes.add(srcNode);
				}
			} else {
				List<String> nodes = new ArrayList<String>();
				nodes.add(srcNode);
				indexNodesMap.put(srcIndex, nodes);
			}

			String tarNode = pe.getTargetNode();
			int tarIndex = getTargetNodeIndex(pe, sortedNodeIDs);
			if (indexNodesMap.containsKey(tarIndex)) {
				List<String> nodes = indexNodesMap.get(tarIndex);
				if (!nodes.contains(tarNode)) {
					nodes.add(tarNode);
				}
			} else {
				List<String> nodes = new ArrayList<String>();
				nodes.add(tarNode);
				indexNodesMap.put(tarIndex, nodes);
			}
		}
		Set<Integer> indexSet = indexNodesMap.keySet();
		for (Iterator<Integer> it = indexSet.iterator(); it.hasNext();) {
			Integer index = it.next();
			List<String> nodes = indexNodesMap.get(index);
			if (nodes.size() > 1) {
				// sort on min node index
				Map<Integer, String> minIndexNodeMap = new TreeMap<Integer, String>();
				for (String node : nodes) {
					int minIndex = getMinNodeIndex(node, sortedNodeIDs);
					minIndexNodeMap.put(minIndex, node);
				}

				Set<Integer> minSet = minIndexNodeMap.keySet();
				for (Iterator<Integer> ite = minSet.iterator(); ite.hasNext();) {
					Integer min = ite.next();
					String node = minIndexNodeMap.get(min);
					edgeNodeList.add(node);
				}
			} else {
				edgeNodeList.addAll(nodes);
			}
		}

		return edgeNodeList;
	}

	private static String canonicalizeNode(String[] nodeIDs,
			List<String> sortedNodeIDs) throws NotationException {
		List<String> ids = Arrays.asList(nodeIDs);
		List<String> results = new ArrayList<String>();
		for (String id : sortedNodeIDs) {
			if (ids.contains(id)) {
				results.add(id);
			}
		}

		int delta = ids.size() - results.size();

		if (delta != 0) {
			throw new NotationException(
					"Invalid polymer node ID found in Generic connection");
		}

		StringBuilder sb = new StringBuilder();
		for (String result : results) {
			if (sb.length() > 0) {
				sb.append(PolymerEdge.NODE_CONCATENATOR);
			}
			sb.append(result);
		}

		return sb.toString();

	}

	private static int getTargetNodeIndex(PolymerEdge polymerEdge,
			List<String> sortedNodeIDs) throws NotationException {
		String[] nodes = polymerEdge.getTargetNodes();
		return getNodeIndex(nodes, sortedNodeIDs);
	}

	private static int getSourceNodeIndex(PolymerEdge polymerEdge,
			List<String> sortedNodeIDs) throws NotationException {
		String[] nodes = polymerEdge.getSourceNodes();
		return getNodeIndex(nodes, sortedNodeIDs);
	}

	private static int getNodeIndex(String[] nodes, List<String> sortedNodeIDs)
			throws NotationException {
		int index = 0;

		for (String node : nodes) {
			boolean found = false;
			for (int i = 0; i < sortedNodeIDs.size(); i++) {
				String tmp = sortedNodeIDs.get(i);
				if (tmp.equals(node)) {
					index = index + i;
					found = true;
					break;
				}
			}
			if (!found) {
				throw new NotationException(
						"Invalid polymer node ID found in Generic connection");
			}
		}

		return index;
	}

	private static int getNodeIndex(String node, List<String> sortedNodeIDs)
			throws NotationException {
		String[] nodes = node.split(PolymerEdge.NODE_CONCATENATOR_REGEX, -1);
		return getNodeIndex(nodes, sortedNodeIDs);
	}

	private static int getTargetNodeMinIndex(PolymerEdge polymerEdge,
			List<String> sortedNodeIDs) throws NotationException {
		String[] nodes = polymerEdge.getTargetNodes();
		return getMinNodeIndex(nodes, sortedNodeIDs);
	}

	private static int getSourceNodeMinIndex(PolymerEdge polymerEdge,
			List<String> sortedNodeIDs) throws NotationException {
		String[] nodes = polymerEdge.getSourceNodes();
		return getMinNodeIndex(nodes, sortedNodeIDs);
	}

	private static int getMinNodeIndex(String[] nodes,
			List<String> sortedNodeIDs) throws NotationException {
		int index = 1000;

		for (String node : nodes) {
			boolean found = false;
			for (int i = 0; i < sortedNodeIDs.size(); i++) {
				String tmp = sortedNodeIDs.get(i);
				if (tmp.equals(node)) {
					if (i < index) {
						index = i;
					}
					found = true;
					break;
				}
			}
			if (!found) {
				throw new NotationException(
						"Invalid polymer node ID found in Generic connection");
			}
		}

		return index;
	}

	private static int getMinNodeIndex(String node, List<String> sortedNodeIDs)
			throws NotationException {
		String[] nodes = node.split(PolymerEdge.NODE_CONCATENATOR_REGEX, -1);
		return getMinNodeIndex(nodes, sortedNodeIDs);
	}

	private static class PolymerEdgeComparator implements Comparator {

		public int compare(Object o1, Object o2) {
			String s1 = ((PolymerEdge) o1).toString();
			String s2 = ((PolymerEdge) o2).toString();
			return s1.compareTo(s2);
		}
	}
}
