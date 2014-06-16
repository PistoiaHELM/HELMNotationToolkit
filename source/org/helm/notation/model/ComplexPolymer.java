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

import java.util.List;
import java.util.Map;

/**
 * This is the data model for ComplexPolymer. It contains the lists for
 * PolymerNode, connected PolymerEdge, and base-paired PolymerEdge.
 * 
 * @author zhangtianhong
 */
public class ComplexPolymer {
	private List<PolymerNode> polymerNodeList;
	private List<PolymerEdge> polymerEdgeList;
	private List<PolymerEdge> basePairList;
	private Map<String, String> polymerNodeAnnotationMap;

	public List<PolymerNode> getPolymerNodeList() {
		return polymerNodeList;
	}

	public void setPolymerNodeList(List<PolymerNode> polymerNodeList) {
		this.polymerNodeList = polymerNodeList;
	}

	public List<PolymerEdge> getPolymerEdgeList() {
		return polymerEdgeList;
	}

	public void setPolymerEdgeList(List<PolymerEdge> polymerEdgeList) {
		this.polymerEdgeList = polymerEdgeList;
	}

	public List<PolymerEdge> getBasePairList() {
		return basePairList;
	}

	public void setBasePairList(List<PolymerEdge> basePairList) {
		this.basePairList = basePairList;
	}

	/**
	 * @return the polymerNodeAnnotationMap
	 */
	public Map<String, String> getPolymerNodeAnnotationMap() {
		return polymerNodeAnnotationMap;
	}

	/**
	 * @param polymerNodeAnnotationMap
	 *            the polymerNodeAnnotationMap to set
	 */
	public void setPolymerNodeAnnotationMap(
			Map<String, String> polymerNodeAnnotationMap) {
		this.polymerNodeAnnotationMap = polymerNodeAnnotationMap;
	}
}
