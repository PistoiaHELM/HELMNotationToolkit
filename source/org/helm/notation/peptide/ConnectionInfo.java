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

/**
 * 
 * @author zhangtianhong
 */
class ConnectionInfo {

	private int listIndex;
	private int monomerIndex;
	private String rgroup;
	private String monomerID;

	public ConnectionInfo() {
	}

	public ConnectionInfo(int listIndex, int monomerIndex, String monomerID,
			String rgroup) {
		this.listIndex = listIndex;
		this.monomerIndex = monomerIndex;
		this.monomerID = monomerID;
		this.rgroup = rgroup;
	}

	public int getListIndex() {
		return listIndex;
	}

	public void setListIndex(int listIndex) {
		this.listIndex = listIndex;
	}

	public String getMonomerID() {
		return monomerID;
	}

	public void setMonomerID(String monomerID) {
		this.monomerID = monomerID;
	}

	public int getMonomerIndex() {
		return monomerIndex;
	}

	public void setMonomerIndex(int monomerIndex) {
		this.monomerIndex = monomerIndex;
	}

	public String getRgroup() {
		return rgroup;
	}

	public void setRgroup(String rgroup) {
		this.rgroup = rgroup;
	}
}
