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

import chemaxon.struc.Molecule;
import java.util.List;
import java.util.Map;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class PeptideFragment {

	private int id;
	private int level;
	private Molecule molecule;
	private List<List<String>> monomerIDsList;
	private List<List<Map<String, String>>> connectionMapsList;

	public PeptideFragment() {
	}

	public PeptideFragment(int level, Molecule molecule) {
		this.level = level;
		this.molecule = molecule;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public int getLevel() {
		return level;
	}

	public void setLevel(int level) {
		this.level = level;
	}

	public Molecule getMolecule() {
		return molecule;
	}

	public void setMolecule(Molecule molecule) {
		this.molecule = molecule;
	}

	public List<List<Map<String, String>>> getConnectionMapsList() {
		return connectionMapsList;
	}

	public void setConnectionMapList(
			List<List<Map<String, String>>> connectionMapsList) {
		this.connectionMapsList = connectionMapsList;
	}

	public List<List<String>> getMonomerIDsList() {
		return monomerIDsList;
	}

	public void setMonomerIDsList(List<List<String>> monomerIDsList) {
		this.monomerIDsList = monomerIDsList;
	}
}
