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
package org.helm.notation.demo.tools;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

public class UniqueSmilesTest {

	public static void main(String[] args) {

		try {

			String smiles1 = "C[C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|";
			toUniqueSmiles(smiles1);

			String smiles2 = "C[C@H](N[*])C([*])=O";
			toUniqueSmiles(smiles2);

			String smiles3 = "C[C@H](N*)C(*)=O";
			toUniqueSmiles(smiles3);

			// monomer
			String smiles5 = "[*]NCC([*])=O |$_R1;;;;_R2;$|";
			toCXSmiles(smiles5);

			String smiles4 = "[*]C(=O)CN[*] |$_R1;;;;;_R2$|";
			toCXSmiles(smiles4);

			String smiles6 = "[*]C(=O)CN[*] |$_R4;;;;;_R3$|";
			toCXSmiles(smiles6);

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static void toCXSmiles(String smiles) throws IOException {
		System.out.println("Input:\t" + smiles);
		InputStream is = new ByteArrayInputStream(smiles.getBytes());
		MolImporter importer = new MolImporter(is);
		Molecule molecule = importer.read();

		String smilesU = molecule.toFormat("cxsmiles:u");
		System.out.println("Output:\t" + smilesU);
	}

	private static void toUniqueSmiles(String smiles) throws IOException {
		System.out.println("Input:\t" + smiles);
		InputStream is = new ByteArrayInputStream(smiles.getBytes());
		MolImporter importer = new MolImporter(is, "smiles");
		Molecule molecule = importer.read();

		String smilesU = molecule.toFormat("smiles:u");
		System.out.println("Output:\t" + smilesU);
	}
}
