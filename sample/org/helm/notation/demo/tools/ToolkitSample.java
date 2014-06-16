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

import chemaxon.marvin.plugin.PluginException;
import org.helm.notation.MonomerException;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.tools.*;
import org.helm.notation.MonomerFactory;
import org.helm.notation.model.MoleculeInfo;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jdom.JDOMException;

/**
 * This is an example on how to use the toolkit API for validation
 * 
 * @author ZHANGTIANHONG
 */
public class ToolkitSample {

	public static void main(String[] args) {
		String seq;
		String notation;
		String smiles;
		MoleculeInfo mi;

		try {

			// String testSmiles =
			// "[H]OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)[nH]1cnc2c1nc(N)[nH]c2=O)[nH]1ccc(=O)[nH]c1=O)[nH]1cnc2c(N)ncnc12";
			// boolean isValid = StructureParser.validateSmiles(testSmiles);
			//
			// //SMILES validation
			// String validSmiles =
			// "CC(C)[C@H](N[*])C([*])=O |r,$;;;;;_R1;;_R2;$|";
			// String invalidSmiles = "Cc7ccccc7Cc1(c7ccccc7)nnnc2cc(C(=O))c21";
			//
			// if (StructureParser.validateSmiles(validSmiles)) {
			// System.out.println("Valid SMILES is valid");
			// } else {
			// System.out.println("Valid SMILES is invalid");
			// }
			// if (StructureParser.validateSmiles(invalidSmiles)) {
			// System.out.println("Invalid SMILES is valid");
			// } else {
			// System.out.println("invalid SMILES is invalid");
			// }
			//
			// MoleculeInfo mi = StructureParser.getMoleculeInfo(validSmiles);
			// System.out.println("MF: " + mi.getMolecularFormula());
			// System.out.println("MW: " + mi.getMolecularWeight());

			// Complex polymer notation validation and conversion to SMILES
			// merge routine can be optimized for long sequence, maybe broken
			// into small chunks, rather than walk the sequence from begining to
			// end?
			// 21mer + 21mer ~2s 21mer ~1s
			// String notation =
			// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P}"
			// +
			// "|RNA2{R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P}$$$$";

			// 21mer + 21mer + 21mer + 21mer ~15s
			// String notation =
			// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P}"+
			// "|RNA2{R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P}"+
			// "|RNA3{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P}"+
			// "|RNA4{R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P}";

			// 42mer + 42mer ~10s 42mer ~4s
			// String notation =
			// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P.R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P}";
			// "|RNA2{R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P.R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P}";

			// 84mer + 84mer ~80s 84mer ~30s
			// String notation =
			// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P.R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P.R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P.R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P}";
			// "|RNA2{R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P.R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P.R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P.R(G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P}";

			// siRNA with Alexa dye
			// String notation =
			// "RNA1{R(G)P.R(G)P.R(U)P.R(C)P.R(G)P.R(G)P.R(A)P.R(G)P.R(U)P.R(C)P.R(A)P.R(A)P.R(C)P.R(G)P.R(G)P.R(A)P.R(U)P.R(U)P.R(U)P.[dR](T)P.[dR](T)P}"+
			// "|RNA2{R(A)P.R(A)P.R(A)P.R(U)P.R(C)P.R(C)P.R(G)P.R(U)P.R(U)P.R(G)P.R(A)P.R(C)P.R(U)P.R(C)P.R(C)P.R(G)P.R(A)P.R(C)P.R(C)P.[dR](T)P.[dR](T)P}"+
			// "|CHEM1{Alexa}$RNA1,CHEM1,-$$$";

			seq = "AGAUGCACAACAAGAUCUACG";
			// seq = "A";
			System.out.println("Time: " + System.currentTimeMillis());
			notation = NucleotideSequenceParser.getSirnaNotation(seq, null);
			System.out.println("Noation: " + notation);

			// 680 long RNA
			notation = "RNA1{R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)}$$$$";
			System.out.println("Time: " + System.currentTimeMillis());
			if (ComplexNotationParser.validateComplexNotation(notation)) {
				System.out.println("Notation is valid");
			} else {
				System.out.println("Notation is invalid");
			}

			// cyclic RNA
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P.R(A)P.R(A)P}$RNA1,RNA1,27:R2-1:R1$$$";
			testMoleculeProperty(notation);

			// linear RNA
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P.R(A)P.R(A)P}$$$$";
			testMoleculeProperty(notation);

			// cyclic peptide
			notation = "PEPTIDE1{A.G.G.K.K}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$";
			testMoleculeProperty(notation);

			// noncyclic peptide
			notation = "PEPTIDE1{A.G.G.K.K}$$$$";
			testMoleculeProperty(notation);

			// siRNA
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|CHEM1{A6P}$RNA1,CHEM1,1:R1-1:R2$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,8:pair-20:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,20:pair-8:pair$RNA3{ss}$";
			testDecompose(notation);

			// cycle + linear mixture
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}$RNA1,RNA1,21:R2-1:R1$$$"; // cycle
																																					// +
																																					// linear
																																					// mixture
			testDecompose(notation);

			// multi connections
			notation = "RNA1{R(A)}|RNA2{R(A)}|CHEM1{CovX-3}|CHEM2{CovX-2}$RNA1,CHEM2,1:R1-1:R1|RNA1,CHEM1,1:R2-1:R1$$$";
			testDecompose(notation);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void testMoleculeProperty(String notation) {
		try {
			System.out.println("Time: " + System.currentTimeMillis());
			String smiles = ComplexNotationParser
					.getComplexPolymerSMILES(notation);
			System.out.println("SMILES: " + smiles);
			System.out.println("Time: " + System.currentTimeMillis());
			MoleculeInfo mi = StructureParser.getMoleculeInfo(smiles);
			System.out.println("MF: " + mi.getMolecularFormula());
			System.out.println("MW: " + mi.getMolecularWeight());
			System.out.println("Time: " + System.currentTimeMillis());
		} catch (PluginException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (IOException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (NotationException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (MonomerException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (StructureException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (JDOMException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		}
	}

	private static void testDecompose(String notation) {
		System.out.println("input: ");
		System.out.println(notation);
		System.out.println("output:");
		try {
			String[] notations = ComplexNotationParser.decompose(notation);
			for (int i = 0; i < notations.length; i++) {
				System.out.println(notations[i]);
			}
		} catch (NotationException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (MonomerException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (JDOMException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (IOException ex) {
			Logger.getLogger(ToolkitSample.class.getName()).log(Level.SEVERE,
					null, ex);
		}
	}
}
