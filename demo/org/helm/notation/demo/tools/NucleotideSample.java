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

import org.helm.notation.tools.*;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import java.util.Map;

/**
 * 
 * @author zhangtianhong
 */
public class NucleotideSample {

	public static void main(String[] args) {
		try {
			// [H]OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)[nH]1cnc2c1nc(N)[nH]c2=O)[nH]1ccc(=O)[nH]c1=O)[nH]1cnc2c(N)ncnc12

			String notation = NucleotideSequenceParser.getNotation("AUG");
			System.out.println(notation);

			System.out.println(SimpleNotationParser.getSimplePolymerSMILES(
					notation, "RNA"));

			String complexNotation = "RNA1{" + notation + "}$$$$";
			String smiles = ComplexNotationParser
					.getComplexPolymerSMILES(complexNotation);
			System.out.println(smiles);

			boolean opp = NucleotideSequenceParser.isInOppositeDirection(
					"AUGTAU", "AUACAU");
			if (opp) {
				System.out.println("opposite direction");
			} else {
				System.out.println("same direction");
			}

			opp = NucleotideSequenceParser.isInOppositeDirection("ATGCAAA",
					"UUUGCAU");
			if (opp) {
				System.out.println("opposite direction");
			} else {
				System.out.println("same direction");
			}

			String senseSeq = "CGAAAUGUUCAUACUGUUGdTdT";
			String antiSenseSeq = "UUACAAUUUGGACUUUCCGdTdT";

			String siNotation = NucleotideSequenceParser.getSirnaNotation(
					senseSeq, antiSenseSeq);
			System.out.println(siNotation);

			siNotation = NucleotideSequenceParser.getSirnaNotation(senseSeq,
					antiSenseSeq,
					NucleotideSequenceParser.RNA_DESIGN_TUSCHL_19_PLUS_2);
			System.out.println(siNotation);

			// String[] sequences = new
			// String[]{"5'-UAU GUC UCC AGA AUG UAG CdTdT-3'",
			// "5'-mGCU ACA UUC UGG AGA CAU AdTdT-3'",
			// "agcuagggu"
			// };
			//
			// for (int i = 0; i < sequences.length; i++) {
			// String sequence = sequences[i];
			// System.out.println("In Sequence: " + sequence);
			// System.out.println("Complement Sequence(Normal): " +
			// NucleotideSequenceParser.getNormalComplementSequence(sequence));
			// System.out.println("Complement Sequence(Reverse): " +
			// NucleotideSequenceParser.getReverseComplementSequence(sequence));
			// System.out.println("Notation: " +
			// NucleotideSequenceParser.getNotation(sequence));
			// }
			//
			// Map<String, Map<String, String>> templates =
			// NucleotideFactory.getInstance().getNucleotideTemplates();
			// String templatesXML =
			// NucleotideSequenceParser.getNucleotideTemplatesXML(templates);
			// System.out.println(templatesXML);
			// NucleotideFactory.getInstance().saveNucleotideTemplates();

			notation = "[LR]([5meC])[sP]";
			Nucleotide n = new Nucleotide("ls5C", notation);
			System.out.println(n.getNaturalAnalog());

			String base = notation.substring(notation.indexOf("(") + 1,
					notation.indexOf(")"));
			base = base.replaceAll("\\[|\\]", "");
			Monomer baseMonomer = MonomerFactory.getInstance().getMonomerDB()
					.get(Monomer.NUCLIEC_ACID_POLYMER_TYPE).get(base);
			System.out.println(baseMonomer.getNaturalAnalog());

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
