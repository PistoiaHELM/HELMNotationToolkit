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

/**
 * 
 * @author ZHANGTIANHONG
 */
public class CalculatorSample {

	public static void main(String[] args) {
		try {
			ExtinctionCoefficientCalculator calculator = ExtinctionCoefficientCalculator
					.getInstance();

			String input = "AGGDDDDDDDDDDDDDDDDDDFFFFFFFFFFFFF";
			float result = calculator.calculateFromAminoAcidSequence(input);
			System.out.println("AA Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "AGGCFFFFFFFFFF";
			result = calculator.calculateFromAminoAcidSequence(input);
			System.out.println("AA Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "AGGYEEEEEEEEEEEEEEEEEEE";
			result = calculator.calculateFromAminoAcidSequence(input);
			System.out.println("AA Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "AGGWEEEEEEEEEEEEEEEEEEE";
			result = calculator.calculateFromAminoAcidSequence(input);
			System.out.println("AA Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "A.G.G.W.E.E.E.E.E.W";
			result = calculator.calculateFromPeptidePolymerNotation(input);
			System.out.println("Peptide Notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$";
			result = calculator.calculateFromComplexNotation(input,
					ExtinctionCoefficientCalculator.PEPTIDE_UNIT_TYPE);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$";
			result = calculator.calculateFromComplexNotation(input);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getDefaultUnit());

			input = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$";
			result = calculator.calculateFromComplexNotation(input,
					ExtinctionCoefficientCalculator.PEPTIDE_UNIT_TYPE);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "ACGTACGT";
			result = calculator.calculateFromNucleotideSequence(input);
			System.out.println("Nucleotide Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getRnaUnit());

			input = "ACGmTACmGT";
			result = calculator.calculateFromModifiedNucleotideSequence(input);
			System.out.println("Modified Nucleotide Sequence\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getRnaUnit());

			input = "P.R(A)P.R(C)P.R(G)P.[mR](A)";
			result = calculator.calculateFromRnaPolymerNotation(input);
			System.out.println("Simple RNA notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getRnaUnit());

			input = "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}$$$$";
			result = calculator.calculateFromComplexNotation(input);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getDefaultUnit());

			input = "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}$$$$";
			result = calculator.calculateFromComplexNotation(input,
					ExtinctionCoefficientCalculator.PEPTIDE_UNIT_TYPE);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getPeptideUnit());

			input = "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}|CHEM1{PEG3}|PEPTIDE1{A.G.G.W.E.E.E.E.E.W}|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$";
			result = calculator.calculateFromComplexNotation(input);
			System.out.println("Complex notation\n" + input);
			System.out.println("Result: " + result + " "
					+ calculator.getDefaultUnit());

		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}
}
