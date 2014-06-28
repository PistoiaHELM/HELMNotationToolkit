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
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * 
 * @author zhangtianhong
 */
public class SimpleNotationSample {

	public static void main(String[] args) {
		String notation = null;
		try {
			// initialize monomer and nucleotide database
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();

			notation = "P.R(A)[sP].RP.R(G)P.[LR]([5meC])";
			testGetComplexNotation(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE);
			testGetCanonicalNotation(notation,
					Monomer.NUCLIEC_ACID_POLYMER_TYPE);
			testValidation(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE);
			testGetMonomerCount(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE);
			testGetMoleculeInfo(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE);
			testReplaceMonomer(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE,
					"P", "sP");
			testGetSMILES(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE);

			notation = "G.G.K.A.A.[seC]";
			testGetComplexNotation(notation, Monomer.PEPTIDE_POLYMER_TYPE);
			testGetCanonicalNotation(notation, Monomer.PEPTIDE_POLYMER_TYPE);
			testValidation(notation, Monomer.PEPTIDE_POLYMER_TYPE);
			testGetMonomerCount(notation, Monomer.PEPTIDE_POLYMER_TYPE);
			testGetMoleculeInfo(notation, Monomer.PEPTIDE_POLYMER_TYPE);
			testReplaceMonomer(notation, Monomer.PEPTIDE_POLYMER_TYPE, "A", "Q");
			testGetSMILES(notation, Monomer.PEPTIDE_POLYMER_TYPE);

			notation = "PEG2";
			testGetComplexNotation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetCanonicalNotation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testValidation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetMonomerCount(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetMoleculeInfo(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testReplaceMonomer(notation, Monomer.CHEMICAL_POLYMER_TYPE, "PEG2",
					"SS3");
			testGetSMILES(notation, Monomer.CHEMICAL_POLYMER_TYPE);

			notation = "[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|";
			testGetComplexNotation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetCanonicalNotation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testValidation(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetMonomerCount(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetMoleculeInfo(notation, Monomer.CHEMICAL_POLYMER_TYPE);
			testGetSMILES(notation, Monomer.CHEMICAL_POLYMER_TYPE);

			notation = "R(C)P.R(G)P.R(A)P.R(U)P.R(A)P.R(U)P.R(G)P.R(G)P.R(G)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.[dR](U)P.[dR](U)";
			System.out.println("getNucleotideList Start: "
					+ System.currentTimeMillis());
			List<Nucleotide> nucleotideList = SimpleNotationParser
					.getNucleotideList(notation);
			System.out.println("getNucleotideList End: "
					+ System.currentTimeMillis());

			System.out.println("getNucleotideList wo validation Start : "
					+ System.currentTimeMillis());
			nucleotideList = SimpleNotationParser.getNucleotideList(notation,
					false);
			System.out.println("getNucleotideList wo validation End: "
					+ System.currentTimeMillis());

			for (int i = 0; i < nucleotideList.size(); i++) {
				Nucleotide nuc = nucleotideList.get(i);
				System.out.println("Symbol: " + nuc.getSymbol());
				System.out.println("Modified: " + nuc.isModified());
				System.out.println("Notation: " + nuc.getNotation());
			}

			System.out.println("Sequence: "
					+ SimpleNotationParser.getNucleotideSequence(notation));

			notation = "K.C.C.C.W.K.[seC]";
			System.out.println("Peptide Sequence: "
					+ SimpleNotationParser.getPeptideSequence(notation));
			System.out.println("Modified Peptide Sequence: "
					+ SimpleNotationParser.getModifiedPeptideSequence(notation,
							"|"));

			notation = "R(A)[sP].R(G)P.R(C).PEG.[LR]([5meC])[sP].R(PEG)";
			// valid = SimpleNotationParser.getStrictNucleotideList(notation,
			// true);
			System.out.println("Input sequence " + notation);
			nucleotideList = SimpleNotationParser.getStrictNucleotideList(
					notation, false);
			StringBuffer outputList = new StringBuffer();
			for (int i = 0; i < nucleotideList.size(); i++) {
				Nucleotide nuc = nucleotideList.get(i);
				System.out.println("Symbol: " + nuc.getSymbol());
				System.out.println("Modified: " + nuc.isModified());
				System.out.println("Notation: " + nuc.getNotation());
				outputList.append(nuc.getSymbol());
			}
			System.out.println("Output List: " + outputList.toString());

		} catch (Exception e) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, e);
		}
	}

	private static void testReplaceMonomer(String notation, String polymerType,
			String oldMonomerID, String newMonomerID) {
		try {
			System.out.println("replacing monomer " + oldMonomerID + " with "
					+ newMonomerID + " in " + notation + " as " + polymerType);

			String result = SimpleNotationParser.replaceMonomer(notation,
					polymerType, oldMonomerID, newMonomerID);
			System.out.println("Result: " + result);
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetMonomerCount(String notation, String polymerType) {
		try {
			System.out.println("Getting monomer count for " + notation + " as "
					+ polymerType);
			int count = SimpleNotationParser.getMonomerCount(notation,
					polymerType);
			System.out.println(count);
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetComplexNotation(String notation,
			String polymerType) {
		try {
			System.out.println("Getting complex notation for " + notation
					+ " as " + polymerType);
			String complexNotation = SimpleNotationParser.getComplexNotation(
					notation, polymerType);
			System.out.println(complexNotation);
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testValidation(String notation, String polymerType) {
		try {
			System.out.println("Validating " + notation + " as " + polymerType);
			boolean valid = SimpleNotationParser.validateSimpleNotation(
					notation, polymerType);
			System.out.println(valid);
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetMoleculeInfo(String notation, String polymerType) {
		try {
			System.out.println("Getting molecule info for " + notation + " as "
					+ polymerType);
			MoleculeInfo mi = SimpleNotationParser.getMoleculeInfo(notation,
					polymerType);
			System.out.println("MW: " + mi.getMolecularWeight());
			System.out.println("MF: " + mi.getMolecularFormula());
			System.out.println("Mass: " + mi.getExactMass());
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetCanonicalNotation(String notation,
			String polymerType) {
		try {
			System.out.println("Getting canonical notation for " + notation
					+ " as " + polymerType);
			String canNotation = SimpleNotationParser
					.getSimpleCanonicalNotation(notation, polymerType);
			System.out.println(canNotation);
			Map.Entry entry = SimpleNotationParser
					.getSimpleCanonicalNotationMapEntry(notation, polymerType);
			System.out.println("Position: " + entry.getKey().toString()
					+ "  Notation: " + entry.getValue());
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetSMILES(String notation, String polymerType) {
		try {
			System.out.println("Getting SMILES for " + notation + " as "
					+ polymerType);
			String smiles = SimpleNotationParser.getSimplePolymerSMILES(
					notation, polymerType);
			System.out.println(smiles);
			System.out.println("***************");
		} catch (Exception ex) {
			Logger.getLogger(SimpleNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}
}
