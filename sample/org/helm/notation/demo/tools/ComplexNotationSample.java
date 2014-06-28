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
import org.helm.notation.model.PolymerNode;
import org.helm.notation.model.RNAPolymerNode;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class ComplexNotationSample {

	public static void main(String[] args) {
		try {
			// make sure both monomer db and nucleotide db are in memory
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();

			System.out.println("-----Generic Connection Tests-------");
			testGenericConnection();
			System.out.println();

			System.out.println("-----Canonicalization Tests-------");
			testCanonicalization();
			System.out.println();

			System.out.println("-----Molecule Info Tests-------");
			testMoleculeInfo();
			System.out.println();

			System.out.println("-----Protein Molecule Info Tests-------");
			testProteinMoleculeInfoCalculation();
			System.out.println();

			System.out.println("-----Notation Manipulation Tests-------");
			testNotationManipulation();
			System.out.println();

			System.out.println("-----Comprehensive Tests-------");
			comprehensiveTest();
			System.out.println();

		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void comprehensiveTest() {
		// backbone cyclic peptide
		String notation = "PEPTIDE1{A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$";
		testAll(notation);

		// branch cyclic RNA
		notation = "RNA1{R(C)P.RP.R(A)P.RP.R(A)P.R(U)P}$RNA1,RNA1,4:R3-9:R3$$$";
		testAll(notation);

		// backbone cyclic RNA
		notation = "RNA1{R(C)P.RP.R(A)P.RP.R(A)P.R(U)P}$RNA1,RNA1,1:R1-16:R2$$$";
		testAll(notation);

		// backbone and branch cyclic RNA
		notation = "RNA1{R(C)P.RP.R(A)P.RP.R(A)P.R(U)P}$RNA1,RNA1,4:R3-9:R3|RNA1,RNA1,1:R1-16:R2$$$";
		testAll(notation);

		// cyclic chem
		notation = "CHEM1{SS3}|CHEM2{SS3}$CHEM1,CHEM2,1:R1-1:R1|CHEM1,CHEM2,1:R2-1:R2|$$$";
		testAll(notation);

		// peptide-chem cycles
		notation = "PEPTIDE1{H.H.E.E.E}|CHEM1{SS3}|CHEM2{EG}$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		testAll(notation);

		// multiple peptide-chem cycles
		notation = "PEPTIDE1{E.E.E.E.E}|PEPTIDE2{E.D.D.I.A.C.D.E}|CHEM1{SS3}|CHEM2{SS3}|CHEM3{SS3}$PEPTIDE2,CHEM2,8:R2-1:R1|PEPTIDE1,CHEM3,5:R2-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1|PEPTIDE2,CHEM3,1:R1-1:R1|CHEM1,CHEM2,1:R2-1:R2$$$";
		testAll(notation);

		// conjugate
		notation = "RNA1{R(A)P.R(A)P}|CHEM1{PEG2}$RNA1,CHEM1,6:R2-1:R1$$$";
		testAll(notation);

		// simple RNA
		notation = "RNA1{R(A)P.R(G)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P.R(A)P}$$$$";
		testAll(notation);

		// siRNA with base pair
		notation = "RNA3{R(A)P.R(U)P.R(U)P.R(C)P.R(G)P.R(C)P}|RNA2{R(A)P.R(U)P.R(U)P.R(C)P.R(C)P}|RNA1{R(A)P.R(U)P.R(U)P.R(C)P.R(C)P}$$RNA1,RNA1,2:pair-5:pair|RNA2,RNA2,2:pair-8:pair|RNA3,RNA3,11:pair-14:pair$$";
		testAll(notation);

		// simple CHEM
		notation = "CHEM1{PEG2}$$$$";
		testAll(notation);

		// adhoc CHEM
		notation = "CHEM1{[*]OCCOCCO[*] |$_R1;;;;;;;;_R2$|}$$$$";
		testAll(notation);

		// conjugate with adhoc chem modifier
		notation = "RNA1{R(A)P.R(A)P}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$RNA1,CHEM1,6:R2-1:R1$$$";
		testAll(notation);

		// ADC with dynamic chemical structure and generic edges
		notation = "PEPTIDE1{A.C.A.C.G.K}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1,CHEM1,generic:K-1:R1$$$";
		testAll(notation);

		// ADC with dynamic chemical structure and generic edge
		notation = "CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|PEPTIDE1{A.C.A.C.G.K}$CHEM1,PEPTIDE1,1:R1-C:1.5$$$";
		testAll(notation);

		// ADC with dynamic chemical structure and ambiguous load ratio,
		// nonexistence fuzzy monomer
		notation = "CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|PEPTIDE1{A.E.A.D.G.K}$CHEM1,PEPTIDE1,1:R1-C:1.5$$$";
		testAll(notation);

		// generic connection between peptides
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|PEPTIDE2{A.K.D.C.A}$PEPTIDE1,CHEM1,generic:C-1:R1|PEPTIDE1,PEPTIDE2,generic:K-3:R3$$$";
		testAll(notation);

		// ADC with two fuzzy connection, via one branch, invalid
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1,CHEM1,C:2.0-1:R1|PEPTIDE2,CHEM1,K:0.75-1:R3$$$";
		testAll(notation);

		// ADC with two grouped generic edge, one standard edge
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|CHEM2{PEG2}$PEPTIDE1+PEPTIDE2,CHEM1,generic:Q-1:R1|PEPTIDE1+PEPTIDE2,CHEM2,generic:Q1+Q2-1:R1|CHEM1,CHEM2,1:R3-1:R2$$$";
		testAll(notation);

		// notation =
		// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)}$$$$";

		// notation =
		// "RNA1{R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)}$$$$";

	}

	private static void testAll(String notation) {
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetSmiles(notation);
		testGetCanonicalSmiles(notation);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);
	}

	private static void testCanonicalization() {
		String notation = "PEPTIDE1{H.H.E.E.E}|CHEM1{SS3}|CHEM2{EG}$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		// change node order
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		// change edge order
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM2,5:R2-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		// change edge direction
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM1,1:R1-1:R2|CHEM2,PEPTIDE1,1:R2-5:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		// change node id
		notation = "CHEM4{SS3}|PEPTIDE2{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM4,1:R1-1:R2|CHEM2,PEPTIDE2,1:R2-5:R2|PEPTIDE2,CHEM4,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		// generic connection between peptide group and chem
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K2+K5-1:R1$$$";
		testGetCanonicalNotation(notation);

		// backbone cyclic peptide
		notation = "PEPTIDE1{K.A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$";
		testGetCanonicalNotation(notation);

		// annotated peptides
		notation = "PEPTIDE1{K.A.A.G.K}|PEPTIDE2{K.A.A.G.K}|RNA1{R(A)P.R(G)}|CHEM1{Alexa}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$PEPTIDE1{hc}|PEPTIDE2{lc}$";
		testGetCanonicalNotation(notation);

		// chemical connection
		notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|CHEM1{CovX-2}$RNA1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);

		notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|PEPTIDE1{A.G.G.G.K.K.K.K}|CHEM1{CovX-2}|CHEM2{3Bio}$PEPTIDE1,CHEM2,8:R2-1:R1|RNA1,CHEM1,1:R1-1:R1$$$";
		testGetCanonicalNotation(notation);
	}

	private static void testGenericConnection() {
		// generic connection between peptide and chem
		String notation = "PEPTIDE1{A.C.A.C.G.K}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1,CHEM1,generic:K-1:R1$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// generic connection between peptides: one generic, one standard
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}$PEPTIDE1,PEPTIDE2,generic:K-3:R3$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// generic connection between peptides: two generics
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}$PEPTIDE1,PEPTIDE2,generic:K-generic:E$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// generic connection between peptide group and chem
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K2+K5-1:R1$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// generic connections between one peptide group and two chem modifiers
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|CHEM2{PEG2}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K-1:R1|PEPTIDE1+PEPTIDE2,CHEM2,generic:Q1+Q2-1:R1$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// generic connection between two peptide groups
		notation = "PEPTIDE1{A.C.A.K.C.G.K.E.E}|PEPTIDE2{G.A.C.A.C.G.K.E.E}|PEPTIDE3{A.C.A.C.G.K.E.E}|PEPTIDE4{A.C.A.E.D.C.G.K.E.E}$PEPTIDE1+PEPTIDE4,PEPTIDE3+PEPTIDE2,generic:K-generic:D$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// two unknown connections
		notation = "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}|CHEM1{C1=CC=CC=C1 |c:0,2,4|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:?$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

		// one unknown connection
		notation = "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}|CHEM1{[*]C1=CC=CC=C1 |$_R8;;;;;;$,c:3,5,t:1|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:R8$$$";
		testComplexNotationValidity(notation);
		testGetCanonicalNotation(notation);
		testGetMoleculeInfo(notation);

	}

	private static void testMoleculeInfo() {
		// siRNA
		String notation = "RNA1{R(A)P.R(U)P.R(C)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)}|RNA2{R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.[dR](A)P.[dR](T)}$$RNA1,RNA2,2:pair-74:pair|RNA1,RNA2,5:pair-71:pair|RNA1,RNA2,8:pair-68:pair|RNA1,RNA2,11:pair-65:pair|RNA1,RNA2,14:pair-62:pair|RNA1,RNA2,17:pair-59:pair|RNA1,RNA2,20:pair-56:pair|RNA1,RNA2,23:pair-53:pair|RNA1,RNA2,26:pair-50:pair|RNA1,RNA2,29:pair-47:pair|RNA1,RNA2,32:pair-44:pair|RNA1,RNA2,35:pair-41:pair|RNA1,RNA2,38:pair-38:pair|RNA1,RNA2,41:pair-35:pair|RNA1,RNA2,44:pair-32:pair|RNA1,RNA2,47:pair-29:pair|RNA1,RNA2,50:pair-26:pair|RNA1,RNA2,53:pair-23:pair|RNA1,RNA2,56:pair-20:pair|RNA1,RNA2,59:pair-17:pair|RNA1,RNA2,62:pair-14:pair|RNA1,RNA2,65:pair-11:pair|RNA1,RNA2,68:pair-8:pair|RNA1,RNA2,71:pair-5:pair|RNA1,RNA2,74:pair-2:pair$$";
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// conjugate
		notation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);
	}

	private static void testProteinMoleculeInfoCalculation() {
		// 100 monomers
		String notation = createPeptideNotation(100);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 150 monomers
		notation = createPeptideNotation(150);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 199 monomers
		notation = createPeptideNotation(199);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 400 monomers
		notation = createPeptideNotation(400);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 800 monomers
		notation = createPeptideNotation(800);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 1600 monomers
		notation = createPeptideNotation(1600);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);

		// 3200 monomers
		notation = createPeptideNotation(3200);
		testGetMoleculeInfo(notation);
		testGetMoleculeInfoViaSmiles(notation);
	}

	private static void testNotationManipulation() {
		try {
			String notation = "RNA1{R(G)P.[fR](U)P.R(A)P.R(A)P.R(G)P.R(A)P.[fR](C)P.[fR](U)P.[fR](A)P.R(C)P.R(U)P.R(C)P.R(A)P.[fR](U)P.R(G)P.R(A)P.[fR](U)P.[fR](C)P.[fR](C)P.R(T)[sP].R(T)}|RNA2{R(G)P.R(G)P.R(A)P.[fR](U)P.[fR](C)P.R(A)P.[fR](U)P.[fR](G)P.[fR](A)P.[fR](G)P.R(U)P.R(A)P.R(G)P.[fR](U)P.[fR](C)P.[fR](U)P.[fR](U)P.R(A)P.[fR](C)P.R(T)[sP].R(T)}$$RNA2,RNA1,29:pair-29:pair|RNA2,RNA1,17:pair-41:pair|RNA2,RNA1,38:pair-20:pair|RNA2,RNA1,47:pair-11:pair|RNA2,RNA1,26:pair-32:pair|RNA2,RNA1,11:pair-47:pair|RNA2,RNA1,56:pair-2:pair|RNA2,RNA1,2:pair-56:pair|RNA2,RNA1,8:pair-50:pair|RNA2,RNA1,41:pair-17:pair|RNA2,RNA1,44:pair-14:pair|RNA2,RNA1,35:pair-23:pair|RNA2,RNA1,5:pair-53:pair|RNA2,RNA1,14:pair-44:pair|RNA2,RNA1,32:pair-26:pair|RNA2,RNA1,50:pair-8:pair|RNA2,RNA1,20:pair-38:pair|RNA2,RNA1,53:pair-5:pair|RNA2,RNA1,23:pair-35:pair$RNA1{as}|RNA2{ss}$";
			System.out.println("getRNAPolymerNodeList start: "
					+ System.currentTimeMillis());
			List<RNAPolymerNode> l = ComplexNotationParser
					.getRNAPolymerNodeList(notation);
			System.out.println("getRNAPolymerNodeList end: "
					+ System.currentTimeMillis());
			for (int i = 0; i < l.size(); i++) {
				RNAPolymerNode node = l.get(i);
				System.out.println("Strand Annotation: " + node.getAnotation());
				System.out.println("Sequence: " + node.getSequence());
				System.out.println("Modified Sequence: "
						+ node.getModifiedSequence());
			}

			notation = "RNA1{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}$$RNA2,RNA1,20:pair-2:pair|RNA2,RNA1,5:pair-17:pair|RNA2,RNA1,2:pair-20:pair|RNA2,RNA1,11:pair-11:pair|RNA2,RNA1,17:pair-5:pair|RNA2,RNA1,14:pair-8:pair|RNA2,RNA1,8:pair-14:pair$RNA1{ss}$";
			// notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)}$$$$";
			// notation =
			// "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}|RNA2{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}$$RNA1,RNA2,5:pair-17:pair|RNA1,RNA2,11:pair-11:pair|RNA1,RNA2,20:pair-2:pair|RNA1,RNA2,2:pair-20:pair|RNA1,RNA2,14:pair-8:pair$$";
			String[] formatedSeqs = ComplexNotationParser
					.getFormatedSirnaSequences(notation, "*", "|");
			for (int i = 0; i < formatedSeqs.length; i++) {
				System.out.println(formatedSeqs[i]);
			}

			formatedSeqs = ComplexNotationParser
					.getFormatedSirnaSequences(notation);
			for (int i = 0; i < formatedSeqs.length; i++) {
				System.out.println(formatedSeqs[i]);
			}

			String notation1 = "RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(G)P.R(G)}|CHEM1{PEG3}$RNA3,CHEM1,4:R2-1:R2|RNA1,CHEM1,10:R2-1:R1$$$";
			String notation2 = "CHEM1{PEG3}|RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(G)P.R(G)}$CHEM1,RNA1,1:R1-10:R2|CHEM1,RNA3,1:R2-4:R2$$$";
			String combinedNotation = ComplexNotationParser
					.getCombinedComlexNotation(notation1, notation2);
			System.out.println(combinedNotation);

			System.out.println("Hybridization Sample:");
			notation = "RNA1{P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}|RNA2{P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$";
			String hybridizedNotation = ComplexNotationParser
					.hybridize(notation);
			System.out.println(hybridizedNotation);

			System.out.println("Decomposition Sample:");
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|CHEM1{A6P}$RNA1,CHEM1,1:R1-1:R2$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,8:pair-20:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,20:pair-8:pair$RNA3{ss}$";
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}$RNA1,RNA1,21:R2-1:R1$$$"; // cycle
																																					// +
																																					// linear
																																					// mixture
			notation = "RNA1{R(A)}|RNA2{R(A)}|CHEM1{CovX-3}|CHEM2{CovX-2}$RNA1,CHEM2,1:R1-1:R1|RNA1,CHEM1,1:R2-1:R1$$$"; // multi
																															// connections
			notation = "PEPTIDE1{A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE2{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE3{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE4{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C}$PEPTIDE2,PEPTIDE4,18:R3-18:R3$$PEPTIDE1{hc}$";
			String[] notations = ComplexNotationParser.decompose(notation);
			for (int i = 0; i < notations.length; i++) {
				System.out.println(notations[i]);
			}

			System.out
					.println("Parse complex notation to sense and antisense notations:");
			String ssNotation = null;
			String asNotation = null;
			boolean duplex = false;
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$";
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA1{ss}$";
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}$$$RNA1{as}$";
			String allNodeNotation = ComplexNotationParser
					.getAllNodeString(notation);
			List<PolymerNode> pnList = ComplexNotationParser
					.getPolymerNodeList(allNodeNotation);
			if (pnList.size() == 0 || pnList.size() > 2) {
				System.out.println("more than two strands");
			} else {
				List<RNAPolymerNode> rnaList = ComplexNotationParser
						.getRNAPolymerNodeList(notation);
				if (rnaList.size() == pnList.size()) {
					for (RNAPolymerNode rna : rnaList) {
						String annotation = rna.getAnotation();
						if (null == annotation
								|| annotation.trim().length() == 0) {
							System.out
									.println("tell users about missing annotation");
						} else {
							if (annotation
									.equalsIgnoreCase(SimpleNotationParser.RNA_ANTISENSE_STRAND_ANNOTATION)) {
								asNotation = rna.getLabel();
							} else if (annotation
									.equalsIgnoreCase(SimpleNotationParser.RNA_SENSE_STRAND_ANNOTATION)) {
								ssNotation = rna.getLabel();
							} else {
								System.out.println("notify user");
							}
						}
					}
				} else {
					System.out.println("there are none RNA polymers");
				}
			}

			if (null != ssNotation && null != asNotation) {
				duplex = true;
			}

			System.out.println("SS:" + ssNotation);
			System.out.println("AS:" + asNotation);
			System.out.println("Duplex:" + duplex);

			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$";
			System.out.println("Replace P with sP for: " + notation);
			String result = ComplexNotationParser.replaceMonomer(notation,
					Monomer.NUCLIEC_ACID_POLYMER_TYPE, "P", "sP");
			System.out.println("Result: " + result);

			String[] rtcNotations = { "R(A)P.R(G)", "A.A.G.L", "PEG2",
					"RNA1{P.R(A)P.R(G)}$$$$" };
			for (String rtcNotation : rtcNotations) {
				String cn = ComplexNotationParser.standardize(rtcNotation);
				System.out.println(rtcNotation + " => " + cn);
			}

			System.out.println("Hybridization bug fix:");
			notation = "RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)[sP].R(U)}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}$$$RNA1{ss}|RNA2{as}$";
			hybridizedNotation = ComplexNotationParser.hybridize(notation);
			System.out.println(hybridizedNotation);

			System.out.println("Test for modification:");
			notation = "RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)[sP].R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)[sP].R(U)}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}$$$RNA1{ss}|RNA2{as}$";
			List<PolymerNode> polyNodeList = ComplexNotationParser
					.getPolymerNodeList(notation);
			boolean hasMod = ComplexNotationParser
					.hasNucleotideModification(polyNodeList);
			System.out.println("hasMod? of " + notation + " = " + hasMod);

			System.out.println("Test for modification:");
			notation = "RNA1{R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)P.R(U)P}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.R(C)P}$$$RNA1{ss}|RNA2{as}$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser
					.hasNucleotideModification(polyNodeList);
			System.out.println("hasMod? of " + notation + " = " + hasMod);

			System.out.println("Test for modification:");
			notation = "RNA1{R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(C)}|RNA2{R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)}$$RNA1,RNA2,2:pair-56:pair|RNA1,RNA2,5:pair-53:pair|RNA1,RNA2,8:pair-50:pair|RNA1,RNA2,11:pair-47:pair|RNA1,RNA2,14:pair-44:pair|RNA1,RNA2,17:pair-41:pair|RNA1,RNA2,20:pair-38:pair|RNA1,RNA2,23:pair-35:pair|RNA1,RNA2,26:pair-32:pair|RNA1,RNA2,29:pair-29:pair|RNA1,RNA2,32:pair-26:pair|RNA1,RNA2,35:pair-23:pair|RNA1,RNA2,38:pair-20:pair|RNA1,RNA2,41:pair-17:pair|RNA1,RNA2,44:pair-14:pair|RNA1,RNA2,47:pair-11:pair|RNA1,RNA2,50:pair-8:pair|RNA1,RNA2,53:pair-5:pair|RNA1,RNA2,56:pair-2:pair$RNA1{ss}|RNA2{as}$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser
					.hasNucleotideModification(polyNodeList);
			System.out.println("hasMod? of " + notation + " = " + hasMod);

			System.out.println("Test for modification:");
			notation = "RNA1{R(A)P.R(U)P.R(C)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)}|RNA2{R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.[dR](A)P.[dR](T)}$$RNA1,RNA2,2:pair-74:pair|RNA1,RNA2,5:pair-71:pair|RNA1,RNA2,8:pair-68:pair|RNA1,RNA2,11:pair-65:pair|RNA1,RNA2,14:pair-62:pair|RNA1,RNA2,17:pair-59:pair|RNA1,RNA2,20:pair-56:pair|RNA1,RNA2,23:pair-53:pair|RNA1,RNA2,26:pair-50:pair|RNA1,RNA2,29:pair-47:pair|RNA1,RNA2,32:pair-44:pair|RNA1,RNA2,35:pair-41:pair|RNA1,RNA2,38:pair-38:pair|RNA1,RNA2,41:pair-35:pair|RNA1,RNA2,44:pair-32:pair|RNA1,RNA2,47:pair-29:pair|RNA1,RNA2,50:pair-26:pair|RNA1,RNA2,53:pair-23:pair|RNA1,RNA2,56:pair-20:pair|RNA1,RNA2,59:pair-17:pair|RNA1,RNA2,62:pair-14:pair|RNA1,RNA2,65:pair-11:pair|RNA1,RNA2,68:pair-8:pair|RNA1,RNA2,71:pair-5:pair|RNA1,RNA2,74:pair-2:pair$$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser
					.hasNucleotideModification(polyNodeList);
			System.out.println("hasMod? of " + notation + " = " + hasMod);
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetPolymerNodeList(String notation) {
		try {
			System.out.println("Testing getPolymerNodeList for: " + notation);
			System.out.println("getPolymerNodeList start: "
					+ System.currentTimeMillis());
			List<PolymerNode> pnl = ComplexNotationParser
					.getPolymerNodeList(notation);
			System.out.println("getPolymerNodeList end: "
					+ System.currentTimeMillis());
			for (int i = 0; i < pnl.size(); i++) {
				PolymerNode node = pnl.get(i);
				System.out.println("ID: " + node.getId());
				System.out.println("Label: " + node.getLabel());
			}
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testComplexNotationValidity(String notation) {
		try {
			System.out.println("Testing Notation Validity for: " + notation);
			System.out.println("Notation Format Valid? "
					+ ComplexNotationParser.validateNotationFormat(notation));
			System.out.println("Node String: "
					+ ComplexNotationParser.getAllNodeString(notation));
			System.out.println("Edge String: "
					+ ComplexNotationParser.getAllEdgeString(notation));
			System.out.println("Base Pair String: "
					+ ComplexNotationParser.getAllBasePairString(notation));
			System.out.println("Node Label String: "
					+ ComplexNotationParser.getAllNodeLabelString(notation));
			System.out.println("Others String: "
					+ ComplexNotationParser.getOtherString(notation));
			System.out.println("Notation Valid? "
					+ ComplexNotationParser.validateComplexNotation(notation));
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetCanonicalNotation(String notation) {
		try {
			System.out.println("Testing getCanonicalNotation for: " + notation);
			System.out.println("Canonical Notation: "
					+ ComplexNotationParser.getCanonicalNotation(notation));
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetSmiles(String notation) {
		try {
			System.out.println("Testing getSmiles for: " + notation);
			System.out.println("SMILES: "
					+ ComplexNotationParser.getComplexPolymerSMILES(notation));
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetCanonicalSmiles(String notation) {
		try {
			System.out.println("Testing getCanonicalSmiles for: " + notation);
			System.out.println("Canonical SMILES: "
					+ ComplexNotationParser
							.getComplexPolymerCanonicalSmiles(notation));
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetMoleculeInfo(String notation) {
		System.out.println("testGetMoleculeInfo...");
		try {
			long start = System.currentTimeMillis();
			System.out.println("Start Time: " + start);
			MoleculeInfo mi = ComplexNotationParser.getMoleculeInfo(notation);
			long end = System.currentTimeMillis();
			System.out.println("End Time: " + end);
			System.out.println("MW = " + mi.getMolecularWeight());
			System.out.println("MF = " + mi.getMolecularFormula());
			System.out.println("Mass = " + mi.getExactMass());
			System.out.println("Time = " + (end - start));
			System.out.println("*********************");
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static void testGetMoleculeInfoViaSmiles(String notation) {
		System.out.println("testGetMoleculeInfoViaSmiles...");
		try {
			long start = System.currentTimeMillis();
			System.out.println("Start Time: " + start);
			String smiles = ComplexNotationParser
					.getComplexPolymerSMILES(notation);
			MoleculeInfo mi = StructureParser.getMoleculeInfo(smiles);
			long end = System.currentTimeMillis();
			System.out.println("End Time: " + end);
			System.out.println("MW = " + mi.getMolecularWeight());
			System.out.println("MF = " + mi.getMolecularFormula());
			System.out.println("Mass = " + mi.getExactMass());
			System.out.println("Time = " + (end - start));
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	private static String createPeptideNotation(int monomerCount) {
		StringBuilder sb = new StringBuilder();
		sb.append("PEPTIDE1{");
		for (int i = 0; i < monomerCount; i++) {
			sb.append("A");
			if (i < monomerCount - 1) {
				sb.append(".");
			}
		}
		sb.append("}$$$$");
		return sb.toString();
	}
}
