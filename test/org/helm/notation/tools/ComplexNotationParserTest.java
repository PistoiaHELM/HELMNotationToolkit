/*******************************************************************************
 * Copyright C 2014, The Pistoia Alliance
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

import org.helm.notation.tools.*;
import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.StructureException;
import org.helm.notation.demo.tools.ComplexNotationSample;
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;
import org.helm.notation.model.RNAPolymerNode;

import java.io.IOException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.*;

import org.jdom.JDOMException;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ComplexNotationParserTest {


	@Before
	public void init() {
		try {
			MonomerFactory.finalizeMonomerCache();
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();	
			SimpleNotationParser.resetSeed();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@After
	public void finish() {
		MonomerFactory.finalizeMonomerCache();
	}
	
	@Test
	public void testGenericConnection() {
		//generic connection between peptide and chem
		String notation = "PEPTIDE1{A.C.A.C.G.K}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1,CHEM1,generic:K-1:R1$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//generic connection between peptides: one generic, one standard
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}$PEPTIDE1,PEPTIDE2,generic:K-3:R3$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//generic connection between peptides: two generics
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}$PEPTIDE1,PEPTIDE2,generic:K-generic:E$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//generic connection between peptide group and chem
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K2+K5-1:R1$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//generic connections between one peptide group and two chem modifiers
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|CHEM2{PEG2}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K-1:R1|PEPTIDE1+PEPTIDE2,CHEM2,generic:Q1+Q2-1:R1$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//generic connection between two peptide groups
		notation = "PEPTIDE1{A.C.A.K.C.G.K.E.E}|PEPTIDE2{G.A.C.A.C.G.K.E.E}|PEPTIDE3{A.C.A.C.G.K.E.E}|PEPTIDE4{A.C.A.E.D.C.G.K.E.E}$PEPTIDE1+PEPTIDE4,PEPTIDE3+PEPTIDE2,generic:K-generic:D$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//two unknown connections
		notation = "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}|CHEM1{C1=CC=CC=C1 |c:0,2,4|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:?$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

		//one unknown connection
		notation = "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}|CHEM1{[*]C1=CC=CC=C1 |$_R8;;;;;;$,c:3,5,t:1|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:R8$$$";
		assertTrue(testComplexNotationValidity(notation));
		assertTrue(testGetCanonicalNotation(notation));
		assertTrue(testGetMoleculeInfo(notation));

	}


	@Test
	public void testCanonicalization() {

		String notation = "PEPTIDE1{H.H.E.E.E}|CHEM1{SS3}|CHEM2{EG}$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));

		//change node order
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));                     

		//change edge order
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM2,5:R2-1:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));

		//change edge direction
		notation = "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM1,1:R1-1:R2|CHEM2,PEPTIDE1,1:R2-5:R2|PEPTIDE1,CHEM1,1:R1-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));

		//change node id
		notation = "CHEM4{SS3}|PEPTIDE2{H.H.E.E.E}|CHEM2{EG}$CHEM2,CHEM4,1:R1-1:R2|CHEM2,PEPTIDE2,1:R2-5:R2|PEPTIDE2,CHEM4,1:R1-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));  

		//generic connection between peptide group and chem
		notation = "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1+PEPTIDE2,CHEM1,generic:K2+K5-1:R1$$$";
		assertTrue(testGetCanonicalNotation(notation));

		//backbone cyclic peptide
		notation = "PEPTIDE1{K.A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$";
		assertTrue(testGetCanonicalNotation(notation));

		//annotated peptides
		notation = "PEPTIDE1{K.A.A.G.K}|PEPTIDE2{K.A.A.G.K}|RNA1{R(A)P.R(G)}|CHEM1{Alexa}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$PEPTIDE1{hc}|PEPTIDE2{lc}$";
		assertTrue(testGetCanonicalNotation(notation));

		//chemical connection; CovX is not registered
		notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|CHEM1{CovX-2}$RNA1,CHEM1,1:R1-1:R1$$$";
		assertFalse(testGetCanonicalNotation(notation));

		notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|PEPTIDE1{A.G.G.G.K.K.K.K}|CHEM1{CovX-2}|CHEM2{3Bio}$PEPTIDE1,CHEM2,8:R2-1:R1|RNA1,CHEM1,1:R1-1:R1$$$";
		assertFalse(testGetCanonicalNotation(notation));
	}


	@Test
	public void testMoleculeInfo() {

		//siRNA
		String notation = "RNA1{R(A)P.R(U)P.R(C)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)}|RNA2{R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.[dR](A)P.[dR](T)}$$RNA1,RNA2,2:pair-74:pair|RNA1,RNA2,5:pair-71:pair|RNA1,RNA2,8:pair-68:pair|RNA1,RNA2,11:pair-65:pair|RNA1,RNA2,14:pair-62:pair|RNA1,RNA2,17:pair-59:pair|RNA1,RNA2,20:pair-56:pair|RNA1,RNA2,23:pair-53:pair|RNA1,RNA2,26:pair-50:pair|RNA1,RNA2,29:pair-47:pair|RNA1,RNA2,32:pair-44:pair|RNA1,RNA2,35:pair-41:pair|RNA1,RNA2,38:pair-38:pair|RNA1,RNA2,41:pair-35:pair|RNA1,RNA2,44:pair-32:pair|RNA1,RNA2,47:pair-29:pair|RNA1,RNA2,50:pair-26:pair|RNA1,RNA2,53:pair-23:pair|RNA1,RNA2,56:pair-20:pair|RNA1,RNA2,59:pair-17:pair|RNA1,RNA2,62:pair-14:pair|RNA1,RNA2,65:pair-11:pair|RNA1,RNA2,68:pair-8:pair|RNA1,RNA2,71:pair-5:pair|RNA1,RNA2,74:pair-2:pair$$";
		assertTrue(testGetMoleculeInfo(notation));
		assertTrue(testGetMoleculeInfoViaSmiles(notation));

		//conjugate
		notation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		assertTrue(testGetMoleculeInfo(notation));
		assertTrue(testGetMoleculeInfoViaSmiles(notation));
	}


	@Test
	public void testProteinMoleculeInfoCalculation() {

		//100 monomers
		String notation = createPeptideNotation(100);
		System.out.println("100");
		assertTrue(testGetMoleculeInfo(notation));
		assertTrue(testGetMoleculeInfoViaSmiles(notation));

		//150 monomers
		System.out.println("150");
		notation = createPeptideNotation(150);
		assertTrue(testGetMoleculeInfo(notation));
		assertTrue(testGetMoleculeInfoViaSmiles(notation));

		//199 monomers
		System.out.println("199");
		notation = createPeptideNotation(199);
		assertTrue(testGetMoleculeInfo(notation));
		assertTrue(testGetMoleculeInfoViaSmiles(notation));

		//400 monomers
		System.out.println("400");
		notation = createPeptideNotation(400);
		assertTrue(testGetMoleculeInfo(notation));
		assertFalse(testGetMoleculeInfoViaSmiles(notation));

		//800 monomers
		System.out.println("800");
		notation = createPeptideNotation(800);
		assertTrue(testGetMoleculeInfo(notation));
		assertFalse(testGetMoleculeInfoViaSmiles(notation));

		//1600 monomers
		System.out.println("1600");
		notation = createPeptideNotation(1600);
		assertTrue(testGetMoleculeInfo(notation));
		assertFalse(testGetMoleculeInfoViaSmiles(notation));

		//3200 monomers
		System.out.println("3200");
		notation = createPeptideNotation(3200);
		assertTrue(testGetMoleculeInfo(notation));
		assertFalse(testGetMoleculeInfoViaSmiles(notation));
	}


	@Test
	public void testNotationManipulation() {
		try {
			String notation = "RNA1{R(G)P.[fR](U)P.R(A)P.R(A)P.R(G)P.R(A)P.[fR](C)P.[fR](U)P.[fR](A)P.R(C)P.R(U)P.R(C)P.R(A)P.[fR](U)P.R(G)P.R(A)P.[fR](U)P.[fR](C)P.[fR](C)P.R(T)[sP].R(T)}|RNA2{R(G)P.R(G)P.R(A)P.[fR](U)P.[fR](C)P.R(A)P.[fR](U)P.[fR](G)P.[fR](A)P.[fR](G)P.R(U)P.R(A)P.R(G)P.[fR](U)P.[fR](C)P.[fR](U)P.[fR](U)P.R(A)P.[fR](C)P.R(T)[sP].R(T)}$$RNA2,RNA1,29:pair-29:pair|RNA2,RNA1,17:pair-41:pair|RNA2,RNA1,38:pair-20:pair|RNA2,RNA1,47:pair-11:pair|RNA2,RNA1,26:pair-32:pair|RNA2,RNA1,11:pair-47:pair|RNA2,RNA1,56:pair-2:pair|RNA2,RNA1,2:pair-56:pair|RNA2,RNA1,8:pair-50:pair|RNA2,RNA1,41:pair-17:pair|RNA2,RNA1,44:pair-14:pair|RNA2,RNA1,35:pair-23:pair|RNA2,RNA1,5:pair-53:pair|RNA2,RNA1,14:pair-44:pair|RNA2,RNA1,32:pair-26:pair|RNA2,RNA1,50:pair-8:pair|RNA2,RNA1,20:pair-38:pair|RNA2,RNA1,53:pair-5:pair|RNA2,RNA1,23:pair-35:pair$RNA1{as}|RNA2{ss}$";
			//System.out.println("getRNAPolymerNodeList start: " + System.currentTimeMillis());
			List<RNAPolymerNode> l = ComplexNotationParser.getRNAPolymerNodeList(notation);

			assertEquals("RNAPolymerNodeCount", l.size(),2);
			RNAPolymerNode node1 = l.get(0);
			assertEquals("Strand Annotation 1", node1.getAnotation(), "as");
			assertEquals("Sequence1", node1.getSequence(), "GUAAGACUACUCAUGAUCCTT");
			assertEquals("Modified Sequence1", node1.getModifiedSequence(), "GfUAAGAfCfUfACUCAfUGAfUfCfCTT");

			RNAPolymerNode node2 = l.get(1);
			assertEquals("Strand Annotation 2", node2.getAnotation(), "ss");
			assertEquals("Sequence1", node2.getSequence(), "GGAUCAUGAGUAGUCUUACTT");
			assertEquals("Modified Sequence1", node2.getModifiedSequence(), "GGAfUfCAfUfGfAfGUAGfUfCfUfUAfCTT");

			
			notation = "RNA1{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}$$RNA2,RNA1,20:pair-2:pair|RNA2,RNA1,5:pair-17:pair|RNA2,RNA1,2:pair-20:pair|RNA2,RNA1,11:pair-11:pair|RNA2,RNA1,17:pair-5:pair|RNA2,RNA1,14:pair-8:pair|RNA2,RNA1,8:pair-14:pair$RNA1{ss}$";
			String[] formatedSeqs = ComplexNotationParser.getFormatedSirnaSequences(notation, "*", "|");
			assertEquals("FormatedSirnaSequence 0", formatedSeqs[0],  "**UUAAGCUTT");
			assertEquals("FormatedSirnaSequence 1", formatedSeqs[1],  "**|||||||**");
			assertEquals("FormatedSirnaSequence 2", formatedSeqs[2],  "TTAAUUCGA**");
			

			formatedSeqs = ComplexNotationParser.getFormatedSirnaSequences(notation);
			assertEquals("FormatedSirnaSequence 0", formatedSeqs[0],  "  UUAAGCUTT");
			assertEquals("FormatedSirnaSequence 1", formatedSeqs[1],  "  |||||||  ");
			assertEquals("FormatedSirnaSequence 2", formatedSeqs[2],  "TTAAUUCGA  ");
			

			String notation1 = "RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(G)P.R(G)}|CHEM1{PEG3}$RNA3,CHEM1,4:R2-1:R2|RNA1,CHEM1,10:R2-1:R1$$$";
			String notation2 = "CHEM1{PEG3}|RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(G)P.R(G)}$CHEM1,RNA1,1:R1-10:R2|CHEM1,RNA3,1:R2-4:R2$$$";
			String combinedNotation = ComplexNotationParser.getCombinedComlexNotation(notation1, notation2);
			assertEquals("CombinedComLexNotation", combinedNotation, "RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(G)P.R(G)}|CHEM1{PEG3}|CHEM2{PEG3}|RNA4{R(A)P.R(G)P.R(C)P.R(U)}|RNA5{R(A)P.R(G)P.R(C)P.R(U)}|RNA6{R(G)P.R(G)}$RNA3,CHEM1,4:R2-1:R2|RNA1,CHEM1,10:R2-1:R1|CHEM2,RNA4,1:R1-10:R2|CHEM2,RNA6,1:R2-4:R2$$$");
			

			notation = "RNA1{P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}|RNA2{P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$";
			String hybridizedNotation = ComplexNotationParser.hybridize(notation);
			assertEquals("Hybridize", hybridizedNotation, "RNA1{P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}|RNA2{P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,3:pair-18:pair|RNA1,RNA2,6:pair-15:pair|RNA1,RNA2,9:pair-12:pair|RNA1,RNA2,12:pair-9:pair|RNA1,RNA2,15:pair-6:pair|RNA1,RNA2,18:pair-3:pair$$");

			
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}|RNA3{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|CHEM1{A6P}$RNA1,CHEM1,1:R1-1:R2$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,8:pair-20:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,20:pair-8:pair$RNA3{ss}$";
			String[] notations = ComplexNotationParser.decompose(notation);
			//assertEquals("Decomp0", notations[0],"RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|CHEM1{CM#2}$RNA1,CHEM1,1:R1-1:R2$$$");
			assertEquals("Decomp1", notations[1],"RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$");
			assertEquals("Decomp2", notations[2],"RNA3{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}$$$RNA3{ss}$");
			
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}$RNA1,RNA1,21:R2-1:R1$$$"; //cycle + linear mixture
			notations = ComplexNotationParser.decompose(notation);
			assertEquals("Decomp0", "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}$RNA1,RNA1,21:R2-1:R1$$$", notations[0]);
			assertEquals("Decomp1", "RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}$$$$", notations[1]);
			
			notation = "RNA1{R(A)}|RNA2{R(A)}|CHEM1{CovX-3}|CHEM2{CovX-2}$RNA1,CHEM2,1:R1-1:R1|RNA1,CHEM1,1:R2-1:R1$$$"; //cycle + linear mixture
			notations = ComplexNotationParser.decompose(notation);
			//assertEquals("Decomp0", "RNA1{R(A)}|CHEM2{CM#4}|CHEM1{CM#3}$RNA1,CHEM2,1:R1-1:R1|RNA1,CHEM1,1:R2-1:R1$$$", notations[0]);
			assertEquals("Decomp1", "RNA2{R(A)}$$$$", notations[1]);
			
			notation = "PEPTIDE1{A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE2{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE3{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE4{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C}$PEPTIDE2,PEPTIDE4,18:R3-18:R3$$PEPTIDE1{hc}$"; //cycle + linear mixture
			notations = ComplexNotationParser.decompose(notation);
			assertEquals("Decomp0", "PEPTIDE1{A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}$$$PEPTIDE1{hc}$", notations[0]);
			assertEquals("Decomp1", "PEPTIDE2{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}|PEPTIDE4{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C}$PEPTIDE2,PEPTIDE4,18:R3-18:R3$$$", notations[1]);
			assertEquals("Decomp2", "PEPTIDE3{A.A.A.A.A.A.L.L.[seC].K.K.K.K.L.A.A.K.C.C.A.A.A.A.A.A.L.L.K.K.K.K.L.A.A.K.C.C}$$$$", notations[2]);
			

			
			//Parse complex notation to sense and antisense notations
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$";
			String allNodeNotation = ComplexNotationParser.getAllNodeString(notation);
			List<PolymerNode> pnList = ComplexNotationParser.getPolymerNodeList(allNodeNotation);
			assertEquals("Count", 2, pnList.size());
			List<RNAPolymerNode> rnaList = ComplexNotationParser.getRNAPolymerNodeList(notation);
			assertEquals("Compare", pnList.size(), rnaList.size());
			RNAPolymerNode rna = rnaList.get(0);
			String annotation = rna.getAnotation();
			assertEquals("Check annotation", annotation, annotation.trim());
			assertEquals("SS notation", "R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)", rna.getLabel());
			rna = rnaList.get(1);
			annotation = rna.getAnotation();
			assertEquals("Check annotation", annotation, annotation.trim());
			assertEquals("AS notation", "R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)", rna.getLabel());
			
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA1{ss}$";
			allNodeNotation = ComplexNotationParser.getAllNodeString(notation);
			pnList = ComplexNotationParser.getPolymerNodeList(allNodeNotation);
			assertEquals("Count", 2, pnList.size());
			rnaList = ComplexNotationParser.getRNAPolymerNodeList(notation);
			assertEquals("Compare", pnList.size(), rnaList.size());
			rna = rnaList.get(0);
			annotation = rna.getAnotation();
			assertEquals("Check annotation", annotation, annotation.trim());
			assertEquals("SS notation", "R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)", rna.getLabel());
			rna = rnaList.get(1);
			annotation = rna.getAnotation();
			assertNull("Check annotation",annotation);
			
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}$$$RNA1{as}$";
			allNodeNotation = ComplexNotationParser.getAllNodeString(notation);
			pnList = ComplexNotationParser.getPolymerNodeList(allNodeNotation);
			assertEquals("Count", 1, pnList.size());
			rnaList = ComplexNotationParser.getRNAPolymerNodeList(notation);
			assertEquals("Compare", pnList.size(), rnaList.size());
			rna = rnaList.get(0);
			annotation = rna.getAnotation();
			assertEquals("Check annotation", annotation, annotation.trim());
			assertEquals("SS notation", "R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)", rna.getLabel());
			

			
			//replace monomers
			notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$";
			System.out.println("Replace P with sP for: " + notation);
			String result = ComplexNotationParser.replaceMonomer(notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE, "P", "sP");
			System.out.println("Result: " + result);

			String[] rtcNotations = {"R(A)P.R(G)", "A.A.G.L", "PEG2", "RNA1{P.R(A)P.R(G)}$$$$"};
			for (String rtcNotation : rtcNotations) {
				String cn = ComplexNotationParser.standardize(rtcNotation);
				System.out.println(rtcNotation + " => " + cn);
			}

			
			//System.out.println("Hybridization bug fix:");
			notation = "RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)[sP].R(U)}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}$$$RNA1{ss}|RNA2{as}$";
			hybridizedNotation = ComplexNotationParser.hybridize(notation);
			assertEquals("Hybridization bug fix:","RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)[sP].R(U)}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}$$RNA1,RNA2,4:pair-62:pair|RNA1,RNA2,7:pair-59:pair|RNA1,RNA2,10:pair-56:pair|RNA1,RNA2,13:pair-53:pair|RNA1,RNA2,16:pair-50:pair|RNA1,RNA2,19:pair-47:pair|RNA1,RNA2,22:pair-44:pair|RNA1,RNA2,25:pair-41:pair|RNA1,RNA2,28:pair-38:pair|RNA1,RNA2,31:pair-35:pair|RNA1,RNA2,34:pair-32:pair|RNA1,RNA2,37:pair-29:pair|RNA1,RNA2,40:pair-26:pair|RNA1,RNA2,43:pair-23:pair|RNA1,RNA2,46:pair-20:pair|RNA1,RNA2,49:pair-17:pair|RNA1,RNA2,52:pair-14:pair|RNA1,RNA2,55:pair-11:pair|RNA1,RNA2,58:pair-8:pair|RNA1,RNA2,61:pair-5:pair|RNA1,RNA2,64:pair-2:pair$RNA1{ss}|RNA2{as}$", hybridizedNotation);
			
			//System.out.println("Test for modification:");
			notation = "RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)[sP].R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)[sP].R(U)}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}$$$RNA1{ss}|RNA2{as}$";
			List<PolymerNode> polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			boolean hasMod = ComplexNotationParser.hasNucleotideModification(polyNodeList);
			assertTrue(hasMod);
			
			
			//System.out.println("Test for modification:");
			notation = "RNA1{R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P.R(A)P.R(U)P}|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P.R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P.R(A)P.R(C)P}$$$RNA1{ss}|RNA2{as}$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser.hasNucleotideModification(polyNodeList);
			assertFalse(hasMod);

			
			//System.out.println("Test for modification:");
			notation = "RNA1{R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(C)}|RNA2{R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)}$$RNA1,RNA2,2:pair-56:pair|RNA1,RNA2,5:pair-53:pair|RNA1,RNA2,8:pair-50:pair|RNA1,RNA2,11:pair-47:pair|RNA1,RNA2,14:pair-44:pair|RNA1,RNA2,17:pair-41:pair|RNA1,RNA2,20:pair-38:pair|RNA1,RNA2,23:pair-35:pair|RNA1,RNA2,26:pair-32:pair|RNA1,RNA2,29:pair-29:pair|RNA1,RNA2,32:pair-26:pair|RNA1,RNA2,35:pair-23:pair|RNA1,RNA2,38:pair-20:pair|RNA1,RNA2,41:pair-17:pair|RNA1,RNA2,44:pair-14:pair|RNA1,RNA2,47:pair-11:pair|RNA1,RNA2,50:pair-8:pair|RNA1,RNA2,53:pair-5:pair|RNA1,RNA2,56:pair-2:pair$RNA1{ss}|RNA2{as}$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser.hasNucleotideModification(polyNodeList);
			assertTrue(hasMod);

			
			//System.out.println("Test for modification:");
			notation = "RNA1{R(A)P.R(U)P.R(C)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)}|RNA2{R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.[dR](A)P.[dR](T)}$$RNA1,RNA2,2:pair-74:pair|RNA1,RNA2,5:pair-71:pair|RNA1,RNA2,8:pair-68:pair|RNA1,RNA2,11:pair-65:pair|RNA1,RNA2,14:pair-62:pair|RNA1,RNA2,17:pair-59:pair|RNA1,RNA2,20:pair-56:pair|RNA1,RNA2,23:pair-53:pair|RNA1,RNA2,26:pair-50:pair|RNA1,RNA2,29:pair-47:pair|RNA1,RNA2,32:pair-44:pair|RNA1,RNA2,35:pair-41:pair|RNA1,RNA2,38:pair-38:pair|RNA1,RNA2,41:pair-35:pair|RNA1,RNA2,44:pair-32:pair|RNA1,RNA2,47:pair-29:pair|RNA1,RNA2,50:pair-26:pair|RNA1,RNA2,53:pair-23:pair|RNA1,RNA2,56:pair-20:pair|RNA1,RNA2,59:pair-17:pair|RNA1,RNA2,62:pair-14:pair|RNA1,RNA2,65:pair-11:pair|RNA1,RNA2,68:pair-8:pair|RNA1,RNA2,71:pair-5:pair|RNA1,RNA2,74:pair-2:pair$$";
			polyNodeList = ComplexNotationParser.getPolymerNodeList(notation);
			hasMod = ComplexNotationParser.hasNucleotideModification(polyNodeList);
			assertTrue(hasMod);
			
			

		} catch  (Exception ex){
			fail("testNotationManipulation");
		}
	}




	private static boolean testComplexNotationValidity(String notation) {
		try {
			System.out.println("Testing Notation Validity for: " + notation);
			System.out.println("Notation Format Valid? " + ComplexNotationParser.validateNotationFormat(notation));
			System.out.println("Node String: " + ComplexNotationParser.getAllNodeString(notation));
			System.out.println("Edge String: " + ComplexNotationParser.getAllEdgeString(notation));
			System.out.println("Base Pair String: " + ComplexNotationParser.getAllBasePairString(notation));
			System.out.println("Node Label String: " + ComplexNotationParser.getAllNodeLabelString(notation));
			System.out.println("Others String: " + ComplexNotationParser.getOtherString(notation));
			System.out.println("Notation Valid? " + ComplexNotationParser.validateComplexNotation(notation));
			System.out.println("*********************");
			return true;
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
			return false;
		}
	}

	private static boolean testGetCanonicalNotation(String notation) {
		try {
			System.out.println("Testing getCanonicalNotation for: " + notation);
			System.out.println("Canonical Notation: " + ComplexNotationParser.getCanonicalNotation(notation));
			System.out.println("*********************");
			return true;
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
			return false;
		}
	}

	private static boolean testGetMoleculeInfo(String notation) {
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
			return true;
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
			return false;
		}
	}



	private static boolean testGetMoleculeInfoViaSmiles(String notation) {
		System.out.println("testGetMoleculeInfoViaSmiles...");
		try {
			long start = System.currentTimeMillis();
			System.out.println("Start Time: " + start);
			String smiles = ComplexNotationParser.getComplexPolymerSMILES(notation);
			MoleculeInfo mi = StructureParser.getMoleculeInfo(smiles);
			long end = System.currentTimeMillis();
			System.out.println("End Time: " + end);
			System.out.println("MW = " + mi.getMolecularWeight());
			System.out.println("MF = " + mi.getMolecularFormula());
			System.out.println("Mass = " + mi.getExactMass());
			System.out.println("Time = " + (end - start));
			return true;
		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
			return false;
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


	@Test
	public void testInlineNotation() throws NotationException, IOException, MonomerException, StructureException, JDOMException  {
		
		String notation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		String smiles = ComplexNotationParser.getComplexPolymerSMILES(notation);

		//replaced A with Smiles String
		notation = "PEPTIDE1{[C[C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|].G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		String smilesInline = ComplexNotationParser.getComplexPolymerSMILES(notation);
		assertEquals(smiles,smilesInline);
		
		
		boolean valid=ComplexNotationParser.validateComplexNotation(notation);
		assertTrue(valid);
		valid=ComplexNotationParser.validateNotationFormat(notation);
		assertTrue(valid);
		
		
		//replaced A with slightly modified A
		notation = "PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		valid=ComplexNotationParser.validateNotationFormat(notation);
		assertTrue(valid);
		valid=ComplexNotationParser.validateComplexNotation(notation);
		assertTrue(valid);
		
		smilesInline = ComplexNotationParser.getComplexPolymerSMILES(notation);
		assertEquals("[H]NCCCC[C@H](NC(=O)[C@H](CS[H])NC(=O)[C@H](CS[H])NC(=O)CNC(=O)CNC(=O)CNCC(=O)[C@H](C)N[H])C(=O)N[C@@H](CCCCN[H])C(=O)N[C@@H](CCCCN[H])C(=O)N[C@@H](CCCCNC(=O)C1CCC(CN2C(=O)C=CC2=O)CC1)C(O)=O",smilesInline);
		
		
	}
	
	
	@Test
	public void testReplaceSmiles() throws NotationException, MonomerException, JDOMException, IOException{
	

		
		String helmNotation="PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
		String notationNoSmiles=ComplexNotationParser.getNotationByReplacingSmiles(helmNotation, MonomerFactory.getInstance().getMonomerStore());
		
		assertEquals(notationNoSmiles, "PEPTIDE1{[AM#1].G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$");
		
		
	}

}
