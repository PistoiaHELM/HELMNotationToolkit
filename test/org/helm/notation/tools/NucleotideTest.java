package org.helm.notation.tools;

import static org.junit.Assert.*;

import java.io.IOException;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import org.jdom.JDOMException;
import org.junit.Test;

public class NucleotideTest {

	@Test
	public void testGetNotation() throws NotationException, IOException,
			JDOMException, StructureException, MonomerException {
		String notation = NucleotideSequenceParser.getNotation("AUG");
		assertEquals("R(A)P.R(U)P.R(G)P", notation);

		String smiles = SimpleNotationParser.getSimplePolymerSMILES(notation,
				"RNA");
		assertEquals(
				"Nc1nc2n(cnc2c(=O)[nH]1)[C@H]1O[C@H](COP(O)(=O)O[C@@H]2[C@@H](COP(O)(=O)O[C@@H]3[C@@H](CO[*])O[C@@H]([C@@H]3O)n3cnc4c(N)ncnc34)O[C@@H]([C@@H]2O)n2ccc(=O)[nH]c2=O)[C@@H](OP(O)([*])=O)[C@H]1O |$;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;_R1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;_R2;;;$|",
				smiles);

		String complexNotation = "RNA1{" + notation + "}$$$$";
		smiles = ComplexNotationParser.getComplexPolymerSMILES(complexNotation);
		assertEquals(
				"[H]OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c1nc(N)[nH]c2=O)n1ccc(=O)[nH]c1=O)n1cnc2c(N)ncnc12",
				smiles);

	}

	@Test
	public void testIsInOppositeDirection() throws NotationException,
			IOException, JDOMException {
		boolean opp = NucleotideSequenceParser.isInOppositeDirection("AUGTAU",
				"AUACAU");
		assertFalse(opp);
		if (opp) {
			System.out.println("opposite direction");
		} else {
			System.out.println("same direction");
		}

		opp = NucleotideSequenceParser.isInOppositeDirection("ATGCAAA",
				"UUUGCAU");
		assertFalse(opp);
		if (opp) {
			System.out.println("opposite direction");
		} else {
			System.out.println("same direction");
		}

	}

	@Test
	public void testGetSirnaNotation() throws NotationException, IOException,
			JDOMException {
		String senseSeq = "CGAAAUGUUCAUACUGUUGdTdT";
		String antiSenseSeq = "UUACAAUUUGGACUUUCCGdTdT";

		String siNotation = NucleotideSequenceParser.getSirnaNotation(senseSeq,
				antiSenseSeq);
		System.out.println(siNotation);
		assertEquals(
				"RNA1{R(C)P.R(G)P.R(A)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(U)P.R(C)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(G)P.R(U)P.R(U)P.R(G)P.[dR](T)P.[dR](T)}|RNA2{R(U)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.R(A)P.R(C)P.R(U)P.R(U)P.R(U)P.R(C)P.R(C)P.R(G)P.[dR](T)P.[dR](T)}$$RNA1,RNA2,5:pair-50:pair|RNA1,RNA2,8:pair-47:pair|RNA1,RNA2,11:pair-44:pair|RNA1,RNA2,14:pair-41:pair$RNA1{ss}|RNA2{as}$",
				siNotation);

		siNotation = NucleotideSequenceParser.getSirnaNotation(senseSeq,
				antiSenseSeq,
				NucleotideSequenceParser.RNA_DESIGN_TUSCHL_19_PLUS_2);
		assertEquals(
				"RNA1{R(C)P.R(G)P.R(A)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(U)P.R(C)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(G)P.R(U)P.R(U)P.R(G)P.[dR](T)P.[dR](T)}|RNA2{R(U)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.R(A)P.R(C)P.R(U)P.R(U)P.R(U)P.R(C)P.R(C)P.R(G)P.[dR](T)P.[dR](T)}$$RNA1,RNA2,2:pair-56:pair|RNA1,RNA2,5:pair-53:pair|RNA1,RNA2,11:pair-47:pair|RNA1,RNA2,14:pair-44:pair|RNA1,RNA2,20:pair-38:pair|RNA1,RNA2,23:pair-35:pair|RNA1,RNA2,29:pair-29:pair|RNA1,RNA2,32:pair-26:pair|RNA1,RNA2,38:pair-20:pair|RNA1,RNA2,44:pair-14:pair|RNA1,RNA2,47:pair-11:pair|RNA1,RNA2,50:pair-8:pair$RNA1{ss}|RNA2{as}$",
				siNotation);
		System.out.println(siNotation);

	}

	@Test
	public void testGetNaturalAnalog() throws MonomerException, IOException,
			JDOMException {
		String notation = "[LR]([5meC])[sP]";
		Nucleotide n = new Nucleotide("ls5C", notation);
		
		assertEquals("C", n.getNaturalAnalog());
		String base = notation.substring(notation.indexOf("(") + 1,
				notation.indexOf(")"));
		base = base.replaceAll("\\[|\\]", "");
		Monomer baseMonomer = MonomerFactory.getInstance().getMonomerDB()
				.get(Monomer.NUCLIEC_ACID_POLYMER_TYPE).get(base);
		
		assertEquals("C", baseMonomer.getNaturalAnalog());
		
		
		n=new Nucleotide("R(A)P",Nucleotide.MIDDLE_POSITION_TYPE);
    	
    	assertEquals("A", n.getNaturalAnalog());
    	n=new Nucleotide("R([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",Nucleotide.MIDDLE_POSITION_TYPE);    	
    	assertEquals("A", n.getNaturalAnalog());
    	
    	n=new Nucleotide("[CO[C@H]1[C@H]([*])O[C@H](CO[*])[C@H]1O[*] |$;;;;_R3;;;;;_R1;;;_R2$|]([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",Nucleotide.MIDDLE_POSITION_TYPE);    	
    	assertEquals("A", n.getNaturalAnalog());
    	n=new Nucleotide("R([Cc1nc2c(nc(N)[nH]c2=O)n1[*] |$;;;;;;;;;;;;_R1$|])P",Nucleotide.MIDDLE_POSITION_TYPE);    	
    	assertEquals("X", n.getNaturalAnalog());
    	n=new Nucleotide("RP",Nucleotide.MIDDLE_POSITION_TYPE);    	
    	assertEquals("X", n.getNaturalAnalog());
    	
    	n=new Nucleotide("R",Nucleotide.ENDING_POSITION_TYPE);    	
    	
    	assertEquals("X", n.getNaturalAnalog());
    	n=new Nucleotide("([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",Nucleotide.STARTING_POSITION_TYPE);    	
    	assertEquals("A", n.getNaturalAnalog());
    	
    
	}

	@Test
	public void testGetComplementSequence() throws NotationException, IOException, JDOMException {

		String sequence = "5'-UAU GUC UCC AGA AUG UAG CdTdT-3'";

		String normal = NucleotideSequenceParser
				.getNormalComplementSequence(sequence);
		assertEquals("5'-AAGCUACAUUCUGGAGACAUA-3'", normal);
		String reverse = NucleotideSequenceParser
				.getReverseComplementSequence(sequence);
		assertEquals("3'-AUACAGAGGUCUUACAUCGAA-5'", reverse);
		String notation = NucleotideSequenceParser.getNotation(sequence);
		assertEquals(
				"R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P.R(C)P.[dR](T)P.[dR](T)P",
				notation);

	
		sequence = "5'-mGCU ACA UUC UGG AGA CAU AdTdT-3'";

		normal = NucleotideSequenceParser.getNormalComplementSequence(sequence);
		assertEquals("5'-AAUAUGUCUCCAGAAUGUAGC-3'", normal);
		reverse = NucleotideSequenceParser
				.getReverseComplementSequence(sequence);
		assertEquals("3'-CGAUGUAAGACCUCUGUAUAA-5'", reverse);
		notation = NucleotideSequenceParser.getNotation(sequence);
		assertEquals(
				"[mR](G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P.R(A)P.[dR](T)P.[dR](T)P",
				notation);

		sequence = "agcuagggu";

		normal = NucleotideSequenceParser.getNormalComplementSequence(sequence);
		assertEquals("5'-ACCCUAGCU-3'", normal);
		reverse = NucleotideSequenceParser
				.getReverseComplementSequence(sequence);
		assertEquals("3'-UCGAUCCCA-5'", reverse);
		notation = NucleotideSequenceParser.getNotation(sequence);
		assertEquals(
				"R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(G)P.R(G)P.R(U)P",
				notation);



//		Map<String, Map<String, String>> templates = NucleotideFactory
//				.getInstance().getNucleotideTemplates();
//		String templatesXML = NucleotideSequenceParser
//				.getNucleotideTemplatesXML(templates);
//		System.out.println(templatesXML);
//		NucleotideFactory.getInstance().saveNucleotideTemplates();

	}
}
