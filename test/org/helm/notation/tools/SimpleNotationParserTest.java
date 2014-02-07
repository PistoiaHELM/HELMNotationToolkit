package org.helm.notation.tools;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import org.jdom.JDOMException;
import org.junit.Test;

import chemaxon.marvin.plugin.PluginException;
import junit.framework.*;

public class SimpleNotationParserTest extends TestCase {

	
	public String getSimpleRNANotation(){		
		return "P.R(A)[sP].RP.R(G)P.[LR]([5meC])";
	}
	
	public String getSimplePeptideNotation(){		
		return "G.G.K.A.A.[seC]";
	}
	
	public String getSimpleChemNotation(){
		return "PEG2";		
	}
	
	public String getSmilesNotation(){
		return "[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|";		
	}
	
	public void setUp() throws MonomerException, IOException, JDOMException, NotationException{
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();	
	}
	
	
	@Test
	public void testMonomerCount() throws NotationException {
		
		    
        int count=SimpleNotationParser.getMonomerCount(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
        assertEquals(11,count);
        
        count=SimpleNotationParser.getMonomerCount(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
        assertEquals(6,count);
        
        count=SimpleNotationParser.getMonomerCount(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals(1,count);
        
        count=SimpleNotationParser.getMonomerCount(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals(1,count);
        
        
	    
	}
	@Test
	public void testCanonicalNotation() throws NotationException, MonomerException, StructureException, JDOMException, IOException {
		
		String canonicalNotation=SimpleNotationParser.getSimpleCanonicalNotation(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
        assertEquals("R(A)[sP].RP.R(G)P.[LR]([5meC])P",canonicalNotation);
        Map.Entry entry = SimpleNotationParser.getSimpleCanonicalNotationMapEntry(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
        assertEquals("1",entry.getKey().toString());
        assertEquals("R(A)[sP].RP.R(G)P.[LR]([5meC])P",entry.getValue());
        
        
        
		canonicalNotation=SimpleNotationParser.getSimpleCanonicalNotation(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
        assertEquals("A.A.[seC].G.G.K",canonicalNotation);
        entry = SimpleNotationParser.getSimpleCanonicalNotationMapEntry(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
        assertEquals("3",entry.getKey().toString());
        assertEquals("A.A.[seC].G.G.K",entry.getValue());
        
        
        canonicalNotation=SimpleNotationParser.getSimpleCanonicalNotation(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("PEG2",canonicalNotation);
        entry = SimpleNotationParser.getSimpleCanonicalNotationMapEntry(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("0",entry.getKey().toString());
        assertEquals("PEG2",entry.getValue());

        
        canonicalNotation=SimpleNotationParser.getSimpleCanonicalNotation(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|",canonicalNotation);
        entry = SimpleNotationParser.getSimpleCanonicalNotationMapEntry(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("0",entry.getKey().toString());
        assertEquals("CM#1",entry.getValue());
        
        
        
	}
	@Test
	public void testGetComplexNotation() throws NotationException, MonomerException, StructureException, JDOMException, IOException{
		
		String complexNotation=SimpleNotationParser.getComplexNotation(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
        assertEquals("RNA1{P.R(A)[sP].RP.R(G)P.[LR]([5meC])}$$$$",complexNotation);
        
		
        complexNotation=SimpleNotationParser.getComplexNotation(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
        assertEquals("PEPTIDE1{G.G.K.A.A.[seC]}$$$$",complexNotation);
        
        complexNotation=SimpleNotationParser.getComplexNotation(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("CHEM1{PEG2}$$$$",complexNotation);
        
        
        complexNotation=SimpleNotationParser.getComplexNotation(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals("CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$$$$",complexNotation);
        
        
	}

	@Test
	public void testReplaceMonomer() throws MonomerException, IOException, JDOMException, NotationException{
		
		String result=SimpleNotationParser.replaceMonomer(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE, "P", "sP");
		assertEquals("[sP].R(A)[sP].R[sP].R(G)[sP].[LR]([5meC])",result);
		
		result=SimpleNotationParser.replaceMonomer(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE,"A", "Q");
		assertEquals("G.G.K.Q.Q.[seC]",result);
		
		result=SimpleNotationParser.replaceMonomer(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE, "PEG2", "SS3");
		assertEquals("SS3",result);
		
		
	}
	
	@Test
	public void testGetSmiles() throws IOException, NotationException, StructureException, MonomerException, JDOMException{
		
		String result=SimpleNotationParser.getSimplePolymerSMILES(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
		assertEquals("Cc1cn([C@H]2O[C@@]3(COP(O)(=O)O[C@@H]4[C@@H](COP(O)(=O)O[C@@H]5[C@@H](COP(S)(=O)O[C@@H]6[C@@H](COP(O)([*])=O)O[C@@H]([C@@H]6O)n6cnc7c(N)ncnc67)O[C@@H]([*])[C@@H]5O)O[C@@H]([C@@H]4O)n4cnc5c4nc(N)[nH]c5=O)CO[C@@H]2[C@@H]3O[*])c(=O)nc1N |r,$;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;_R1;;;;;;;;;;;;;;;;;;_R3;;;;;;;;;;;;;;;;;;;;;;;_R2;;;;;$|",result);
		
		result=SimpleNotationParser.getSimplePolymerSMILES(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
		assertEquals("C[C@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCN[*])NC(=O)CNC(=O)CN[*])C(=O)N[C@@H](C[SeH])C([*])=O |r,$;;;;;;;;;;;;;;;;_R3;;;;;;;;;;_R1;;;;;;;;_R2;$|",result);
		
		result=SimpleNotationParser.getSimplePolymerSMILES(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
		assertEquals("[*]OCCOCCO[*] |$_R1;;;;;;;;_R2$|",result);
		
		result=SimpleNotationParser.getSimplePolymerSMILES(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
		assertEquals(getSmilesNotation(),result);
		
	}
	@Test
	public void testValidation() throws IOException, NotationException, MonomerException, StructureException, JDOMException{
		
		assertTrue(SimpleNotationParser.validateSimpleNotation(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE));
		
		assertTrue(SimpleNotationParser.validateSimpleNotation(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE));
		
		assertTrue(SimpleNotationParser.validateSimpleNotation(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE));
		
		assertTrue(SimpleNotationParser.validateSimpleNotation(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE));		
		
	}
	
	@Test
	public void testGetMoleculeInfo() throws NotationException, MonomerException, IOException, JDOMException, PluginException, StructureException{
				
        MoleculeInfo mi=SimpleNotationParser.getMoleculeInfo(getSimpleRNANotation(), Monomer.NUCLIEC_ACID_POLYMER_TYPE);
        assertEquals(1251.805,mi.getMolecularWeight(),1e-15);
        assertEquals("C36H49N13O27P4S1",mi.getMolecularFormula());
        assertEquals(1251.153200165,mi.getExactMass(),1e-15);
        
        mi=SimpleNotationParser.getMoleculeInfo(getSimplePeptideNotation(), Monomer.PEPTIDE_POLYMER_TYPE);
        assertEquals(552.48,mi.getMolecularWeight(),1e-15);
        assertEquals("C19H35N7O7Se1",mi.getMolecularFormula());
        assertEquals(553.176318337,mi.getExactMass(),1e-15);
        
        
        mi=SimpleNotationParser.getMoleculeInfo(getSimpleChemNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals(106.1204,mi.getMolecularWeight(),1e-15);
        assertEquals("C4H10O3",mi.getMolecularFormula());
        assertEquals(106.062994186,mi.getExactMass(),1e-15);
        
        
        mi=SimpleNotationParser.getMoleculeInfo(getSmilesNotation(), Monomer.CHEMICAL_POLYMER_TYPE);
        assertEquals(149.16506,mi.getMolecularWeight(),1e-15);
        assertEquals("C6H13O4",mi.getMolecularFormula());
        assertEquals(149.081383904,mi.getExactMass(),1e-15);
        
        
	}
	
	
	public void testNucleotideFunctions() throws NotationException, MonomerException, IOException, JDOMException, StructureException{
      String notation = "R(C)P.R(G)P.R(A)P.R(U)P.R(A)P.R(U)P.R(G)P.R(G)P.R(G)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.[dR](U)P.[dR](U)";
      System.out.println("getNucleotideList Start: " + System.currentTimeMillis());
      List<Nucleotide> nucleotideList = SimpleNotationParser.getNucleotideList(notation);
      System.out.println("getNucleotideList End: " + System.currentTimeMillis());

      System.out.println("getNucleotideList wo validation Start : " + System.currentTimeMillis());
      nucleotideList = SimpleNotationParser.getNucleotideList(notation, false);
      System.out.println("getNucleotideList wo validation End: " + System.currentTimeMillis());

      for (int i = 0; i < nucleotideList.size(); i++) {
          Nucleotide nuc = nucleotideList.get(i);
          System.out.println("Symbol: " + nuc.getSymbol());
          System.out.println("Modified: " + nuc.isModified());
          System.out.println("Notation: " + nuc.getNotation());
      }
      
      assertEquals(21,nucleotideList.size());
      
      String nucleotideSeq =SimpleNotationParser.getNucleotideSequence(notation);
      System.out.println("Sequence: " + nucleotideSeq);
      
      assertEquals("CGAUAUGGGCUGAAUACAAUU",nucleotideSeq);

      notation = "K.C.C.C.W.K.[seC]";
      System.out.println("Peptide Sequence: " + SimpleNotationParser.getPeptideSequence(notation));
      assertEquals("KCCCWKC",SimpleNotationParser.getPeptideSequence(notation));
      System.out.println("Modified Peptide Sequence: " + SimpleNotationParser.getModifiedPeptideSequence(notation, "|"));
      assertEquals("K|C|C|C|W|K|seC",SimpleNotationParser.getModifiedPeptideSequence(notation, "|"));


      notation = "R(A)[sP].R(G)P.R(C).PEG.[LR]([5meC])[sP].R(PEG)";
      //valid = SimpleNotationParser.getStrictNucleotideList(notation, true);
      System.out.println("Input sequence " + notation);
      nucleotideList = SimpleNotationParser.getStrictNucleotideList(notation, false);
      StringBuffer outputList = new StringBuffer();
      for (int i = 0; i < nucleotideList.size(); i++) {
          Nucleotide nuc = nucleotideList.get(i);
          System.out.println("Symbol: " + nuc.getSymbol());
          System.out.println("Modified: " + nuc.isModified());
          System.out.println("Notation: " + nuc.getNotation());
          outputList.append(nuc.getSymbol());
      }
      assertEquals(6,nucleotideList.size());
      System.out.println("Output List: " + outputList.toString());
      assertEquals("modAGmodCXmodCendX",outputList.toString());
		
	}
	
	
	 
	    


	    

}
