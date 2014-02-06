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
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.demo.tools.ComplexNotationSample;
import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;
import org.helm.notation.model.RNAPolymerNode;

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.*;

import org.junit.Test;

public class ComplexNotationParserTest {

	@Test
	public void testGenericConnection() {
		try {
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();
		} catch  (Exception ex) {
            fail("Initialization");
        }
        
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
		
		try {
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();
		} catch  (Exception ex) {
            fail("Initialization");
        }
		
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

        //chemical connection
        notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|CHEM1{CovX-2}$RNA1,CHEM1,1:R1-1:R1$$$";
        assertTrue(testGetCanonicalNotation(notation));

        notation = "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|PEPTIDE1{A.G.G.G.K.K.K.K}|CHEM1{CovX-2}|CHEM2{3Bio}$PEPTIDE1,CHEM2,8:R2-1:R1|RNA1,CHEM1,1:R1-1:R1$$$";
        assertTrue(testGetCanonicalNotation(notation));
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
        } catch (Exception ex) {
            fail("testComplexNotationValidity: " + notation);
        }
        return true;
    }
	
    private static boolean testGetCanonicalNotation(String notation) {
        try {
            System.out.println("Testing getCanonicalNotation for: " + notation);
            System.out.println("Canonical Notation: " + ComplexNotationParser.getCanonicalNotation(notation));
            System.out.println("*********************");
        } catch (Exception ex) {
            fail("testGetCanonicalNotation: " + notation);
        }
        return true;
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
        } catch (Exception ex) {
            Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
        }
        return true;
    }



}
