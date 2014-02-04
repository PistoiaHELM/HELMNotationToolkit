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

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.helm.notation.MonomerFactory;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;
import org.helm.notation.tools.ComplexNotationParser;
import org.helm.notation.tools.MonomerParser;
import org.helm.notation.tools.SimpleNotationParser;
import org.jdom.output.XMLOutputter;

public class xHELMSample {

	public static void main(String[] args) {
		// Load HELM and convert to xHELM
		try {
			MonomerFactory.getInstance();

			String notation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
			//String notation = "RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}|PEPTIDE1{A.C.G.K.E.D.K.R}|CHEM1{SMCC}$PEPTIDE1,CHEM1,2:R3-1:R2|RNA1,CHEM1,1:R1-1:R1$$$";
			//String notation = "PEPTIDE1{Q.A.Q.L.Q.E.S.G.P.G.L.A.K.P.S.E.T.L.S.L.T.C.T.A.S.G.G.S.I.S.G.Y.Y.W.S.W.I.R.Q.P.V.G.K.G.L.E.W.I.G.R.I.Y.T.S.G.S.T.N.Y.N.P.S.L.K.S.R.A.T.M.S.A.D.T.S.K.N.Q.F.S.L.K.L.S.S.A.T.V.V.D.T.V.A.Y.Y.C.V.R.G.R.F.T.Y.F.D.Y.W.G.Q.G.T.L.A.T.A.S.S.V.S.T.K.G.P.S.A.F.P.L.V.P.S.S.K.S.T.S.G.G.T.V.V.L.G.C.L.A.K.D.Y.F.P.E.P.A.T.A.S.W.N.S.G.V.L.T.S.G.A.H.T.F.P.V.A.L.Q.S.S.G.L.Y.S.L.S.S.A.A.T.A.P.S.S.S.L.G.T.Q.T.Y.I.C.N.A.N.H.K.P.S.N.T.K.A.D.K.K.A.E.P.K.S.C.D.K.T.H.T.C.P.P.C.P.V.P.E.L.L.G.G.P.S.A.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.A.T.C.A.A.A.D.A.S.H.E.D.P.E.A.K.F.N.W.Y.A.D.G.A.E.A.H.N.V.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.A.A.S.A.L.T.A.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.A.S.N.K.V.L.P.V.P.I.E.K.T.I.S.K.V.K.G.Q.P.R.E.P.Q.A.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.A.S.L.T.C.L.A.K.G.F.Y.P.S.D.I.V.A.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.A.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.A.D.K.S.R.W.Q.Q.G.N.A.F.S.C.S.A.M.H.E.V.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE2{Q.A.Q.L.Q.E.S.G.P.G.L.A.K.P.S.E.T.L.S.L.T.C.T.A.S.G.G.S.I.S.G.Y.Y.W.S.W.I.R.Q.P.V.G.K.G.L.E.W.I.G.R.I.Y.T.S.G.S.T.N.Y.N.P.S.L.K.S.R.A.T.M.S.A.D.T.S.K.N.Q.F.S.L.K.L.S.S.A.T.V.V.D.T.V.A.Y.Y.C.V.R.G.R.F.T.Y.F.D.Y.W.G.Q.G.T.L.A.T.A.S.S.V.S.T.K.G.P.S.A.F.P.L.V.P.S.S.K.S.T.S.G.G.T.V.V.L.G.C.L.A.K.D.Y.F.P.E.P.A.T.A.S.W.N.S.G.V.L.T.S.G.A.H.T.F.P.V.A.L.Q.S.S.G.L.Y.S.L.S.S.A.A.T.A.P.S.S.S.L.G.T.Q.T.Y.I.C.N.A.N.H.K.P.S.N.T.K.A.D.K.K.A.E.P.K.S.C.D.K.T.H.T.C.P.P.C.P.V.P.E.L.L.G.G.P.S.A.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.A.T.C.A.A.A.D.A.S.H.E.D.P.E.A.K.F.N.W.Y.A.D.G.A.E.A.H.N.V.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.A.A.S.A.L.T.A.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.A.S.N.K.V.L.P.V.P.I.E.K.T.I.S.K.V.K.G.Q.P.R.E.P.Q.A.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.A.S.L.T.C.L.A.K.G.F.Y.P.S.D.I.V.A.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.A.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.A.D.K.S.R.W.Q.Q.G.N.A.F.S.C.S.A.M.H.E.V.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE3{E.I.A.L.T.Q.S.P.V.T.L.S.L.S.P.G.E.R.V.T.L.S.C.R.V.S.Q.I.A.S.S.V.Y.L.V.W.Y.Q.Q.K.P.G.Q.V.P.R.L.L.M.F.G.S.S.S.R.V.T.G.I.P.D.R.F.S.G.S.G.S.G.T.D.F.T.L.T.I.S.R.L.E.P.E.D.F.V.A.Y.Y.C.Q.Q.Y.G.S.S.Q.G.T.F.G.P.G.T.K.A.D.I.K.R.T.A.V.V.P.S.A.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.V.S.A.A.C.L.L.N.N.F.Y.P.R.E.V.K.A.Q.W.K.A.D.N.V.L.Q.S.G.N.S.Q.E.S.A.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.V.D.Y.E.K.H.K.A.Y.V.C.E.A.T.H.Q.G.L.S.S.P.A.T.K.S.F.N.R.G.E.C}|PEPTIDE4{E.I.A.L.T.Q.S.P.V.T.L.S.L.S.P.G.E.R.V.T.L.S.C.R.V.S.Q.I.A.S.S.V.Y.L.V.W.Y.Q.Q.K.P.G.Q.V.P.R.L.L.M.F.G.S.S.S.R.V.T.G.I.P.D.R.F.S.G.S.G.S.G.T.D.F.T.L.T.I.S.R.L.E.P.E.D.F.V.A.Y.Y.C.Q.Q.Y.G.S.S.Q.G.T.F.G.P.G.T.K.A.D.I.K.R.T.A.V.V.P.S.A.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.V.S.A.A.C.L.L.N.N.F.Y.P.R.E.V.K.A.Q.W.K.A.D.N.V.L.Q.S.G.N.S.Q.E.S.A.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.V.D.Y.E.K.H.K.A.Y.V.C.E.A.T.H.Q.G.L.S.S.P.A.T.K.S.F.N.R.G.E.C}$PEPTIDE1,PEPTIDE1,22:R3-95:R3|PEPTIDE1,PEPTIDE1,143:R3-199:R3|PEPTIDE1,PEPTIDE3,219:R3-215:R3|PEPTIDE1,PEPTIDE2,225:R3-225:R3|PEPTIDE1,PEPTIDE2,228:R3-228:R3|PEPTIDE1,PEPTIDE1,260:R3-320:R3|PEPTIDE1,PEPTIDE1,366:R3-424:R3|PEPTIDE2,PEPTIDE2,22:R3-95:R3|PEPTIDE2,PEPTIDE2,143:R3-199:R3|PEPTIDE2,PEPTIDE4,219:R3-215:R3|PEPTIDE2,PEPTIDE2,260:R3-320:R3|PEPTIDE2,PEPTIDE2,366:R3-424:R3|PEPTIDE3,PEPTIDE3,23:R3-89:R3|PEPTIDE3,PEPTIDE3,135:R3-195:R3|PEPTIDE4,PEPTIDE4,23:R3-89:R3|PEPTIDE4,PEPTIDE4,135:R3-195:R3$$$";
			System.out.println("Edge String: " + ComplexNotationParser.getAllEdgeString(notation));

			//get all simple polymers
			Set<Monomer> set = new HashSet<Monomer>();
			List<PolymerNode> simplePolymers = ComplexNotationParser.getPolymerNodeList(notation);
			for (PolymerNode node: simplePolymers) {
				System.out.println(node.getId()  + ": " + node.getLabel() + "(" + node.getType() + ")");
				List<Monomer> monomers = SimpleNotationParser.getMonomerList(node.getLabel(), node.getType());
				for (Monomer subnode: monomers) {
					System.out.println(subnode.getName() +" | " + subnode.getAlternateId());
					
					set.add(subnode);
				}

			}
			//Distinct monomers
			System.out.println("\n Distinct monomers");
			for (Monomer distinctmonomer: set) {
				System.out.println(distinctmonomer.getName() +" | " + distinctmonomer.getAlternateId());
				XMLOutputter outp = new XMLOutputter();
				String s = outp.outputString(MonomerParser.getMonomerElement(distinctmonomer)); 
				System.out.println(s);
				
			}




		} catch (Exception ex) {
			Logger.getLogger(ComplexNotationSample.class.getName()).log(Level.SEVERE, null, ex);
		}

	}

}
