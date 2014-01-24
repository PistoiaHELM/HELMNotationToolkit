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
