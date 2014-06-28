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
import org.helm.notation.model.Monomer;
import java.io.InputStream;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class MonomerParserSample {

	/**
	 * @param args
	 *            the command line arguments
	 */
	public static void main(String[] args) {

		try {
			// monomer validation
			InputStream in = MonomerParserSample.class
					.getResourceAsStream("resource/BadMonomerSample.xml");
			SAXBuilder builder = new SAXBuilder();
			Document doc = builder.build(in);
			Element monomerElement = doc.getRootElement();
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());
			System.out.println("XML read in:");
			out.output(monomerElement, System.out);

			Monomer monomer = MonomerParser.getMonomer(monomerElement);
			if (MonomerParser.validateMonomer(monomer)) {
				System.out.println("Monomer is valid");
			} else {
				System.out.println("Monomer is invalid");
			}
			System.out.println("Molfile in monomer:");
			System.out.println(monomer.getMolfile());

			// Monomer to Element
			System.out.println("Convert Monomer to Element:");
			Element mElement = MonomerParser.getMonomerElement(monomer);
			out.output(mElement, System.out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
