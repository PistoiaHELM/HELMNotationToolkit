package org.helm.notation.tools;

import static org.junit.Assert.*;

import java.io.FileInputStream;
import java.io.IOException;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

public class xHelmNotationParserTest {
	

	@Test
	public void testParseXHelmNotation() throws JDOMException, IOException, MonomerException,
			NotationException, StructureException {

		FileInputStream in = new FileInputStream(
				"samples/PeptideLinkerNucleotide.xhelm");
		SAXBuilder builder = new SAXBuilder();
		Document doc = builder.build(in);
		Element xHELMElement = doc.getRootElement();
			String helmString = xHelmNotationParser.extractComplexNotationString(
				xHELMElement, true, true);
		
		assertEquals("RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}|PEPTIDE1{[aaa].C.G.K.E.D.K.R}|CHEM1{SMCC}$PEPTIDE1,CHEM1,2:R3-1:R2|RNA1,CHEM1,1:R1-1:R1$$$",helmString);
		
	}
	
	@BeforeClass
	public static void init() {
		try {
			MonomerFactory.finalizeMonomerCache();
			MonomerFactory.getInstance();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@AfterClass
	public static void finish() {
		MonomerFactory.finalizeMonomerCache();
		
	}

}
