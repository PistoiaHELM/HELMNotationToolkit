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

	private FileInputStream in;
	private SAXBuilder builder;
	private Document doc;
	private Element xHELMElement;
	private String helmString;

	@Test
	public void testParseXHelmNotation() throws JDOMException, IOException,
			MonomerException, NotationException, StructureException {

		in = new FileInputStream("samples/PeptideLinkerNucleotide.xhelm");
		builder = new SAXBuilder();
		doc = builder.build(in);
		xHELMElement = doc.getRootElement();
		String helmString = xHelmNotationParser.extractComplexNotationString(
				xHELMElement, true);

		assertEquals(
				"RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}|PEPTIDE1{[aaa].C.G.K.E.D.K.R}|CHEM1{SMCC}$PEPTIDE1,CHEM1,2:R3-1:R2|RNA1,CHEM1,1:R1-1:R1$$$",
				helmString);

	}
	
	@Test
	public void testXHelmValidation() throws JDOMException, IOException,
			MonomerException, NotationException, StructureException {

		boolean exception=false;
		try {
			in = new FileInputStream("samples/bad.xhelm");
			builder = new SAXBuilder();
			doc = builder.build(in);
			xHELMElement = doc.getRootElement();
			helmString = xHelmNotationParser.extractComplexNotationString(
					xHELMElement, true);			
		}

		catch (Exception e) {
			exception=true;
		}		
		assertTrue("Validation should have thrown an error",exception);

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
