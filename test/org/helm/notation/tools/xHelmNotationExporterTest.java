package org.helm.notation.tools;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerStore;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;
import org.junit.Test;

import chemaxon.marvin.plugin.PluginException;

public class xHelmNotationExporterTest {

	@Test
	public void testExport() throws JDOMException, IOException, MonomerException {
		
		//import
		InputStream in = this.getClass().getResourceAsStream("resources/PeptideLinkerNucleotide.xhelm");
		SAXBuilder builder = new SAXBuilder();
		Document doc = builder.build(in);
		
		
		XMLOutputter outputter=new XMLOutputter(Format.getPrettyFormat()) ;
		String expectedDocString=outputter.outputString(doc);
		
		String helmString = xHelmNotationParser.getComplexNotationString(doc.getRootElement());
		MonomerStore monomerStore = xHelmNotationParser.getMonomerStore(doc.getRootElement());	
		
		//export
		String output=xHelmNotationExporter.writeXHELM(helmString, monomerStore);
	
				
		assertEquals(expectedDocString.length(), output.length());
		
	
		
		
	}

}
