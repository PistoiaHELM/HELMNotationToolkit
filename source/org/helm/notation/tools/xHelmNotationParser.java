package org.helm.notation.tools;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.Monomer;
import org.jdom.JDOMException;
import org.jdom.Element;

public class xHelmNotationParser {	
	public static String extractComplexNotationString(Element rootElement, boolean validate) throws MonomerException, IOException, JDOMException, NotationException, StructureException  {
		//extract helm string
		Element helmNotationElement = rootElement.getChild("HelmNotation");
		
		String helmNotation = helmNotationElement.getText();
			
		
		//add monomers to externalMonomerDB
		Element monomerListElement = rootElement.getChild("Monomers");
		@SuppressWarnings("unchecked")
		List<Element> elementList = (List<Element>) monomerListElement
				.getChildren("Monomer");

		MonomerFactory factory = MonomerFactory.getInstance();
		
		for (Element monomerElement : elementList) {
			Monomer m = MonomerParser.getMonomer(monomerElement);
			factory.addExternalMonomer(m);
		}
		
		
		//validate
		if ( validate) {
			ComplexNotationParser.validateComplexNotation(helmNotation);
		}
	
		return helmNotation;
	}

}
