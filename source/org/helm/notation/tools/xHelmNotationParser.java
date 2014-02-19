package org.helm.notation.tools;

import java.io.IOException;
import java.util.List;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerStore;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.ComplexPolymer;
import org.helm.notation.model.Monomer;
import org.jdom.JDOMException;
import org.jdom.Element;

public class xHelmNotationParser {	
	public static String getComplexNotationString(Element rootElement) throws MonomerException, IOException, JDOMException, NotationException, StructureException  {
		Element helmNotationElement = rootElement.getChild("HelmNotation");
		return  helmNotationElement.getText();
	}

	//read monomers to monomerStore
	public static MonomerStore getMonomerStore(Element rootElement)
			throws MonomerException, IOException, JDOMException {
		MonomerStore monomerStore = new MonomerStore();
		Element monomerListElement = rootElement.getChild("Monomers");
		@SuppressWarnings("unchecked")
		List<Element> elementList = (List<Element>) monomerListElement
				.getChildren("Monomer");

		for (Element monomerElement : elementList) {
			Monomer m = MonomerParser.getMonomer(monomerElement);
			monomerStore.addMonomer(m);
		}
		return monomerStore;
	}
	
	

}
