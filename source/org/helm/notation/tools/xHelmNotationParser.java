package org.helm.notation.tools;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerStore;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.ComplexPolymer;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;
import org.jdom.JDOMException;
import org.jdom.Element;
import org.jdom.output.XMLOutputter;

/**
 * Class to parse the XHELM XML format
 * 
 * @author maisel
 *
 */
public class xHelmNotationParser {
	
	/**
	 * Extracts the complex notation string from the root node of the XHELM document
	 * 
	 * @param rootElement
	 * @return the complex notation string
	 */
	public static String getComplexNotationString(Element rootElement)
	{
		Element helmNotationElement = rootElement.getChild("HelmNotation");
		return helmNotationElement.getText();
	}

	/**
	 * Generates the monomer store from a given XHELM document
	 * 
	 * @param rootElement
	 * @return a monomer store
	 * @throws MonomerException
	 * @throws IOException
	 */
	public static MonomerStore getMonomerStore(Element rootElement)
			throws MonomerException, IOException {
		MonomerStore monomerStore = new MonomerStore();
		Element monomerListElement = rootElement.getChild("Monomers");
		if ( monomerListElement != null) {
			@SuppressWarnings("unchecked")
			List<Element> elementList = (List<Element>) monomerListElement
					.getChildren("Monomer");
	
			for (Element monomerElement : elementList) {
				Monomer m = MonomerParser.getMonomer(monomerElement);
				monomerStore.addMonomer(m);
			}
		}
		return monomerStore;
	}



}
