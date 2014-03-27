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

public class xHelmNotationParser {
	public static String getComplexNotationString(Element rootElement)
	{
		Element helmNotationElement = rootElement.getChild("HelmNotation");
		return helmNotationElement.getText();
	}

	// read monomers to monomerStore
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

	public static String writeXHELM(String xHelmCode, MonomerStore store) {
		String xml = "";
		try {

			xml += "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
					+ System.getProperty("line.separator");;
			xml += "<Xhelm>" + System.getProperty("line.separator");;
			xml += "	<HelmNotation>" + xHelmCode + "</HelmNotation>"
					+ System.getProperty("line.separator");
			xml += "<Monomers>" + System.getProperty("line.separator");

			// get all simple polymers
			Set<Monomer> set = new HashSet<Monomer>();
			List<PolymerNode> simplePolymers = ComplexNotationParser
					.getPolymerNodeList(xHelmCode);
			for (PolymerNode node : simplePolymers) {
				List<Monomer> monomers = SimpleNotationParser.getMonomerList(
						node.getLabel(), node.getType(), store);
				for (Monomer subnode : monomers) {
					set.add(subnode);
				}

			}
			// Distinct monomers
			for (Monomer distinctmonomer : set) {
				XMLOutputter outp = new XMLOutputter();
				String s = outp.outputString(MonomerParser
						.getMonomerElement(distinctmonomer));
				xml += "	" + s + System.getProperty("line.separator")
						+ System.getProperty("line.separator");

			}
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}

		xml += "</Monomers>" + System.getProperty("line.separator");
		xml += "</Xhelm>" + System.getProperty("line.separator");
		return xml;
	}

}
