package org.helm.notation.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerStore;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;


/**
 * Class to export to XHELM XML format
 * 
 * @author maisel
 *
 */
public class xHelmNotationExporter {

	public static final String XHELM_ELEMENT 			= "Xhelm";
	public static final String MONOMER_LIST_ELEMENT		= "Monomers";
    public static final String MONOMER_ELEMENT 			= "Monomer";
  
    public static final String HELM_NOTATION_ELEMENT 	= "HelmNotation";
    

	/**
	 * Writes the XHELM XML format for the given notation and the corresponding monomer store
	 * 
	 * @param helmNotation
	 * @param store
	 * @return the XHELM format
	 */
	public static String writeXHELM(String helmNotation, MonomerStore store) {
		Element root = new Element(XHELM_ELEMENT);

		Document doc = new Document(root);

		Element helmElement = new Element(HELM_NOTATION_ELEMENT);
		helmElement.setText(helmNotation);

		root.addContent(helmElement);

		Element monomerListElement = new Element(MONOMER_LIST_ELEMENT);

		try {
			// get all simple polymers
			Set<Monomer> set = new HashSet<Monomer>();
			List<PolymerNode> simplePolymers = ComplexNotationParser
					.getPolymerNodeList(helmNotation);
			for (PolymerNode node : simplePolymers) {
				List<Monomer> monomers = SimpleNotationParser.getMonomerList(
						node.getLabel(), node.getType(), store);
				for (Monomer subnode : monomers) {
					if (!subnode.isAdHocMonomer()) {
						set.add(subnode);
					}
				}

			}
			// Distinct monomers
			for (Monomer distinctmonomer : set) {
				Element monomerElement = MonomerParser
						.getMonomerElement(distinctmonomer);
				monomerListElement.getChildren().add(monomerElement);

			}

		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}

		root.addContent(monomerListElement);

		XMLOutputter xmlOutput = new XMLOutputter();
		// display nice
		xmlOutput.setFormat(Format.getPrettyFormat());

		return xmlOutput.outputString(doc);

	}
	 
	 	
	
}
