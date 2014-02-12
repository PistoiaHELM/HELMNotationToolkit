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

	public static final String PEPTIDE_POLYMER_TYPE = "PEPTIDE";
	public static final String CHEM_POLYMER_TYPE = "CHEM";
	public static final String RNA_POLYMER_TYPE = "RNA";

	
	
	
	public static String extractComplexNotationString(Element rootElement, boolean loadExternalMonomerDB, boolean validate) throws MonomerException, IOException, JDOMException, NotationException, StructureException  {
		
			
		Element helmNotationElement = rootElement.getChild("HelmNotation");
		
		if ( loadExternalMonomerDB) {
			loadExternalMonomerDB( rootElement);
		}

		String helmNotation = helmNotationElement.getText();
		
		
		if ( validate) {
			ComplexNotationParser.validateComplexNotation(helmNotation);
		}
	
		return helmNotation;
	}

	private static void loadExternalMonomerDB( Element rootElement) throws MonomerException, IOException, JDOMException {
		MonomerFactory monomerFactory = MonomerFactory.getInstance();
		Element monomerListElement = rootElement.getChild("Monomers");
		@SuppressWarnings("unchecked")
		List<Element> elementList = (List<Element>) monomerListElement
				.getChildren("Monomer");
		Map<String, Map<String, Monomer>> monomerDB = xHelmNotationParser
				.buildMonomerDB(elementList);
		
		monomerFactory.setExternalMonomerDB(monomerDB);

	}
	
	private static Map<String, Map<String, Monomer>> buildMonomerDB(
			List<Element> monomers) throws MonomerException, IOException {
		Map<String, Map<String, Monomer>> monomerDB = new HashMap<String, Map<String, Monomer>>();
		
		Map<String, Monomer> pepMonomerMap = new HashMap<String, Monomer>();
		Map<String, Monomer> rnaMonomerMap = new HashMap<String, Monomer>();
		Map<String, Monomer> chemMonomerMap = new HashMap<String, Monomer>();
		
		
		for (Element monomerElement : monomers) {

			Monomer m = MonomerParser.getMonomer(monomerElement);
			if (MonomerParser.validateMonomer(m)) {
				String polymerType = m.getPolymerType();
				if (polymerType.equalsIgnoreCase(PEPTIDE_POLYMER_TYPE)) {
					pepMonomerMap.put(m.getAlternateId(), m);
				} else if (polymerType.equalsIgnoreCase(RNA_POLYMER_TYPE)) {
					rnaMonomerMap.put(m.getAlternateId(), m);
				} else if (polymerType.equalsIgnoreCase(CHEM_POLYMER_TYPE)) {
					chemMonomerMap.put(m.getAlternateId(), m);
				}

			}
		}
		monomerDB.put(PEPTIDE_POLYMER_TYPE, pepMonomerMap);
		monomerDB.put(RNA_POLYMER_TYPE, rnaMonomerMap);
		monomerDB.put(CHEM_POLYMER_TYPE, chemMonomerMap);

		return monomerDB;
	}


	
}
