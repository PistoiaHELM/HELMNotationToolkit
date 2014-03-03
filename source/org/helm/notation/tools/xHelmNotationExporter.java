package org.helm.notation.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerStore;
import org.helm.notation.model.Monomer;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

public class xHelmNotationExporter {

	public static final String XHELM_ELEMENT 			= "xHelm";
	public static final String MONOMER_LIST_ELEMENT		= "Monomers";
    public static final String MONOMER_ELEMENT 			= "Monomer";
  
    public static final String HELM_NOTATION_ELEMENT 	= "HelmNotation";
    
     
	// URL: http://www.mkyong.com/java/how-to-create-xml-file-in-java-jdom-parser/
	 public static Document buildXHelmDocument(String helmNotation,MonomerStore monomerStore) throws MonomerException {
	        Element root = new Element(XHELM_ELEMENT);
			
	        Document doc = new Document(root);
			doc.setRootElement(root);
	        
	        Element helmElement=new Element(HELM_NOTATION_ELEMENT);
	        helmElement.setText(helmNotation);	        
	        
	        root.addContent(helmElement);
	       
	        Element monomerListElement = new Element(MONOMER_LIST_ELEMENT);
	        Set<String> polymerTypeSet = monomerStore.getPolymerTypeSet();
	        for (Iterator<String> i = polymerTypeSet.iterator();  i.hasNext();) {
	            String polymerType = i.next();
	  
	            Map<String, Monomer> monomerMap = monomerStore.getMonomers(polymerType);
	                 
	            Collection<Monomer> monomers=monomerMap.values();
	            for (Iterator<Monomer> iterator = monomers.iterator(); iterator
						.hasNext();) {
					Monomer monomer = iterator.next();
					Element monomerElement = MonomerParser.getMonomerElement(monomer);
					monomerListElement.getChildren().add(monomerElement);
				}	              
	        }
	     
	        root.addContent(monomerListElement);

	        return doc;
	    }
	 
	 	public static void saveXHelmDocument( Document document, File file) throws IOException {
	 		XMLOutputter xmlOutput = new XMLOutputter();	 
			// display nice 
			xmlOutput.setFormat(Format.getPrettyFormat());
		
	 		FileOutputStream fos = new FileOutputStream(file);
	 		xmlOutput.output(document, fos);
	 	}
	 
	 	public static void saveXHelm(String helmNotation,MonomerStore monomerStore,File file)  throws IOException, MonomerException {
	 		Document xHelm = buildXHelmDocument(helmNotation,monomerStore);
	        saveXHelmDocument( xHelm, file);  
	    }
	
}
