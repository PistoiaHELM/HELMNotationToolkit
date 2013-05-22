/*******************************************************************************
 * Copyright C 2012, The Pistoia Alliance
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/
package org.helm.notation;

import org.helm.notation.tools.DeepCopy;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.MonomerCache;
import org.helm.notation.tools.MonomerParser;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.Attribute;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

/**
 * This is a factory class to build monomer database from MonomerDB.xml document
 * @author zhangtianhong
 */
public class MonomerDBHandler {
    
    public static final String POLYMER_LIST_ELEMENT = "POLYMER_LIST";
    public static final String POLYMER_ELEMENT = "POLYMER";
    public static final String POLYMER_TYPE_ATTRIBUTE = "polymerType";
    public static final String ATTACHMENT_LIST_ELEMENT = "ATTACHMENT_LIST";

    
    private static Logger logger = Logger.getLogger(MonomerDBHandler.class.toString());


    
    public static String buildMonomerDbXMLFromCache(MonomerCache cache) throws MonomerException {
        XMLOutputter outputer = new XMLOutputter(Format.getCompactFormat());

        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?><MONOMER_DB xmlns=\"lmr\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">");

        Map<String, Map<String, Monomer>> mDB = cache.getMonomerDB();
        Element polymerListElement = new Element(POLYMER_LIST_ELEMENT);
        Set<String> polymerTypeSet = mDB.keySet();
        for (Iterator i = polymerTypeSet.iterator();  i.hasNext();) {
            String polymerType = (String) i.next();
            Element polymerElement = new Element(POLYMER_ELEMENT);
            Attribute att = new Attribute(POLYMER_TYPE_ATTRIBUTE, polymerType);
            polymerElement.setAttribute(att);
            polymerListElement.getChildren().add(polymerElement);

            Map<String, Monomer> monomerMap = mDB.get(polymerType);
            Set<String> monomerSet = monomerMap.keySet();

            for(Iterator it = monomerSet.iterator(); it.hasNext();) {
                String monomerID = (String) it.next();
                Monomer m = monomerMap.get(monomerID);
                Element monomerElement = MonomerParser.getMonomerElement(m);
                polymerElement.getChildren().add(monomerElement);
            }
        }
        String polymerListString = outputer.outputString(polymerListElement);
        sb.append(polymerListString);


        Map<String, Attachment> aDB = cache.getAttachmentDB();
        Element attachmentListElement = new Element(ATTACHMENT_LIST_ELEMENT);
        Set<String> attachmentSet = aDB.keySet();
        for (Iterator itr = attachmentSet.iterator(); itr.hasNext();) {
            String attachmentID = (String) itr.next();
            Attachment attachment = aDB.get(attachmentID);
            Element attachmentElement = MonomerParser.getAttachementElement(attachment);
            attachmentListElement.getChildren().add(attachmentElement);
        }
        String attachmentListString = outputer.outputString(attachmentListElement);
        sb.append(attachmentListString);

        sb.append("</MONOMER_DB>");

        return sb.toString();
    }

}
