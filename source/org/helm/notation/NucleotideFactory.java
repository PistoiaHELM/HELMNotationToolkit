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

import org.helm.notation.tools.NucleotideSequenceParser;
import org.helm.notation.tools.SimpleNotationParser;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

/**
 * This is a factory class to build nucleotide templates
 * @author zhangtianhong
 */
public class NucleotideFactory {

    public static final String NOTATION_DIRECTORY = NotationConstant.NOTATION_DIRECTORY;
    public static final String LOCAL_NUCLEOTIDE_TEMPLATE_FILE_NAME = "NucleotideTemplates.xml";
    public static final String LOCAL_NUCLEOTIDE_TEMPLATE_FILE_PATH = NOTATION_DIRECTORY + System.getProperty("file.separator") + LOCAL_NUCLEOTIDE_TEMPLATE_FILE_NAME;
    public static final String NUCLEOTIDE_TEMPLATE_XML_RESOURCE = "resources/NucleotideTemplates.xml";
    public static final String NUCLEOTIDE_TEMPLATE_SCHEMA_RESOURCE = "resources/NucleotideTemplateSchema.xsd";
    public static final String XML_SCHEMA_VALIDATION_FEATURE = "http://apache.org/xml/features/validation/schema";
    public static final String EXTERNAL_SCHEMA_LOCATION_KEY = "http://apache.org/xml/properties/schema/external-schemaLocation";
    public static final String DEFAULT_NAME_SPACE = "lmr";
    private static NucleotideFactory instance;
    /**
     * First key is notation source, such as "HELM Notation"
     * Second key is nucleotide symbol, such as "A"
     */
    private static Map<String, Map<String, String>> nucleotideTemplates;
    private static Map<String, String> reverseNucleotideMap;
    private static SAXBuilder builder;
    private static Logger logger = Logger.getLogger(NucleotideFactory.class.toString());

    public synchronized Map<String, Map<String, String>> getNucleotideTemplates() {
        return nucleotideTemplates;
    }

    public synchronized void setNucleotideTemplates(Map<String, Map<String, String>> newNucleotideTemplates) {
        nucleotideTemplates = newNucleotideTemplates;
        reverseNucleotideMap = getReverseNucleotideTemplateMap(NotationConstant.NOTATION_SOURCE);
    }

    private static void setupBuilder() {
        URL schema = MonomerFactory.class.getResource(NUCLEOTIDE_TEMPLATE_SCHEMA_RESOURCE);
        builder = new SAXBuilder(true); //checks both well-formedness and validity
        builder.setFeature(XML_SCHEMA_VALIDATION_FEATURE, true);
        builder.setProperty(EXTERNAL_SCHEMA_LOCATION_KEY, DEFAULT_NAME_SPACE + " " + schema.toString());
    }

    private NucleotideFactory() {
    }

    public static NucleotideFactory getInstance() throws IOException, JDOMException, NotationException {
        if (null == instance) {
            initializeNucleotideTemplates();
            instance = new NucleotideFactory();
        }
        return instance;
    }

    private static Map<String, Map<String, String>> buildNucleotideTemplates(InputStream templatesInputStream) throws IOException, JDOMException {
        if (null == builder) {
            setupBuilder();
        }
        Document doc = builder.build(templatesInputStream);
        Element root = doc.getRootElement();
        Map<String, Map<String, String>> templates = NucleotideSequenceParser.getNucleotideTemplates(root);
        return templates;
    }

    public Map<String, Map<String, String>> buildNucleotideTemplatesFromXML(String nucleotideTemplatesXML) throws IOException, JDOMException {
        ByteArrayInputStream bais = new ByteArrayInputStream(nucleotideTemplatesXML.getBytes());
        return buildNucleotideTemplates(bais);
    }

    /**
     * This method is called during startup, use local version if exist, otherwise use XML version in jar
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     * @throws org.jdom.JDOMException
     */
    private static void initializeNucleotideTemplates() throws IOException, JDOMException, NotationException {

        InputStream in = null;
        File localFile = new File(LOCAL_NUCLEOTIDE_TEMPLATE_FILE_PATH);
        Map<String, Map<String, String>> templates = null;

        if (localFile.exists()) {
            try {
                in = new FileInputStream(localFile);
                templates = buildNucleotideTemplates(in);
                validate(templates);
                logger.log(Level.INFO, LOCAL_NUCLEOTIDE_TEMPLATE_FILE_PATH + " is used for nucleotide templates initialization");
            } catch (Exception e) {
                logger.log(Level.INFO, "Unable to use local nucleotide templates for initialization");
                localFile.delete();
                logger.log(Level.INFO, "Deleted local nucleotide templates file");
            }
        }

        if (null == templates) {
            in = NucleotideFactory.class.getResourceAsStream(NUCLEOTIDE_TEMPLATE_XML_RESOURCE);
            templates = buildNucleotideTemplates(in);
            validate(templates);
            logger.log(Level.INFO, NUCLEOTIDE_TEMPLATE_XML_RESOURCE + " is used for nucleotide templates initialization");
        }

        nucleotideTemplates = templates;
        reverseNucleotideMap = getReverseNucleotideTemplateMap(NotationConstant.NOTATION_SOURCE);
    }

    /**
     * save Nucleotide Templates to disk file
     * @throws java.io.IOException
     */
    public void saveNucleotideTemplates() throws IOException {
        File f = new File(NOTATION_DIRECTORY);
        if (!f.exists()) {
            f.mkdir();
        }
        String nucleotideTemplatesXML = NucleotideSequenceParser.getNucleotideTemplatesXML(getNucleotideTemplates());
        FileOutputStream fos = new FileOutputStream(LOCAL_NUCLEOTIDE_TEMPLATE_FILE_PATH);
        fos.write(nucleotideTemplatesXML.getBytes());
    }

    /**
     * key is notation such as R(A)P, value is symblo such as A
     * @param notationSource
     * @return notation/nucleotidSymbol map
     */
    private static synchronized Map<String, String> getReverseNucleotideTemplateMap(String notationSource) {
        Map<String, String> map = new HashMap<String, String>();
        Map<String, String> normalMap = nucleotideTemplates.get(notationSource);
        if (normalMap != null) {
            Set<String> set = normalMap.keySet();
            for (Iterator i = set.iterator(); i.hasNext();) {
                String tmp = (String) i.next();
                map.put((String) normalMap.get(tmp), tmp);
            }
        }
        return map;
    }

    /**
     * 
     * @return  reverse nucleotide template map for 'HELM Notation'
     */
    public synchronized Map<String, String> getReverseNucleotideTemplateMap() {
        return reverseNucleotideMap;
    }

    private static boolean validate(Map<String, Map<String, String>> templates) throws IOException, NotationException, JDOMException {
        Map<String, String> nucMap = templates.get(NotationConstant.NOTATION_SOURCE);
        Set<String> symbols = nucMap.keySet();
        for (String symbol : symbols) {
            String notation = nucMap.get(symbol);
            try {
                SimpleNotationParser.validateSimpleNotationForRNA(notation);
            } catch (Exception e) {
                throw new NotationException("Invalid notation " + notation + " for nucleotide " + symbol);
            }
        }
        return true;
    }
}
