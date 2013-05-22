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
package org.helm.notation.tools;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import com.pfizer.pgrd.sdlib.EncoderException;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;

/**
 * This class provides mehtods that handle monomers used in polymer notation
 * @author zhangtianhong
 */
public class MonomerParser {

    public static final String MONOMER_ELEMENT = "MONOMER";
    public static final String MONOMER_ID_ELEMENT = "MONOMER_ID";
    public static final String MONOMER_SMILES_ELEMENT = "MONOMER_SMILES";
    public static final String MONOMER_MOL_FILE_ELEMENT = "MONOMER_MOL_FILE";
    public static final String MONOMER_TYPE_ELEMENT = "MONOMER_TYPE";
    public static final String POLYMER_TYPE_ELEMENT = "POLYMER_TYPE";
    public static final String NATURAL_ANALOG_ELEMENT = "NATURAL_ANALOG";
    public static final String MONOMER_NAME_ELEMENT = "MONOMER_NAME";
    public static final String ATTACHEMENTS_LIST_ELEMENT = "ATTACHMENT_LIST";
    public static final String ATTACHEMENTS_ELEMENT = "ATTACHMENTS";
    public static final String ATTACHEMENT_ELEMENT = "ATTACHMENT";
    public static final String ATTACHEMENT_ID_ELEMENT = "ATTACHMENT_ID";
    public static final String ATTACHEMENT_LABEL_ELEMENT = "ATTACHMENT_LABEL";
    public static final String CAP_GROUP_NAME_ELEMENT = "CAP_GROUP_NAME";
    public static final String CAP_GROUP_SMILES_ELEMENT = "CAP_GROUP_SMILES";
    private static List<String> polymerTypes = new ArrayList<String>();

    static {
        polymerTypes = Arrays.asList(Monomer.SUPPORTED_POLYMER_TYPES);
    }

    /**
     * Convert ATTACHMENT element to Attachment object
     * @param attachment element
     * @return Attachment
     */
    public static Attachment getAttachment(Element attachment) {

        Namespace ns = attachment.getNamespace();

        Attachment att = new Attachment();
        att.setAlternateId(attachment.getChildText(ATTACHEMENT_ID_ELEMENT, ns));
        att.setLabel(attachment.getChildText(ATTACHEMENT_LABEL_ELEMENT, ns));
        att.setCapGroupName(attachment.getChildText(CAP_GROUP_NAME_ELEMENT, ns));
        att.setCapGroupSMILES(attachment.getChildText(CAP_GROUP_SMILES_ELEMENT, ns));

        return att;
    }

    /**
     * This method converts Attachment to ATTACHMENT XML element
     * @param att -- Attachment
     * @return Element 
     */
    public static Element getAttachementElement(Attachment att) {


        Element attachment = new Element(ATTACHEMENT_ELEMENT);

        if (null != att.getAlternateId() && att.getAlternateId().length() >0) {
            Element e = new Element(ATTACHEMENT_ID_ELEMENT);
            e.setText(att.getAlternateId());
            attachment.getChildren().add(e);
        }

        if (null != att.getLabel() && att.getLabel().length() >0) {
            Element e = new Element(ATTACHEMENT_LABEL_ELEMENT);
            e.setText(att.getLabel());
            attachment.getChildren().add(e);
        }

        if (null != att.getCapGroupName() && att.getCapGroupName().length() >0) {
            Element e = new Element(CAP_GROUP_NAME_ELEMENT);
            e.setText(att.getCapGroupName());
            attachment.getChildren().add(e);
        }

        if (null != att.getCapGroupSMILES() && att.getCapGroupSMILES().length() >0) {
            Element e = new Element(CAP_GROUP_SMILES_ELEMENT);
            e.setText(att.getCapGroupSMILES());
            attachment.getChildren().add(e);
        }

        return attachment;
    }

    /**
     * This method validates Attachment by the following rules<br>
     * <li>Attachment must have unique ID<br>
     * <li>cap group SMILES must be valid<br>
     * <li>cap group SMILES must contain one R group<br>
     * <li>R group in SMILES must match R group label<br>
     * @param attachment
     * @return true or false
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     */
    public static boolean validateAttachement(Attachment attachment) throws MonomerException, IOException {

        String alternateId = attachment.getAlternateId();
        if (null == alternateId) {
            throw new MonomerException("Attachment must have unique ID");
        }

        String smiles = attachment.getCapGroupSMILES();
        if (null != smiles) {

            if (!StructureParser.validateSmiles(smiles)) {
                throw new MonomerException("Attachment cap group SMILES is invalid");
            }

            List<String> labels = getAttachmentLabels(smiles);
            if (null == labels || labels.size() != 1) {
                throw new MonomerException("Attachment must have one R group in SMILES");
            }

            if (!(labels.get(0).equals(attachment.getLabel()))) {
                throw new MonomerException("R group in monomer SMILES and R group label must match");
            }
        }
        return true;
    }

    /**
     * Convert monomer element to Monomer object
     * @param monomer element
     * @return Monomer
     */
    public static Monomer getMonomer(Element monomer) throws MonomerException {
        Monomer m = new Monomer();
        Namespace ns = monomer.getNamespace();
        m.setAlternateId(monomer.getChildText(MONOMER_ID_ELEMENT, ns));
        m.setCanSMILES(monomer.getChildText(MONOMER_SMILES_ELEMENT, ns));
        String encodedMolfile = monomer.getChildText(MONOMER_MOL_FILE_ELEMENT, ns);

        String molfile = null;
        try {
            molfile = MolfileEncoder.decode(encodedMolfile);
        } catch (EncoderException ex) {
            throw new MonomerException("Invalid monomer molfile");
        }
        m.setMolfile(molfile);
        m.setMonomerType(monomer.getChildText(MONOMER_TYPE_ELEMENT, ns));
        m.setPolymerType(monomer.getChildText(POLYMER_TYPE_ELEMENT, ns));
        m.setNaturalAnalog(monomer.getChildText(NATURAL_ANALOG_ELEMENT, ns));
        m.setName(monomer.getChildText(MONOMER_NAME_ELEMENT, ns));
        Element attachmentElement = monomer.getChild(ATTACHEMENTS_ELEMENT, ns);

        if (null != attachmentElement) {
            List attachments = attachmentElement.getChildren(ATTACHEMENT_ELEMENT, ns);
            List<Attachment> l = new ArrayList<Attachment>();
            Iterator i = attachments.iterator();
            while (i.hasNext()) {
                Element attachment = (Element) i.next();
                Attachment att = getAttachment(attachment);
                l.add(att);
            }
            m.setAttachmentList(l);
        }
        return m;
    }

    /**
     * This method converts Monomer to MONOMER XML element
     * @param monomer
     * @return Element 
     */
    public static Element getMonomerElement(Monomer monomer) throws MonomerException {
        Element element = new Element(MONOMER_ELEMENT);

        if (null != monomer.getAlternateId()) {
            Element e = new Element(MONOMER_ID_ELEMENT);
            e.setText(monomer.getAlternateId());
            element.getChildren().add(e);
        }

        if (null != monomer.getCanSMILES()) {
            Element e = new Element(MONOMER_SMILES_ELEMENT);
            e.setText(monomer.getCanSMILES());
            element.getChildren().add(e);
        }

        if (null != monomer.getMolfile()) {
            Element e = new Element(MONOMER_MOL_FILE_ELEMENT);
            String encodedMolfile = null;                  
            try {
                encodedMolfile = MolfileEncoder.encode(monomer.getMolfile());
            } catch (EncoderException ex) {
                throw new MonomerException("Invalid monomer molfile");
            }
//            CDATA cdata = new CDATA(monomer.getMolfile());
//            e.setContent(cdata);
            e.setText(encodedMolfile);
            element.getChildren().add(e);
        }

        if (null != monomer.getMonomerType()) {
            Element e = new Element(MONOMER_TYPE_ELEMENT);
            e.setText(monomer.getMonomerType());
            element.getChildren().add(e);
        }

        if (null != monomer.getPolymerType()) {
            Element e = new Element(POLYMER_TYPE_ELEMENT);
            e.setText(monomer.getPolymerType());
            element.getChildren().add(e);
        }

        if (null != monomer.getNaturalAnalog()) {
            Element e = new Element(NATURAL_ANALOG_ELEMENT);
            e.setText(monomer.getNaturalAnalog());
            element.getChildren().add(e);
        }

        if (null != monomer.getName()) {
            Element e = new Element(MONOMER_NAME_ELEMENT);
            e.setText(monomer.getName());
            element.getChildren().add(e);
        }

        List<Attachment> l = monomer.getAttachmentList();
        if (null != l && l.size() > 0) {
            Element attachments = new Element(ATTACHEMENTS_ELEMENT);

            for (int i = 0; i < l.size(); i++) {
                Attachment att = l.get(i);
                Element attachment = getAttachementElement(att);
                attachments.getChildren().add(attachment);
            }
            element.getChildren().add(attachments);
        }

        return element;
    }

    public static List<Monomer> getMonomerList(String monomerXMLString) throws JDOMException, IOException, MonomerException {
        List<Monomer> l = new ArrayList<Monomer>();
        if (null != monomerXMLString && monomerXMLString.length() > 0) {
            SAXBuilder builder = new SAXBuilder();
            ByteArrayInputStream bais = new ByteArrayInputStream(monomerXMLString.getBytes());              
            Document doc = builder.build(bais);
            Element root = doc.getRootElement();

            List monomers = root.getChildren();
            Iterator it = monomers.iterator();
            while (it.hasNext()) {
                Element monomer = (Element) it.next();
                Monomer m = getMonomer(monomer);
                if (MonomerParser.validateMonomer(m)) {
                    l.add(m);
                }
            }
        }
        return l;
    }
    
    public static Monomer getMonomer(String monomerXMLString) throws JDOMException, IOException, MonomerException {
        Monomer m = null;
        if (monomerXMLString != null && monomerXMLString.length() >0) {
        SAXBuilder builder = new SAXBuilder();
            ByteArrayInputStream bais = new ByteArrayInputStream(monomerXMLString.getBytes());              
            Document doc = builder.build(bais);
            Element root = doc.getRootElement();
            m = getMonomer(root); 
        } 
        return m;
            
    }

    /**
     * This methods checks the validity of the monomer based on the following rules<br>
     * <li>monomer cannot be null<br>
     * <li>polymer type cannot be null and must be one of the defined polymer type<br>
     * <li>monomer type cannot be null and must be one of the defined monomer type for a given polymer type<br>
     * <li>Monomer ID cannot be null<br>
     * <li>structure cannot be null for non-chemical type monomer<br>
     * <li>structure SMILES must be valid<br>
     * <li>attachment labels on monomer must be unique<br>
     * <li>Attachment number on SMILES must match attachment List size<br>
     * <li>Each attachment in attachment list must be valid (call validateAttachment())<br>
     * <li>Attachment labels on monomer must match atachment label on attachment list<br>
     * <li>For non-chemical type monomers, modified monomer (ID length greater than 1) must have natural analog<br>
     * <li>All monomers must have at least one attachment
     * @param monomer
     * @return true or false
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     */
    public static boolean validateMonomer(Monomer monomer) throws MonomerException, IOException {

        if (null == monomer) {
            throw new MonomerException("Monomer is null");
        } else {
            String polymerType = monomer.getPolymerType();
            if (null == polymerType) {
                throw new MonomerException("Monomer has no polymer type defined");
            } else if (!polymerTypes.contains(polymerType)) {
                throw new MonomerException("Unknown polymer type '" + polymerType + "'");
            }

            String monomerType = monomer.getMonomerType();
            if (null == monomerType) {
                throw new MonomerException("Monomer has no monomer type defined");
            } else {
                if (polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
                    if (!monomerType.equals(Monomer.UNDEFINED_MOMONER_TYPE)) {
                        throw new MonomerException("Valid monomer type for chemical structures can only be '" + Monomer.UNDEFINED_MOMONER_TYPE + "'");
                    }
                } else {
                    if (!(monomerType.equals(Monomer.BACKBONE_MOMONER_TYPE) || monomerType.equals(Monomer.BRANCH_MOMONER_TYPE))) {
                        throw new MonomerException("Valid monomer type for simple polymer can only be '" + Monomer.BACKBONE_MOMONER_TYPE + "' or '" + Monomer.BRANCH_MOMONER_TYPE + "'");
                    }
                }
            }

            String alternateId = monomer.getAlternateId();
            if (null == alternateId || alternateId.length() ==0) {
                throw new MonomerException("Monomer has no monomerID defined");
            }
            String smiles = monomer.getCanSMILES();
            String molfile = monomer.getMolfile();
            List<Attachment> attachments = monomer.getAttachmentList();

            if (!polymerType.equals(Monomer.CHEMICAL_POLYMER_TYPE)) {
                if (null == smiles || null == molfile || null == attachments || attachments.size() == 0) {
                    throw new MonomerException("Monomers for specific polymer type must have structure info");
                }
            }

            String errorNote = alternateId+" ("+polymerType+")";
            if (null != smiles && smiles.length() >0) {

                boolean validSmiles = StructureParser.validateSmiles(smiles);
                if (!validSmiles) {
                    throw new MonomerException("Monomer SMILES must be valid: "+ errorNote);
                }
                List<String> attachmentLabels = getAttachmentLabels(smiles);
                boolean unique = areAttachmentLabelsUnique(attachmentLabels);
                if (!unique) {
                    throw new MonomerException("Attachment labels on monomer must be unique: "+ errorNote);
                }
                if (attachmentLabels.size() != attachments.size()) {
                    throw new MonomerException("Attachment label number on monomer must match attachment number: "+errorNote);
                }
                for (int i = 0; i < attachments.size(); i++) {
                    Attachment att = attachments.get(i);
                    validateAttachement(att);
                }

                for (int i = 0; i < attachmentLabels.size(); i++) {
                    String label = attachmentLabels.get(i);
                    boolean found = false;
                    for (int j = 0; j < attachments.size(); j++) {
                        Attachment att = attachments.get(j);
                        if (att.getAlternateId().startsWith(label)) {
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        throw new MonomerException("Attachment label in SMILES is not found in attachment list: "+errorNote);
                    }
                }
            }

            if (monomer.getAlternateId().length() > 0 && !(monomer.getPolymerType().equals(Monomer.CHEMICAL_POLYMER_TYPE))) {
                String naturalAnalog = monomer.getNaturalAnalog();

                if (null == naturalAnalog) {
                    throw new MonomerException("Modified monomer must have natural analog defined: "+errorNote);
                } else {
                    if (naturalAnalog.length() != 1) {
                        throw new MonomerException("Natural analog must be single letter: "+errorNote);
                    }
                }
            }

            if (monomer.getAttachmentList() == null || monomer.getAttachmentList().size() == 0) {
                throw new MonomerException("Monomer must have at least one attachment: "+errorNote);
            }

            //make sure R group can only be connected to one atom via single achiral bond
            //MolBond javadoc: getType()Gets the bond type. Possible values: 1 (single), 2 (double), 3 (triple), coordinate, conjugated and query bond types.
            if (null != smiles && smiles.length() >0) {
                Molecule molecule = StructureParser.getMolecule(smiles);
                List<String> attachmentLabels = getAttachmentLabels(smiles);
                for (int i = 0; i < attachmentLabels.size(); i++) {
                    String rgroupId = attachmentLabels.get(i).substring(1);
                    MolAtom atom =null;
                    try {
                        atom = StructureParser.getRgroupAtom(molecule, Integer.parseInt(rgroupId));
                    } catch (StructureException ex) {
                        throw new MonomerException("Unable to find R"+rgroupId+ " in monomer: "+ errorNote);
                    }
                    if (atom.getBondCount() != 1) {
                        throw new MonomerException("R group can only connect with one atom in monomer: "+errorNote);
                    } else {
                        MolBond bond = atom.getBond(0);
                        if (bond.getType() != 1)
                            throw new MonomerException("R group can only connect with another atom via single bond in monomer: "+ errorNote);
                    }
                }
            }
        }

        return true;
    }

    /**
     * Convert the extendedSMILES of the monomer mixture (monomer and default capping groups at attachment points)
     * to a Monomer object. Only structure info for monomer is instantiated.
     * Input: [*][H].O[*].OP([*])([*])=O |$_R2;;;_R1;;;_R1;_R2;$|
     * Output:
     *   OP([*])([*])=O |$;;_R1;_R2;$|
     *   O[*] |$;_R1$|
     *   [*][H] |$_R2;$|
     * @param mixtureExtendedSMILES the extended smiles for the mixture, which include the monomer, and attachment structure
     * @return an instance of Monomer with structure info, still need to fill in non-structure info, such as polymer type, monomer type, alternate ID
     * @throws org.helm.notation.MonomerException
     */
    private static Monomer decomposeMonomerStructure(String mixtureExtendedSMILES) throws StructureException, MonomerException, IOException {

        boolean valid = isValidMonomerSMILES(mixtureExtendedSMILES);
        List<String> smilesExtList = StructureParser.decomposeMixture(mixtureExtendedSMILES);

        Monomer monomer = new Monomer();
        List<Attachment> attachments = new ArrayList<Attachment>();

        //single attachment monomer
        if (smilesExtList.size() == 2) {
            String component1 = smilesExtList.get(0);
            String component2 = smilesExtList.get(1);


            if (component1.length() >= component2.length()) {
                monomer.setCanSMILES(component1);
                Attachment att = new Attachment();
                att.setCapGroupSMILES(component2);
                String label = getAttachmentLabel(component2);
                att.setLabel(label);
                attachments.add(att);
            } else {
                monomer.setCanSMILES(component2);
                Attachment att = new Attachment();
                att.setCapGroupSMILES(component1);
                String label = getAttachmentLabel(component1);
                att.setLabel(label);
                attachments.add(att);
            }

        } else {
            for (int i = 0; i < smilesExtList.size(); i++) {
                String component = smilesExtList.get(i);

                String[] tokens = component.split("\\*");

                //monomer or attachment
                if (tokens.length > 2) {
                    monomer.setCanSMILES(component);
                } else {
                    Attachment att = new Attachment();
                    att.setCapGroupSMILES(component);
                    String label = getAttachmentLabel(component);
                    att.setLabel(label);
                    attachments.add(att);
                }
            }
        }

        monomer.setAttachmentList(attachments);

        return monomer;
    }

    /**
     * This methods checks if the mixture extended SMILES is valid to represent a monomer and its attachments structure
     * @param mixtureExtendedSMILES
     * @return true or false
     * @throws org.helm.notation.StructureException
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     */
    private static boolean isValidMonomerSMILES(String mixtureExtendedSMILES) throws StructureException, MonomerException, IOException {
        List<String> componentSmiles = StructureParser.decomposeMixture(mixtureExtendedSMILES);

        List<String> singles = new ArrayList<String>();
        List<String> multiples = new ArrayList<String>();

        for (int i = 0; i < componentSmiles.size(); i++) {
            int asteriskNumber = componentSmiles.get(i).split("\\*").length - 1;
            if (asteriskNumber > 1) {
                multiples.add(componentSmiles.get(i));
            } else if (asteriskNumber == 1) {
                singles.add(componentSmiles.get(i));
            } else {
                throw new StructureException("Invalid Monomer Structure: component contains no R group");
            }
        }

        if (multiples.size() == 1) {

            int count = multiples.get(0).split("\\*").length - 1;

            if (count != singles.size()) {
                throw new StructureException("Invalid Monomer Structure: number of attachment points on monomer does not match number of attachments");
            }

            List<String> multipleLabels = getAttachmentLabels(multiples.get(0));
            Map<String, Integer> multipleMap = new HashMap<String, Integer>();
            for (int i = 0; i < multipleLabels.size(); i++) {
                multipleMap.put(multipleLabels.get(i), new Integer(i));
            }
            if (multipleLabels.size() != multipleMap.size()) {
                throw new StructureException("Invalid Monomer Structure: R groups on monomer are not unique");
            }
            List<String> singleLabels = new ArrayList<String>();
            for (int i = 0; i < singles.size(); i++) {
                String label = getAttachmentLabel(singles.get(i));
                singleLabels.add(label);
            }
            Map<String, Integer> singleMap = new HashMap<String, Integer>();
            for (int i = 0; i < singleLabels.size(); i++) {
                singleMap.put(singleLabels.get(i), new Integer(i));
            }
            if (singleLabels.size() != singleMap.size()) {
                throw new StructureException("Invalid Monomer Structure: R groups on caping components are not unique");
            }
            Set<String> singleKeySet = singleMap.keySet();
            for (Iterator it = singleKeySet.iterator(); it.hasNext();) {
                String key = (String) it.next();
                if (!(multipleMap.containsKey(key))) {
                    throw new MonomerException("Invalid Monomer Structure: R groups on monomer and R groups on caping components do not match");
                }
            }

        } else if (multiples.size() > 1) {
            throw new MonomerException("Invalid Monomer Structure: more than one component contains multiple R groups");
        } else {
            //we can only have two singles
            if (singles.size() != 2) {
                throw new MonomerException("Invalid Monomer Structure: number of attachment points on monomer does not match number of attachments");
            } else {
                String label1 = getAttachmentLabel(singles.get(0));
                String label2 = getAttachmentLabel(singles.get(1));

                if (!(label1.equals(label2))) {
                    throw new MonomerException("Invalid Monomer Structure: attachment label on monomer does not match attachment label on capping group");
                }
            }
        }

        return true;
    }

    /**
     * This method returns the first R group label in the SMILES string
     * @param smilesExtension
     * @return rgroup label
     * @throws org.helm.notation.MonomerException
     */
    private static String getAttachmentLabel(String smilesExtension) throws MonomerException {
        int rPos = smilesExtension.indexOf("R");
        StringBuffer sb = new StringBuffer();
        while (rPos > 0) {
            rPos++;
            String nextLetter = smilesExtension.substring(rPos, rPos + 1);
            if (nextLetter.matches("[0-9]")) {
                sb.append(nextLetter);
            } else {
                break;
            }
        }
        if (sb.length() > 0) {
            return "R" + sb.toString();
        } else {
            throw new MonomerException("Invalid Monomer Structure: there is no attachment label");
        }
    }

    /**
     * This methods return the list of R groups in the extended SMILES string
     * @param extendedSmiles
     * @return string list
     */
    private static List<String> getAttachmentLabels(String extendedSmiles) {
        List<String> labels = new ArrayList<String>();

        int start = 0;
        int rPos = extendedSmiles.indexOf("R");
        StringBuffer sb = new StringBuffer();
        while (rPos > 0) {
            rPos++;
            String nextLetter = extendedSmiles.substring(rPos, rPos + 1);
            if (nextLetter.matches("[0-9]")) {
                sb.append(nextLetter);
            } else {
                labels.add("R" + sb.toString());
                sb = new StringBuffer();
                start = rPos + 1;
                rPos = extendedSmiles.indexOf("R", start);
            }
        }
        return labels;
    }

    /**
     * This mehtod checks if strings in a list are unique
     * @param labels
     * @return true or fals
     */
    private static boolean areAttachmentLabelsUnique(List<String> labels) {
        Map<String, String> map = new HashMap<String, String>();
        for (int i = 0; i < labels.size(); i++) {
            map.put(labels.get(i), labels.get(i));
        }
        if (labels.size() == map.size()) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * This method checks if attachment label is in the format of R#, where # is a number
     * @param label
     * @throws org.helm.notation.MonomerException
     */
    public static void validateAttachmentLabel(String label) throws MonomerException {

        if (label.equalsIgnoreCase(Attachment.PAIR_ATTACHMENT)) {
            return;
        }
        char[] chars = label.toCharArray();

        if (!(String.valueOf(chars[0])).equals("R")) {
            throw new MonomerException("Invalid Attachment Label format");
        }
        for (int i = 1; i < chars.length; i++) {
            char c = chars[i];
            if (!(String.valueOf(c)).matches("[0-9]")) {
                throw new MonomerException("Invalid Attachment Label format");
            }
        }
    }

    public static void fillAttachmentInfo(Attachment att) throws MonomerException, IOException, JDOMException {
        Map<String, Attachment> attachmentMap = MonomerFactory.getInstance().getAttachmentDB();

        Attachment attach = attachmentMap.get(att.getAlternateId());
        att.setLabel(attach.getLabel());
        att.setCapGroupSMILES(attach.getCapGroupSMILES());
        att.setCapGroupName(attach.getCapGroupName());
    }
}
