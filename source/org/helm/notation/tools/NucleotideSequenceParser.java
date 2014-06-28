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

import org.helm.notation.NotationConstant;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.model.Nucleotide;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jdom.Attribute;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

/**
 * This class provides methods that handle nucleotide sequences
 * 
 * @author ZHANGTIANHONG
 */
public class NucleotideSequenceParser {

	@Deprecated
	public static final String DEFAULT_NOTATION_SOURCE = "HELM Notation";
	public static final int MINUMUM_MATCH_FRAGMENT_LENGTH = 2;
	public static Map<String, String> complementMap = new HashMap<String, String>();

	static {
		complementMap.put("A", "U");
		complementMap.put("G", "C");
		complementMap.put("C", "G");
		complementMap.put("U", "A");
		complementMap.put("T", "A");
		complementMap.put("X", "X");
	}
	public static final String RNA_DESIGN_NONE = "NONE";
	public static final String RNA_DESIGN_TUSCHL_19_PLUS_2 = "TUSCHL_19_PLUS_2";
	// ss 5 1------19--
	// as 3 --19------1
	public static final String RNA_DESIGN_DICER_27_R = "DICER_27_R";
	// ss 5' 1-----------------------25
	// as 3' 123-----------------------27
	public static final String RNA_DESIGN_DICER_27_L = "DICER_27_L";
	// ss 5' 1-------------------------27
	// as 3' 1-----------------------25
	public static final List<String> SUPPORTED_DESIGN_LIST = new ArrayList<String>();

	static {
		SUPPORTED_DESIGN_LIST.add(RNA_DESIGN_NONE);
		SUPPORTED_DESIGN_LIST.add(RNA_DESIGN_TUSCHL_19_PLUS_2);
		SUPPORTED_DESIGN_LIST.add(RNA_DESIGN_DICER_27_L);
		SUPPORTED_DESIGN_LIST.add(RNA_DESIGN_DICER_27_R);
	}
	private static final String NUCLEOTIDE_SYMBOL_ELEMENT = "SYMBOL";
	private static final String NUCLEOTIDE_MONOMER_NOTATION_ELEMENT = "MONOMER_NOTATION";
	private static final String NUCLEOTIDE_ELEMENT = "NUCLEOTIDE";
	private static final String TEMPLATE_ELEMENT = "TEMPLATE";
	private static final String TEMPLATE_NOTATION_SOURCE_ATTRIBUTE = "notationSource";

	/**
	 * returns complement sequence in reverse (3-5) direction
	 * 
	 * @param sequence
	 * @return complement sequence in reverse (3-5) direction
	 * @throws org.helm.notation.NotationException
	 */
	public static String getReverseComplementSequence(String sequence)
			throws NotationException, IOException, JDOMException {
		return getComplementSequence(sequence, false);
	}

	public static String getReverseComplementSequence(String sequence,
			String notationSource) throws NotationException, IOException,
			JDOMException {
		return getComplementSequence(sequence, notationSource, false);
	}

	/**
	 * returns complement sequence in normal (5-3) direction
	 * 
	 * @param sequence
	 * @return complement sequence in normal (5-3) direction
	 * @throws org.helm.notation.NotationException
	 */
	public static String getNormalComplementSequence(String sequence)
			throws NotationException, IOException, JDOMException {
		return getComplementSequence(sequence, true);
	}

	public static String getNormalComplementSequence(String sequence,
			String notationSource) throws NotationException, IOException,
			JDOMException {
		return getComplementSequence(sequence, notationSource, true);
	}

	private static String getComplementSequence(String sequence,
			boolean outNormalDirection) throws NotationException, IOException,
			JDOMException {
		return getComplementSequence(sequence,
				NotationConstant.NOTATION_SOURCE, outNormalDirection);
	}

	/**
	 * This method generates the complement sequence
	 * 
	 * @param sequence
	 *            input sequence, could be 5-3 or 3-5
	 * @param notationSource
	 * @param outNormalDirection
	 *            output direction, true for 5-3, false for 3-5
	 * @return complement sequence
	 * @throws org.helm.notation.NotationException
	 */
	private static String getComplementSequence(String sequence,
			String notationSource, boolean outNormalDirection)
			throws NotationException, IOException, JDOMException {
		String cleanSequence = cleanup(sequence);
		List<Nucleotide> nucList = getNucleotideList(cleanSequence,
				notationSource);
		List<Nucleotide> compList = getComplementList(nucList);

		StringBuffer sb = new StringBuffer();
		boolean inNormalDirection = isNormalDirection(sequence);
		if (inNormalDirection) {
			if (outNormalDirection) {
				sb.append("5'-");
				for (int i = compList.size() - 1; i >= 0; i--) {
					sb.append(compList.get(i).getNaturalAnalog());
				}
				sb.append("-3'");
			} else {
				sb.append("3'-");
				for (int i = 0; i < compList.size(); i++) {
					sb.append(compList.get(i).getNaturalAnalog());
				}
				sb.append("-5'");
			}
		} else {
			if (outNormalDirection) {
				sb.append("5'-");
				for (int i = 0; i < compList.size(); i++) {
					sb.append(compList.get(i).getNaturalAnalog());
				}
				sb.append("-3'");
			} else {
				sb.append("3'-");
				for (int i = compList.size() - 1; i >= 0; i--) {
					sb.append(compList.get(i).getNaturalAnalog());
				}
				sb.append("-5'");
			}

		}

		return sb.toString();
	}

	/**
	 * Assume sequence input in one of the following format<br>
	 * <li>5'-AGCUdT-3'<br> <li>3'-dTUCGA-5'<br> <li>AFCUdT (default 5' to 3')
	 * <br>
	 * 
	 * @param sequence
	 *            RNA nucleotide sequence 5' to 3'
	 * @return notation
	 * @throws org.helm.notation.NotationException
	 */
	public static String getNotation(String sequence) throws NotationException,
			IOException, JDOMException {
		return getNotation(sequence, NotationConstant.NOTATION_SOURCE);
	}

	public static String getNotation(String sequence, String notationSource)
			throws NotationException, IOException, JDOMException {
		StringBuffer sb = new StringBuffer();
		List<Nucleotide> normalNucleotideList = getNormalList(sequence,
				notationSource);

		for (int i = 0; i < normalNucleotideList.size(); i++) {
			Nucleotide nucleotide = normalNucleotideList.get(i);
			if (sb.length() > 0) {
				sb.append(".");
			}
			sb.append(nucleotide.getNotation());
		}

		return sb.toString();
	}

	/**
	 * assume the sequence direction is from 5' to 3', unless started with 3
	 * 
	 * @param sequence
	 * @return true or false
	 */
	public static boolean isNormalDirection(String sequence) {
		if (sequence.startsWith("3")) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * for non-natural ones, use 'X R(X)P"
	 * 
	 * @param inList
	 * @return list of Nucleotide
	 */
	private static List<Nucleotide> getComplementList(List<Nucleotide> inList) {
		List<Nucleotide> outList = new ArrayList<Nucleotide>();
		for (int i = 0; i < inList.size(); i++) {
			Nucleotide nucleotide = inList.get(i);
			String naturalAnalog = nucleotide.getNaturalAnalog();
			String symbol = complementMap.get(naturalAnalog);
			String notation = "R(" + symbol + ")P";
			Nucleotide compNu = new Nucleotide(symbol, notation);
			outList.add(compNu);
		}
		return outList;
	}

	public static List<Nucleotide> getNormalList(String sequence)
			throws NotationException, IOException, JDOMException {
		return getNormalList(sequence, NotationConstant.NOTATION_SOURCE);
	}

	/**
	 * get the noraml nucleotide list regardness the direction of input sequence
	 * 
	 * @param sequence
	 * @param notationSource
	 * @return normal nucleotide list
	 * @throws org.helm.notation.NotationException
	 */
	public static List<Nucleotide> getNormalList(String sequence,
			String notationSource) throws NotationException, IOException,
			JDOMException {
		boolean inNormalDirection = isNormalDirection(sequence);
		String cleanSequence = cleanup(sequence);
		List<Nucleotide> nucList = getNucleotideList(cleanSequence,
				notationSource);

		if (inNormalDirection) {
			return nucList;

		} else {
			return getReverseList(nucList);
		}
	}

	/**
	 * remove non nucleotide related characters
	 * 
	 * @param sequence
	 */
	public static String cleanup(String sequence) {
		String condensedSequence = sequence.replaceAll("\\s", ""); // remove all
																	// white
																	// space
		String[] tokens = condensedSequence.split("-");
		String result = null;
		if (tokens[0].startsWith("5")) {
			result = tokens[1];
		} else if (tokens[0].startsWith("3")) {
			result = tokens[1];
		} else {
			result = tokens[0];
		}

		if (result.equals(result.toLowerCase())) {
			result = result.toUpperCase();
		}
		return result;
	}

	/**
	 * Reverse the order of list element, the last one becomes the first one,
	 * and the first one becomes the last one
	 * 
	 * @param inList
	 * @return outList
	 */
	private static List<Nucleotide> getReverseList(List<Nucleotide> inList) {
		List<Nucleotide> outList = new ArrayList<Nucleotide>();
		for (int i = inList.size() - 1; i >= 0; i--) {
			outList.add(inList.get(i));
		}
		return outList;
	}

	/**
	 * This method converts non-direction nucleotide sequence into a List of
	 * nucleotide
	 * 
	 * @param nonDirectionSeuqence
	 * @return list of nucleotides
	 * @throws org.helm.notation.NotationException
	 */
	private static List<Nucleotide> getNucleotideList(
			String nonDirectionSeuqence) throws NotationException, IOException,
			JDOMException {
		return getNucleotideList(nonDirectionSeuqence,
				NotationConstant.NOTATION_SOURCE);
	}

	/**
	 * This method converts non-direction nucleotide sequence into a List of
	 * nucleotide (sequence contains letters only, after cleanup)
	 * 
	 * @param nonDirectionSeuqence
	 * @return list of nucleotides
	 * @throws org.helm.notation.NotationException
	 */
	private static List<Nucleotide> getNucleotideList(
			String nonDirectionSeuqence, String notationSource)
			throws NotationException, IOException, JDOMException {

		if (null == nonDirectionSeuqence) {
			throw new NotationException("Nucleotide Sequence must be specified");
		}
		if (null == notationSource) {
			throw new NotationException("Notation Source must be specified");
		}
		Map<String, Map<String, String>> templates = NucleotideFactory
				.getInstance().getNucleotideTemplates();
		Map<String, String> nucleotides = null;
		if (templates.containsKey(notationSource)) {
			nucleotides = templates.get(notationSource);
		} else {
			throw new NotationException("Unknown Notation Source ["
					+ notationSource + "]");
		}
		Set<String> keySet = nucleotides.keySet();

		// walk the sequence
		List<Nucleotide> l = new ArrayList<Nucleotide>();
		int pos = 0;
		while (pos < nonDirectionSeuqence.length()) {
			boolean found = false;
			for (Iterator i = keySet.iterator(); i.hasNext();) {
				String symbol = (String) i.next();

				if (nonDirectionSeuqence.startsWith(symbol, pos)) {
					found = true;
					String notation = nucleotides.get(symbol);
					Nucleotide nuc = new Nucleotide(symbol, notation);
					l.add(nuc);
					pos = pos + symbol.length();
					break;
				}
			}
			if (!found) {
				throw new NotationException(
						"Sequence contains unknown nucleotide starting at "
								+ nonDirectionSeuqence.substring(pos));
			}
		}

		return l;
	}

	public static Map<String, Map<String, String>> getNucleotideTemplates(
			Element templatesElement) {
		Map<String, Map<String, String>> map = new HashMap<String, Map<String, String>>();

		List templates = templatesElement.getChildren();
		for (Iterator i = templates.iterator(); i.hasNext();) {
			Element templateElement = (Element) i.next();
			String notationSource = templateElement
					.getAttributeValue(TEMPLATE_NOTATION_SOURCE_ATTRIBUTE);
			Map<String, String> tmpMap = new HashMap<String, String>();
			map.put(notationSource, tmpMap);
			List nucleotides = templateElement.getChildren();
			for (Iterator it = nucleotides.iterator(); it.hasNext();) {
				Element nucleotideElement = (Element) it.next();
				Nucleotide nucleotide = getNucleotide(nucleotideElement);
				tmpMap.put(nucleotide.getSymbol(), nucleotide.getNotation());
			}
		}
		return map;
	}

	public static Nucleotide getNucleotide(Element nucleotideElement) {
		Element symbolE = nucleotideElement.getChild(NUCLEOTIDE_SYMBOL_ELEMENT,
				nucleotideElement.getNamespace());
		Element notationE = nucleotideElement.getChild(
				NUCLEOTIDE_MONOMER_NOTATION_ELEMENT,
				nucleotideElement.getNamespace());
		return new Nucleotide(symbolE.getText(), notationE.getText());
	}

	public static Nucleotide getNucleotide(String nucleotideXML)
			throws JDOMException, IOException {
		Nucleotide nuc = null;
		if (nucleotideXML != null && nucleotideXML.length() > 0) {
			SAXBuilder builder = new SAXBuilder();
			ByteArrayInputStream bais = new ByteArrayInputStream(
					nucleotideXML.getBytes());
			Document doc = builder.build(bais);
			Element root = doc.getRootElement();
			nuc = getNucleotide(root);
		}
		return nuc;
	}

	public static String getNucleotideTemplatesXML(
			Map<String, Map<String, String>> templates) {
		XMLOutputter outputer = new XMLOutputter(Format.getPrettyFormat());

		StringBuilder sb = new StringBuilder();
		sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<NUCLEOTIDE_TEMPLATES xsi:schemaLocation=\"lmr NucleotideTemplateSchema.xsd\" xmlns=\"lmr\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n");

		Set<String> templateSet = templates.keySet();
		for (Iterator i = templateSet.iterator(); i.hasNext();) {
			String template = (String) i.next();
			Element templateElement = new Element(TEMPLATE_ELEMENT);
			Attribute att = new Attribute(TEMPLATE_NOTATION_SOURCE_ATTRIBUTE,
					template);
			templateElement.setAttribute(att);

			Map<String, String> nucMap = templates.get(template);
			Set<String> nucleotideSet = nucMap.keySet();

			for (Iterator it = nucleotideSet.iterator(); it.hasNext();) {
				Element nucleotideElement = new Element(NUCLEOTIDE_ELEMENT);
				templateElement.getChildren().add(nucleotideElement);

				String symbol = (String) it.next();
				Element symbolElement = new Element(NUCLEOTIDE_SYMBOL_ELEMENT);
				symbolElement.setText(symbol);
				nucleotideElement.getChildren().add(symbolElement);

				String notation = (String) nucMap.get(symbol);
				Element notationElement = new Element(
						NUCLEOTIDE_MONOMER_NOTATION_ELEMENT);
				notationElement.setText(notation);
				nucleotideElement.getChildren().add(notationElement);
			}

			String templateString = outputer.outputString(templateElement);
			sb.append(templateString);
		}

		sb.append("\n</NUCLEOTIDE_TEMPLATES>");

		return sb.toString();
	}

	/**
	 * This method converts nucleotide sequences into HELM Notation
	 * 
	 * @param senseSeq
	 *            5-3 nucleotide sequence for Default notation
	 * @param antiSenseSeq
	 *            5-3 nucleotide sequence for default notation
	 * @return HELM notation for siRNA
	 */
	public static String getSirnaNotation(String senseSeq, String antiSenseSeq)
			throws NotationException, IOException, JDOMException {
		return getSirnaNotation(senseSeq, antiSenseSeq, RNA_DESIGN_NONE);
	}

	/**
	 * this method converts nucleotide sequences into HELM notation based on
	 * design pattern
	 * 
	 * @param senseSeq
	 * @param antiSenseSeq
	 * @param rnaDesignType
	 * @return HELM Notation for siRNA
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static String getSirnaNotation(String senseSeq, String antiSenseSeq,
			String rnaDesignType) throws NotationException, IOException,
			JDOMException {
		validateSirnaDesign(senseSeq, antiSenseSeq, rnaDesignType);

		StringBuffer sb = new StringBuffer();
		String[] tokens = new String[] { "", "" };
		int count = 0;

		if (senseSeq != null && senseSeq.length() > 0) {
			count++;
			String ssNotation = getNotation(senseSeq);
			ssNotation = removeLastP(ssNotation);
			tokens[0] = "RNA" + count + "{" + ssNotation + "}";
			tokens[1] = "RNA" + count + "{ss}";
		}

		if (antiSenseSeq != null && antiSenseSeq.length() > 0) {
			count++;
			String asNotation = getNotation(antiSenseSeq);
			asNotation = removeLastP(asNotation);
			if (tokens[0].length() > 0) {
				tokens[0] = tokens[0] + "|";
			}
			tokens[0] = tokens[0] + "RNA" + count + "{" + asNotation + "}";

			if (tokens[1].length() > 0) {
				tokens[1] = tokens[1] + "|";
			}
			tokens[1] = tokens[1] + "RNA" + count + "{as}";
		}

		String basePair = basePair = hybridization(senseSeq, antiSenseSeq,
				rnaDesignType);

		if (tokens[0].length() > 0) {
			sb.append(tokens[0]);
			sb.append("$$");
			if (basePair.length() > 0) {
				sb.append(basePair);
			}
			sb.append("$");
			sb.append(tokens[1]);
			sb.append("$");
		}

		return sb.toString();
	}

	private static String hybridization(String senseSeq, String antiSenseSeq,
			String rnaDesignType) throws NotationException, IOException,
			JDOMException {
		String basePair = "";
		if (senseSeq != null && senseSeq.length() > 0 && antiSenseSeq != null
				&& antiSenseSeq.length() > 0) {
			String analogSeqSS = getNaturalAnalogSequence(senseSeq).replaceAll(
					"T", "U");
			String analogSeqAS = getNaturalAnalogSequence(antiSenseSeq)
					.replaceAll("T", "U");

			if (RNA_DESIGN_NONE.equalsIgnoreCase(rnaDesignType)) {
				String normalCompAS = cleanup(getNormalComplementSequence(analogSeqAS));
				String maxMatch = getMaxMatchFragment(analogSeqSS, normalCompAS);
				if (maxMatch.length() > 0) {
					int ssStart = analogSeqSS.indexOf(maxMatch);
					int normalCompStart = normalCompAS.indexOf(maxMatch);
					int asStart = analogSeqAS.length() - maxMatch.length()
							- normalCompStart;

					for (int i = 0; i < maxMatch.length(); i++) {
						int ssPos = (i + ssStart) * 3 + 2;
						int asPos = (asStart + maxMatch.length() - 1 - i) * 3 + 2;
						if (basePair.length() > 0) {
							basePair = basePair + "|";
						}
						basePair = basePair + "RNA1,RNA2," + ssPos + ":pair-"
								+ asPos + ":pair";
					}
				}
			} else if (RNA_DESIGN_TUSCHL_19_PLUS_2
					.equalsIgnoreCase(rnaDesignType)) {
				int matchLength = 19;
				basePair = hybridizationWithLengthFromStart(analogSeqSS,
						analogSeqAS, matchLength);
			} else if (RNA_DESIGN_DICER_27_R.equalsIgnoreCase(rnaDesignType)) {
				int matchLength = 25;
				basePair = hybridizationWithLengthFromStart(analogSeqSS,
						analogSeqAS, matchLength);
			} else if (RNA_DESIGN_DICER_27_L.equalsIgnoreCase(rnaDesignType)) {
				int matchLength = 25;
				basePair = hybridizationWithLengthFromStart(analogSeqSS,
						analogSeqAS, matchLength);
			} else {
				throw new NotationException("Unknow RNA Design");
			}
		}
		return basePair;
	}

	private static String hybridizationWithLengthFromStart(
			String senseAnalogSeq, String antisenseAnalogSeq,
			int lengthFromStart) {
		String basePair = "";
		for (int i = 0; i < lengthFromStart; i++) {
			int ssPos = i * 3 + 2;
			int asPos = (lengthFromStart - 1 - i) * 3 + 2;
			String ssChar = String.valueOf(senseAnalogSeq.charAt(i));
			String asChar = String.valueOf(antisenseAnalogSeq
					.charAt(lengthFromStart - 1 - i));

			if (complementMap.get(ssChar).equalsIgnoreCase(asChar)) {
				if (basePair.length() > 0) {
					basePair = basePair + "|";
				}
				basePair = basePair + "RNA1,RNA2," + ssPos + ":pair-" + asPos
						+ ":pair";
			}
		}
		return basePair;
	}

	private static boolean validateSirnaDesign(String senseSeq,
			String antisenseSeq, String rnaDesignType)
			throws NotationException, IOException, JDOMException {
		if (RNA_DESIGN_NONE.equalsIgnoreCase(rnaDesignType)) {
			return true;
		}

		if (!SUPPORTED_DESIGN_LIST.contains(rnaDesignType)) {
			throw new NotationException("Unsupported RNA Design Type '"
					+ rnaDesignType + "'");
		}

		List<Nucleotide> senseNucList = getNucleotideList(senseSeq);
		List<Nucleotide> antisenseNucList = getNucleotideList(antisenseSeq);

		if (rnaDesignType.equals(RNA_DESIGN_TUSCHL_19_PLUS_2)) {
			if (senseNucList.size() != 21) {
				throw new NotationException(
						"Sense strand for Tuschl 19+2 design must have 21 nucleotides");
			}
			if (antisenseNucList.size() != 21) {
				throw new NotationException(
						"Antisense strand for Tuschl 19+2 design must have 21 nucleotides");
			}
		} else if (rnaDesignType.equals(RNA_DESIGN_DICER_27_R)) {
			if (senseNucList.size() != 25) {
				throw new NotationException(
						"Sense strand for Dicer 27R design must have 25 nucleotides");
			}
			if (antisenseNucList.size() != 27) {
				throw new NotationException(
						"Antisense strand for Dicer 27R design must have 27 nucleotides");
			}
		} else if (rnaDesignType.equals(RNA_DESIGN_DICER_27_L)) {
			if (senseNucList.size() != 27) {
				throw new NotationException(
						"Sense strand for Dicer 27L design must have 27 nucleotides");
			}
			if (antisenseNucList.size() != 25) {
				throw new NotationException(
						"Antisense strand for Dicer 27L design must have 25 nucleotides");
			}
		}

		return true;
	}

	/**
	 * This method cleans up the last phosphate
	 * 
	 * @param simpleRNANotation
	 *            - the simple notation for single RNA sequence
	 * @return the simple notation for RNA sequence with the last P removed
	 */
	public static String removeLastP(String simpleRNANotation) {
		String result = null;
		if (simpleRNANotation != null && simpleRNANotation.endsWith("P")) {
			int len = simpleRNANotation.length();
			result = simpleRNANotation.substring(0, len - 1);
		} else {
			result = simpleRNANotation;
		}
		return result;
	}

	/**
	 * This method returns the sequence in reverse direction from the input
	 * sequence, ie. AAGCU => UCGAA
	 * 
	 * @param sequence
	 *            - input nucleotide sequence
	 * @return sequence in opposite direction (from right to left)
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static String getReverseSequence(String sequence)
			throws NotationException, IOException, JDOMException {
		List<Nucleotide> l = getNucleotideList(sequence);
		List<Nucleotide> rl = getReverseList(l);
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < rl.size(); i++) {
			sb.append(rl.get(i).getSymbol());
		}
		return sb.toString();
	}

	/**
	 * convert modified sequence to natural analog sequence
	 * 
	 * @param sequence
	 *            modified sequence
	 * @return single letter natural analog sequence
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static String getNaturalAnalogSequence(String sequence)
			throws NotationException, IOException, JDOMException {
		List<Nucleotide> l = getNucleotideList(sequence);
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < l.size(); i++) {
			sb.append(l.get(i).getNaturalAnalog());
		}
		return sb.toString();
	}

	/**
	 * 
	 * @param senseSeq
	 *            - siRNA sense strand sequence, generally in 5-3 direction
	 * @param antiSenseSeq
	 *            - siRNA antisense strand sequence, could be in 5-3 (same
	 *            direction) or 3-5 (opposition direction)
	 * @return true if opposite, false is same
	 * @throws org.helm.notation.NotationException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static boolean isInOppositeDirection(String senseSeq,
			String antiSenseSeq) throws NotationException, IOException,
			JDOMException {
		String upperSS = senseSeq.toUpperCase();
		String upperAS = antiSenseSeq.toUpperCase();
		String normalCompAS = NucleotideSequenceParser
				.getNormalComplementSequence(upperAS);
		String reverseCompAS = NucleotideSequenceParser
				.getReverseComplementSequence(upperAS);

		String normalMax = getMaxMatchFragment(upperSS, normalCompAS);
		String reverseMax = getMaxMatchFragment(upperSS, reverseCompAS);

		if (normalMax.length() < reverseMax.length()) {
			return true;
		} else {
			return false;
		}
	}

	public static String getMaxMatchFragment(String seq1, String seq2)
			throws NotationException {
		return getMaxMatchFragment(seq1, seq2, MINUMUM_MATCH_FRAGMENT_LENGTH);
	}

	/**
	 * This method returns the largest matched fragment between two sequences,
	 * replace T with U before match
	 * 
	 * @param seq1
	 *            single letter, all upper case nucleotide sequence
	 * @param seq2
	 *            single letter, all upper case nucleotide sequence
	 * @param minLength
	 *            - minimum fragment length
	 * @return largest match fragment
	 */
	public static String getMaxMatchFragment(String seq1, String seq2,
			int minLength) throws NotationException {
		if (null == seq1 || null == seq2) {
			throw new NotationException("Both sequences must not be null ");
		}

		if (!seq1.equals(seq1.toUpperCase())
				|| !seq2.equals(seq2.toUpperCase())) {
			throw new NotationException(
					"Both sequences must be natural nucleotide sequence in upper case ");
		}

		String longSeq, shortSeq;

		if (seq1.length() > seq2.length()) {
			longSeq = seq1;
			shortSeq = seq2;
		} else {
			longSeq = seq2;
			shortSeq = seq1;
		}
		// replace T with U
		longSeq = longSeq.replaceAll("T", "U");
		shortSeq = shortSeq.replaceAll("T", "U");

		int min = MINUMUM_MATCH_FRAGMENT_LENGTH;
		if (minLength > min) {
			min = minLength;
		}

		for (int len = shortSeq.length(); len > min; len--) {
			for (int i = 0; i <= shortSeq.length() - len; i++) {
				String tmp = shortSeq.substring(i, i + len);

				if (longSeq.contains(tmp)) {
					return tmp;
				}
			}
		}

		return "";
	}
}
