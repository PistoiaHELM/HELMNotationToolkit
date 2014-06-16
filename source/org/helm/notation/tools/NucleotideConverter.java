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

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.NotationConstant;
import org.helm.notation.NotationException;
import org.helm.notation.NucleotideFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import org.helm.notation.model.PolymerNode;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jdom.JDOMException;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class NucleotideConverter {

	private static NucleotideConverter instance;

	private NucleotideConverter() {
	}

	public static NucleotideConverter getInstance() throws IOException,
			JDOMException, NotationException, MonomerException {
		if (null == instance) {
			MonomerFactory.getInstance();
			NucleotideFactory.getInstance();
			instance = new NucleotideConverter();
		}
		return instance;
	}

	/**
	 * This method returns the nucleotide sequences for nucleic acid polymer
	 * mixture
	 * 
	 * @param complexNotation
	 *            -- HELM notation for nucleic acid polymer, could be mixture
	 * @return nucleotide sequence for nucleic acid, separated with space
	 *         between strands
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */
	public String getNucleotideSequencesFromComplexNotation(
			String complexNotation) throws NotationException, MonomerException,
			IOException, JDOMException, StructureException {
		List<PolymerNode> polymerNodes = ComplexNotationParser
				.getPolymerNodeList(complexNotation);

		StringBuffer sb = new StringBuffer();
		for (PolymerNode node : polymerNodes) {
			String polymerType = node.getType();
			if (!polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
				throw new NotationException(
						"Input complex notation contains non-nucleic acid polymer");
			}
			String simpleNotation = node.getLabel();
			String notation = getNucleotideSequenceFromSimpleRNANotation(simpleNotation);
			if (sb.length() > 0) {
				sb.append(" ");
			}
			sb.append(notation);
		}
		return sb.toString();
	}

	/**
	 * This method returns the nucleotide sequence for a single simple RNA
	 * notation
	 * 
	 * @param simpleRNANotation
	 * @return nucleotide sequence
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */
	public String getNucleotideSequenceFromSimpleRNANotation(
			String simpleRNANotation) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		List<Nucleotide> nucList = SimpleNotationParser.getNucleotideList(
				simpleRNANotation, false);

		int count = 0;
		StringBuffer sb = new StringBuffer();
		Map<String, String> reverseNucMap = NucleotideFactory.getInstance()
				.getReverseNucleotideTemplateMap();
		for (Nucleotide nuc : nucList) {
			String nucleotide = nuc.getNotation();
			String nucleoside = nuc.getNucleosideNotation();
			String linker = nuc.getLinkerNotation();

			// it is ok for the first nucleotide not to have a nucleoside
			if (count == 0 && nucleoside.length() == 0) {
				sb.append(nuc.getPhosphateMonomer().getAlternateId());
				count++;
				continue;
			}

			// it is ok for the last nucleotide not to have a linker
			if (count == nucList.size() - 1 && linker.length() == 0) {
				nucleotide = nucleotide + Monomer.ID_P;
			}

			if (reverseNucMap.containsKey(nucleotide)) {
				sb.append(reverseNucMap.get(nucleotide));
			} else {
				throw new NotationException("Unknown nucleotide found for "
						+ nucleotide + " : missing nucleotide template");
			}

			count++;
		}
		return sb.toString();
	}

	/**
	 * This method return the HELM complex notation for nucleotide sequences,
	 * which could have more than one strand separated by space
	 * 
	 * @param nucleotideSequences
	 * @return complex notatoion
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */
	public String getComplexNotation(String nucleotideSequences)
			throws NotationException, MonomerException, IOException,
			JDOMException, StructureException {
		if (null == nucleotideSequences
				|| nucleotideSequences.trim().length() == 0) {
			return null;
		}

		String result = null;
		String[] singles = nucleotideSequences.split("\\s");
		if (singles.length == 1) {
			result = getComplexNotationFromSingleNucleotideSequence(singles[0]);
		} else {
			result = getComplexNotationFromSingleNucleotideSequence(singles[0]);
			for (int i = 1; i < singles.length; i++) {
				String notation = getComplexNotationFromSingleNucleotideSequence(singles[i]);
				result = ComplexNotationParser.getCombinedComlexNotation(
						result, notation);
			}
		}

		result = ComplexNotationParser.hybridize(result);

		return result;
	}

	/**
	 * This method returns the complex notation for a single oligonucleotide
	 * sequence
	 * 
	 * @param singleOligoSequence
	 * @return complex notation
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 * @throws StructureException
	 */
	public String getComplexNotationFromSingleNucleotideSequence(
			String singleOligoSequence) throws NotationException,
			MonomerException, IOException, JDOMException, StructureException {
		String simpleNotation = getSimpleRNANotationFromSingleNucleotideSequence(singleOligoSequence);
		return SimpleNotationParser.getComplextNotationForRNA(simpleNotation);
	}

	/**
	 * This method returns the simple RNA notation for a single oligonucleotide
	 * sequence
	 * 
	 * @param singleOligoSequence
	 * @return simple RNA notation
	 * @throws NotationException
	 * @throws MonomerException
	 * @throws IOException
	 * @throws JDOMException
	 */
	public String getSimpleRNANotationFromSingleNucleotideSequence(
			String singleOligoSequence) throws NotationException,
			MonomerException, IOException, JDOMException {
		Map<String, Monomer> rnaMonomers = MonomerFactory.getInstance()
				.getMonomerDB().get(Monomer.NUCLIEC_ACID_POLYMER_TYPE);
		Map<String, Monomer> linkerMonomers = new HashMap<String, Monomer>();
		Set<String> rnaSet = rnaMonomers.keySet();
		for (String key : rnaSet) {
			Monomer m = rnaMonomers.get(key);
			if (m.getNaturalAnalog().equals(Monomer.ID_P)) {
				linkerMonomers.put(key, m);
			}
		}

		Set<String> linkerSet = linkerMonomers.keySet();

		StringBuffer sb = new StringBuffer();
		String seq = singleOligoSequence;
		for (String key : linkerSet) {
			if (seq.startsWith(key)) {
				if (key.length() > 1)
					sb.append(SimpleNotationParser.MODIFICATION_START_SYMBOL
							+ key
							+ SimpleNotationParser.MODIFICATION_END_SYMBOL);
				else
					sb.append(key);
				seq = seq.substring(key.length());
				break;
			}
		}

		Map<String, String> nucleotides = NucleotideFactory.getInstance()
				.getNucleotideTemplates().get(NotationConstant.NOTATION_SOURCE);
		Set<String> nucSet = nucleotides.keySet();
		while (seq.length() > 0) {
			if (sb.length() > 0) {
				sb.append(".");
			}

			boolean found = false;
			for (String nuc : nucSet) {
				if (seq.startsWith(nuc)) {
					sb.append(nucleotides.get(nuc));
					seq = seq.substring(nuc.length());
					found = true;
					break;
				}
			}

			if (!found) {
				throw new NotationException(
						"Unknown nucleotide found starting at " + seq);
			}
		}

		return NucleotideSequenceParser.removeLastP(sb.toString());
	}
}
