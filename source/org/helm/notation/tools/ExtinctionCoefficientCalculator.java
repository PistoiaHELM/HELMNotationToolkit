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

import org.helm.notation.CalculationException;
import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.MonomerStore;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.PolymerNode;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import org.jdom.JDOMException;

/**
 *
 * @author ZHANGTIANHONG
 */
public class ExtinctionCoefficientCalculator {

    public static final int RNA_UNIT_TYPE = 1;
    public static final int PEPTIDE_UNIT_TYPE = 2;
    private static Map<String, Float> monoNucleotideMap = new HashMap<String, Float>();
    private static Map<String, Float> diNucleotideMap = new HashMap<String, Float>();
    private static Map<String, Float> aminoAcidMap = new HashMap<String, Float>();
    private static String rnaPropertyFile = "/org/helm/notation/resources/RNAExtinctionCoefficient.properties";
    private static String peptidePropertyFile = "/org/helm/notation/resources/PEPTIDEExtinctionCoefficient.properties";
    private static ExtinctionCoefficientCalculator instance;

    private ExtinctionCoefficientCalculator() {
    }

    public static ExtinctionCoefficientCalculator getInstance() throws CalculationException {
        if (null == instance) {
            instance = new ExtinctionCoefficientCalculator();
            try {
                initMaps();
            } catch (IOException ex) {
                throw new CalculationException("Unable to initialize extinction coefficient property files");
            }
        }
        return instance;
    }

    public String getUnit(int unitType) {
        switch (unitType) {
            case RNA_UNIT_TYPE:
                return getRnaUnit();
            case PEPTIDE_UNIT_TYPE:
                return getPeptideUnit();
            default:
                return null;
        }
    }

    public String getRnaUnit() {
        return "mM-1cm-1";
    }

    public String getPeptideUnit() {
        return "M-1cm-1";
    }
    
    public String getDefaultUnit() {
        return getRnaUnit();
    }
    
    public int getDefaultUnitType() {
        return RNA_UNIT_TYPE;
    }

    private static void initMaps() throws IOException {
        InputStream ris = ExtinctionCoefficientCalculator.class.getResourceAsStream(rnaPropertyFile);
        Properties rp = new Properties();
        rp.load(ris);

        Enumeration re = rp.propertyNames();
        while (re.hasMoreElements()) {
            String key = (String) re.nextElement();
            String value = rp.getProperty(key);
            Float f = new Float(value);
            int len = key.length();
            if (len == 1) {
                monoNucleotideMap.put(key, f);
            } else if (len == 2) {
                diNucleotideMap.put(key, f);
            }
        }
        ris.close();

        InputStream pis = ExtinctionCoefficientCalculator.class.getResourceAsStream(peptidePropertyFile);
        Properties pp = new Properties();
        pp.load(pis);

        Enumeration pe = pp.propertyNames();
        while (pe.hasMoreElements()) {
            String key = (String) pe.nextElement();
            String value = pp.getProperty(key);

            Float f = new Float(value);
            aminoAcidMap.put(key, f);
        }
        pis.close();
    }

    /**
     * This method calculates extinction coefficient for complex polymer notation
     * @param complexNotation - complex polymer notation
     * @return extinction coefficient in RNA UNIT mM-1cm-1
     * @throws org.helm.notation.NotationException
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     * @throws org.jdom.JDOMException
     * @throws org.helm.notation.StructureException
     * @throws org.helm.notation.CalculationException
     */
    public float calculateFromComplexNotation(String complexNotation) throws NotationException, MonomerException, IOException, JDOMException, StructureException, CalculationException {
        return calculateFromComplexNotation(complexNotation, getDefaultUnitType());
    }

    
    public float calculateFromComplexNotation(String complexNotation, int unitType) throws NotationException, MonomerException, IOException, JDOMException, StructureException, CalculationException {
    	MonomerFactory factory = null;
		try {
			factory = MonomerFactory.getInstance();
		} catch (Exception ex) {
			throw new NotationException("Unable to initialize monomer factory",
					ex);
		}
    	return calculateFromComplexNotation(complexNotation,unitType,factory.getMonomerStore());
    }
        
    
    /**
     * This method calculates extinction coefficient for complex polymer notation
     * @param complexNotation
     * @param unitType either RNA or PEPTIDE
     * @return EC with give unit type
     * @throws NotationException
     * @throws MonomerException
     * @throws IOException
     * @throws JDOMException
     * @throws StructureException
     * @throws CalculationException 
     */
    
    public float calculateFromComplexNotation(String complexNotation, int unitType,MonomerStore monomerStore) throws NotationException, MonomerException, IOException, JDOMException, StructureException, CalculationException {
        
    	float result = 0.0f;
        List<PolymerNode> polymerNodes = ComplexNotationParser.getPolymerNodeList(complexNotation,monomerStore);
        for (PolymerNode polymerNode : polymerNodes) {
            String polymerType = polymerNode.getType();
            String notation = polymerNode.getLabel();
            float ext = 0.0f;
            if (polymerType.equals(Monomer.NUCLIEC_ACID_POLYMER_TYPE)) {
                ext = calculateFromRnaPolymerNotation(notation,monomerStore);
                if (unitType == PEPTIDE_UNIT_TYPE) {
                    ext = ext*1000;
                }
            } else if (polymerType.equals(Monomer.PEPTIDE_POLYMER_TYPE)) {
                ext = calculateFromPeptidePolymerNotation(notation,monomerStore);
                if (unitType == RNA_UNIT_TYPE) {
                    ext = ext/1000;
                }
            }
            result = result + ext;
        }
        return result;
    }

    /**
     * This method calculates extinction coefficient for RNA polymer notation
     * @param simpleNotation RNA polymer notation
     * @return extinction coefficient
     * @throws org.helm.notation.NotationException
     * @throws org.helm.notation.MonomerException
     * @throws java.io.IOException
     * @throws org.jdom.JDOMException
     * @throws org.helm.notation.StructureException
     * @throws org.helm.notation.CalculationException
     */
    
    public float calculateFromRnaPolymerNotation(String simpleNotation) throws NotationException, MonomerException, CalculationException, IOException, JDOMException, StructureException {
        MonomerStore store=MonomerFactory.getInstance().getMonomerStore();
        return calculateFromRnaPolymerNotation(simpleNotation,store);
        
    }
    
    public float calculateFromRnaPolymerNotation(String simpleNotation,MonomerStore monomerStore) throws NotationException, MonomerException, CalculationException, IOException, JDOMException, StructureException {
        String naturalSequence = SimpleNotationParser.getTrimmedNucleotideSequence(simpleNotation,monomerStore);
        return calculateFromNucleotideSequence(naturalSequence);
    }

    /**
     * This method calculates extinction coefficient for modified nucleotide sequence
     * @param sequence modified nucleotide sequence
     * @return extinction coefficient
     * @throws org.helm.notation.CalculationException
     * @throws org.helm.notation.NotationException
     * @throws java.io.IOException
     * @throws org.jdom.JDOMException
     */
    public float calculateFromModifiedNucleotideSequence(String sequence) throws CalculationException, NotationException, IOException, JDOMException {
        String naturalSequence = NucleotideSequenceParser.getNaturalAnalogSequence(sequence);
        return calculateFromNucleotideSequence(naturalSequence);
    }

    /**
     * This method calculates extinction coefficient for natural nucleotide sequence
     * @param sequence natural nucletodie seuqnce
     * @return extinction coefficient
     * @throws org.helm.notation.CalculationException
     */
    public float calculateFromNucleotideSequence(String sequence) throws CalculationException {
        float result = 0.0f;
        if (null == sequence || sequence.length() == 0) {
            throw new CalculationException("Input sequence cannot be null");
        } else if (sequence.length() == 1) {
            if (monoNucleotideMap.containsKey(sequence)) {
                result = monoNucleotideMap.get(sequence).floatValue();
            } else {
                throw new CalculationException("Unknown nucleotide [" + sequence + "] found");
            }
        } else {
            Map<String, Integer> monoMap = new HashMap<String, Integer>();
            String mono;
            for (int i = 1; i < sequence.length() - 1; i++) {
                mono = sequence.substring(i, i + 1);
                if (monoNucleotideMap.containsKey(mono)) {
                    if (monoMap.containsKey(mono)) {
                        int newInt = monoMap.get(mono).intValue() + 1;
                        monoMap.put(mono, new Integer(newInt));
                    } else {
                        monoMap.put(mono, new Integer(1));
                    }
                } else {
                    throw new CalculationException("Unknown nucleotide [" + mono + "] found");
                }
            }

            Map<String, Integer> diMap = new HashMap<String, Integer>();
            String di;
            for (int i = 0; i < sequence.length() - 1; i++) {
                di = sequence.substring(i, i + 2);
                if (diNucleotideMap.containsKey(di)) {
                    if (diMap.containsKey(di)) {
                        int newInt = diMap.get(di).intValue() + 1;
                        diMap.put(di, new Integer(newInt));
                    } else {
                        diMap.put(di, new Integer(1));
                    }
                } else {
                    throw new CalculationException("Unknown dinucleotide [" + di + "] found");
                }
            }

            result = calculate(monoMap, diMap);
        }

        return result;

    }

    private float calculate(Map<String, Integer> tmpMonoNucleotideMap, Map<String, Integer> tmpDiNucleotideMap) {
        float monoResult = 0.0f;
        Set<String> monoSet = tmpMonoNucleotideMap.keySet();
        for (Iterator i = monoSet.iterator(); i.hasNext();) {
            String key = (String) i.next();
            Integer count = tmpMonoNucleotideMap.get(key);
            Float value = monoNucleotideMap.get(key);
            monoResult = monoResult + (value.floatValue()) * (count.floatValue());
        }

        float diResult = 0.0f;
        Set<String> diSet = tmpDiNucleotideMap.keySet();
        for (Iterator i = diSet.iterator(); i.hasNext();) {
            String key = (String) i.next();
            Integer count = tmpDiNucleotideMap.get(key);
            Float value = diNucleotideMap.get(key);
            diResult = diResult + (value.floatValue()) * (count.floatValue());
        }

        return 2 * diResult - monoResult;
    }

    /**
     * This method returns the extinction coefficient of single letter amino acid sequence
     * @param sequence - single letter amino acid sequence
     * @return extinction coefficient of natural amino acid sequence
     * @throws CalculationException 
     */
    public float calculateFromAminoAcidSequence(String sequence) throws CalculationException {
        if (null == sequence || sequence.length() == 0) {
            return 0.0f;
        }
        List<String> ids = new ArrayList<String>();
        if (null != sequence) {
            char[] chars = sequence.toCharArray();
            for (int i = 0; i < chars.length; i++) {
                String id = String.valueOf(chars[i]);
                ids.add(id);
            }
        }
        return calculate(ids);
    }
    
    
    public float calculateFromPeptidePolymerNotation(String simpleNotation) throws NotationException, MonomerException, CalculationException, IOException, JDOMException, StructureException {
    	MonomerFactory factory = null;
    	try {
    		factory = MonomerFactory.getInstance();
    	} catch (Exception ex) {
    		throw new NotationException("Unable to initialize monomer factory",
    				ex);
    	}
    	return calculateFromPeptidePolymerNotation(simpleNotation,factory.getMonomerStore());
    }
    

    public float calculateFromPeptidePolymerNotation(String simpleNotation,MonomerStore monomerStore) throws NotationException, MonomerException, CalculationException, IOException, JDOMException, StructureException {
        List<Monomer> monomers = SimpleNotationParser.getMonomerList(simpleNotation, Monomer.PEPTIDE_POLYMER_TYPE,monomerStore);
        List<String> ids = new ArrayList<String>();
        for (Monomer monomer : monomers) {
            String id = monomer.getNaturalAnalog();
            ids.add(id);
        }
        return calculate(ids);
    }

    private float calculate(List<String> aaIDlist) {
        if (null == aaIDlist || aaIDlist.isEmpty()) {
            return 0.0f;
        }

        Map<String, Integer> countMap = new HashMap<String, Integer>();
        for (String id : aaIDlist) {
            if (aminoAcidMap.containsKey(id)) {
                int count = 1;
                if (countMap.containsKey(id)) {
                    count = count + countMap.get(id);
                }
                countMap.put(id, count);
            }
        }

        float result = 0.0f;
        Set<String> keys = countMap.keySet();
        for (Iterator<String> it = keys.iterator(); it.hasNext();) {
            String key = it.next();
            int count = countMap.get(key);
            float factor = aminoAcidMap.get(key);
            result = result + factor * count;
        }

        return result;
    }
}
