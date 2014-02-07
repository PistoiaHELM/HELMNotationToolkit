package org.helm.notation.tools;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.helm.notation.CalculationException;
import org.helm.notation.MonomerException;
import org.helm.notation.NotationException;
import org.helm.notation.StructureException;
import org.jdom.JDOMException;
import org.junit.Before;
import org.junit.Test;


public class CalculatorTest {
	
	private ExtinctionCoefficientCalculator calculator;
	
	@Before
	public void setUp() throws CalculationException{
		calculator = ExtinctionCoefficientCalculator.getInstance();
	}
	
	@Test
	public void testCalculateFromAminoAcidSequence() throws CalculationException {

		String input = "AGGDDDDDDDDDDDDDDDDDDFFFFFFFFFFFFF";
		float result = calculator.calculateFromAminoAcidSequence(input);
		assertEquals(result,0.0,1e-15);
		
		
		input = "AGGCFFFFFFFFFF";
        result = calculator.calculateFromAminoAcidSequence(input);
        assertEquals(result,0.0,62.5);
		
        
        input = "AGGYEEEEEEEEEEEEEEEEEEE";
        result = calculator.calculateFromAminoAcidSequence(input);
        assertEquals(result,1490.0,1e-15);
		
        
        
        input = "AGGWEEEEEEEEEEEEEEEEEEE";
        result = calculator.calculateFromAminoAcidSequence(input);
        assertEquals(result,5500.0,1e-15);
	}
	
	
	public void testCalculateFromPeptidePolymerNotation() throws NotationException, MonomerException, CalculationException, IOException, JDOMException, StructureException{
		String input = "A.G.G.W.E.E.E.E.E.W";
        float result = calculator.calculateFromPeptidePolymerNotation(input);
        assertEquals(result,11000.0,1e-15);
        
		
	}
	
	
	
	

}
