package org.helm.notation.tools;

import static org.junit.Assert.*;

import org.junit.Test;

import com.pfizer.pgrd.sdlib.EncoderException;

public class EncoderTest {

	@Test
	public void testEncode() throws EncoderException {
		            String molfile ="\n  Marvin  08200815002D          \n\n"+
		" 10  9  0  0  0  0            999 V2000\n"+
		"   -2.7541    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -3.4686    0.9102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -3.4686    1.7352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -4.1830    2.1477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -2.0397    1.7351    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"+
		"   -1.3250    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -1.3250    2.9726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -2.0397    0.9101    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -0.6105    1.7350    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"   -2.7542    0.4976    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n"+
		"  3  1  1  0  0  0  0\n  5  1  1  1  0  0  0\n  3  2  1  0  0  0  0\n"+
		"  3  4  1  0  0  0  0\n  5  8  1  0  0  0  0\n  5  6  1  0  0  0  0\n"+
		"  6  7  2  0  0  0  0\n  6  9  1  0  0  0  0\n  8 10  1  0  0  0  0\n"+
		"M  RGP  2   9   2  10   1\nM  END\n";    
		            
		System.out.println(molfile);		            
        String encodedMolfile = MolfileEncoder.encode(molfile);
        System.out.println(encodedMolfile);

        //String encodedString = "H4sIAAAAAAAAAO3UsQrCMBAG4P2e4hed5S612oBLacSpIkVcOjm6ODj4/N61oJbLG9gjhMt34RIyZL/om1Rf6r69PV/3BwFjAnAMJVdSinDCJ4ggAVJo/Wd8I8aIa2BmshWvNeNs1mDaIj/mLnOXMeuWf3sXAcIwT1Sp8Kq08apUelXaepWs6t6dV6XKq1L0GiHsVEnEq9gnkznNawt0x/PwGnaoVa0hglUOp0S00qA3z+wuYecEAAA=";
        
        String decodeString = MolfileEncoder.decode(encodedMolfile);
        System.out.println(decodeString);
        assertEquals(molfile,decodeString);
		                  
		 } 
		
		
}


