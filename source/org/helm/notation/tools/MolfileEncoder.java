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

import com.pfizer.pgrd.sdlib.EncoderException;
import com.pfizer.pgrd.sdlib.SDUnit;
import com.pfizer.pgrd.sdlib.SDEncoder;

public class MolfileEncoder {

	public static String encode(String string) throws EncoderException {
		if (null != string) {
			SDUnit unit = new SDUnit(string);
			return SDEncoder.encode(unit, SDEncoder.TYPE_GZIP_ENCODED);
		} else {
			return null;
		}
	}

	public static String decode(String encodedString) throws EncoderException {
		if (null != encodedString) {
			SDUnit decodedUnit = SDEncoder.decode(encodedString,
					SDEncoder.TYPE_GZIP_ENCODED);
			return decodedUnit.getMACCSData();
		} else {
			return null;
		}
	}
}
