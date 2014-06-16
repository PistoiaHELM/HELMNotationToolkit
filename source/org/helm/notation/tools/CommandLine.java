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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;

/**
 * Commonand line tools for different conversion, this is the manin class in
 * manifest file
 * 
 * @author zhangtianhong
 */
public class CommandLine {
	public static final String[] options = { "seq2notation",
			"notation2property" };

	public static void main(String[] args) {
		try {

			if (args.length != 3) {
				System.out
						.println("Usage: java -jar NotationToolkit.jar conversion_option[seq2notation|notation2property] input_file output_file\n");
				System.out
						.println("Example: java -jar NotationToolkit.jar seq2notation c:/data/seq.txt c:/data/notation.txt\n");
				System.exit(0);
			}

			if (!isValidOption(args[0])) {
				System.out.println("The conversion option [" + args[0]
						+ "] is not supported");
				System.exit(0);
			}

			if (args[0].equalsIgnoreCase(options[0])) {

				File infile = new File(args[1]);
				File outfile = new File(args[2]);

				FileOutputStream outStream = new FileOutputStream(outfile);
				BufferedReader inReader = new BufferedReader(new FileReader(
						infile));
				String line = inReader.readLine();
				int i = 0;
				while (null != line) {
					i++;
					String notation = "";
					try {
						notation = NucleotideSequenceParser.getNotation(line);
					} catch (Exception e) {
						notation = "Invalid Sequence";
					}
					String tmp = line + "\t" + notation + "\n";
					outStream.write(tmp.getBytes());
					System.out.print("" + i + "\t" + tmp);
					line = inReader.readLine();
				}

				inReader.close();
				outStream.close();
			} else {
				System.out.println("The conversion option [" + args[0]
						+ "] is not supported yet");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static boolean isValidOption(String option) {
		for (String choice : options) {
			if (choice.equalsIgnoreCase(option)) {
				return true;
			}
		}
		return false;
	}
}
