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
package org.helm.notation.demo.tools;

import org.helm.notation.tools.*;
import org.helm.notation.MonomerFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class MonomerDBSample {

	/**
	 * @param args
	 *            the command line arguments
	 */
	public static void main(String[] args) {

		try {
			// show all monomers
			Map monomerDB = MonomerFactory.getInstance().getMonomerDB();
			showMonomerDB(monomerDB);

			// show all attachment
			Map attachmentDB = MonomerFactory.getInstance().getAttachmentDB();
			Set keyset = attachmentDB.keySet();
			for (Iterator i = keyset.iterator(); i.hasNext();) {
				String attachID = (String) i.next();
				Attachment m = (Attachment) attachmentDB.get(attachID);
				System.out.println(attachID + ": " + m.getCapGroupSMILES());
			}

			System.out.println("Serializing monomer DB............");
			File f = new File(MonomerFactory.NOTATION_DIRECTORY);
			if (!f.exists()) {
				f.mkdir();
			}
			MonomerFactory.getInstance().saveMonomerCache();
			System.out.println("Serializing monomer DB............Done");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void showMonomerDB(Map monomerDB) throws StructureException,
			IOException {
		Set polymerSet = monomerDB.keySet();
		for (Iterator it = polymerSet.iterator(); it.hasNext();) {
			String polymer = (String) it.next();
			System.out.println("Polymer Type: " + polymer);

			Map nucMap = (Map) monomerDB.get(polymer);

			Set keyset = nucMap.keySet();
			for (Iterator i = keyset.iterator(); i.hasNext();) {
				String monomerID = (String) i.next();
				Monomer m = (Monomer) nucMap.get(monomerID);
				String smi = m.getCanSMILES();
				if (null != smi) {
					String usmi = StructureParser.getUniqueExtendedSMILES(smi);
					System.out.println(monomerID);
					System.out.println(smi);
					System.out.println(usmi);
				}
			}

		}
	}
}
