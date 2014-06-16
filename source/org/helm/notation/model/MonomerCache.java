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
package org.helm.notation.model;

import java.io.Serializable;
import java.util.Map;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class MonomerCache implements Serializable {

	private Map<String, Map<String, Monomer>> monomerDB;

	private Map<String, Attachment> attachmentDB;

	private Map<String, Monomer> smilesMonomerDB;

	public Map<String, Map<String, Monomer>> getMonomerDB() {
		return monomerDB;
	}

	public Map<String, Attachment> getAttachmentDB() {
		return attachmentDB;
	}

	public Map<String, Monomer> getSmilesMonomerDB() {
		return smilesMonomerDB;
	}

	public void setMonomerDB(Map<String, Map<String, Monomer>> monomerDB) {
		this.monomerDB = monomerDB;
	}

	public void setAttachmentDB(Map<String, Attachment> attachmentDB) {
		this.attachmentDB = attachmentDB;
	}

	public void setSmilesMonomerDB(Map<String, Monomer> smilesMonomerDB) {
		this.smilesMonomerDB = smilesMonomerDB;
	}
}
