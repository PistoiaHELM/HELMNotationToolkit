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
import org.helm.notation.NotationException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.Nucleotide;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author lih25
 */
public class DeepCopy {

	public static List copy(List list) throws IOException,
			ClassNotFoundException {
		// serialize ArrayList into byte array

		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(list);
		byte buf[] = baos.toByteArray();
		oos.close();

		// deserialize byte array into ArrayList

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		ArrayList newlist = (ArrayList) ois.readObject();
		ois.close();

		return newlist;
	}

	public static ArrayList copy(ArrayList list) throws IOException,
			ClassNotFoundException {
		// serialize ArrayList into byte array

		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(list);
		byte buf[] = baos.toByteArray();
		oos.close();

		// deserialize byte array into ArrayList

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		ArrayList newlist = (ArrayList) ois.readObject();
		ois.close();

		return newlist;
	}

	public static Monomer copy(Monomer monomer) throws IOException,
			MonomerException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(monomer);
		byte buf[] = baos.toByteArray();
		oos.close();

		// deserialize byte array into ArrayList

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		Monomer newMonomer;
		try {
			newMonomer = (Monomer) ois.readObject();
		} catch (ClassNotFoundException cnfe) {
			throw new MonomerException("Unable to copy Monomer");
		}
		ois.close();

		return newMonomer;
	}

	public static Attachment copy(Attachment attachment) throws IOException,
			MonomerException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(attachment);
		byte buf[] = baos.toByteArray();
		oos.close();

		// deserialize byte array into ArrayList

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		Attachment newAttachment;
		try {
			newAttachment = (Attachment) ois.readObject();
		} catch (ClassNotFoundException cnfe) {
			throw new MonomerException("Unable to copy Attachment");
		}

		ois.close();

		return newAttachment;
	}

	public static Nucleotide copy(Nucleotide nucleotide) throws IOException,
			NotationException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(nucleotide);
		byte buf[] = baos.toByteArray();
		oos.close();

		// deserialize byte array into ArrayList

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		Nucleotide nuc;
		try {
			nuc = (Nucleotide) ois.readObject();
		} catch (ClassNotFoundException cnfe) {
			throw new NotationException("Unable to copy Nucleotide");
		}

		ois.close();

		return nuc;
	}

	public static Serializable copy(Serializable input) throws IOException,
			NotationException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(input);
		byte buf[] = baos.toByteArray();
		oos.close();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		ObjectInputStream ois = new ObjectInputStream(bais);
		Serializable output;
		try {
			output = (Serializable) ois.readObject();
		} catch (ClassNotFoundException cnfe) {
			throw new NotationException("Unable to copy input object");
		}
		ois.close();

		return output;
	}
}
