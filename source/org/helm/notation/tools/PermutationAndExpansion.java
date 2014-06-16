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

import org.helm.notation.tools.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class PermutationAndExpansion {

	public static List<List<String>> linearize(List<List<String[]>> originalList) {
		List<List<String>> result = new ArrayList<List<String>>();
		for (List<String[]> l : originalList) {
			List<String> tmp = new ArrayList<String>();
			for (String[] sa : l) {
				for (int i = 0; i < sa.length; i++) {
					tmp.add(sa[i]);
				}
			}
			result.add(tmp);
		}
		return result;
	}

	public static void expand(List<List<String[]>> parent, List<String[]> child)
			throws IOException, ClassNotFoundException {
		if (parent.size() == 0) {
			parent.add(child);
		} else {
			List<List<String[]>> remove = new ArrayList<List<String[]>>();
			List<List<String[]>> keep = new ArrayList<List<String[]>>();
			for (int i = 0; i < parent.size(); i++) {
				List<String[]> tmp = parent.get(i);
				remove.add(tmp);
				for (int j = 0; j < child.size(); j++) {
					List<String[]> l = DeepCopy.copy(tmp);
					l.add(child.get(j));
					keep.add(l);
				}
			}
			parent.removeAll(remove);
			parent.addAll(keep);
		}
	}

	public static void permutate(List<String[]> l, List<String> ls) {
		String[] a = ls.toArray(new String[0]);
		int n = ls.size();
		permutate(l, a, n);
	}

	public static void permutate(List<String[]> l, String[] a, int n) {

		if (n == 1) {
			String[] result = new String[a.length];
			for (int i = 0; i < a.length; i++) {
				result[i] = a[i];
			}
			l.add(result);
		}
		for (int i = 0; i < n; i++) {
			swap(a, i, n - 1);
			permutate(l, a, n - 1);
			swap(a, i, n - 1);
		}
	}

	// swap the String at indices i and j
	private static void swap(String[] a, int i, int j) {
		String c = a[i];
		a[i] = a[j];
		a[j] = c;
	}
}
