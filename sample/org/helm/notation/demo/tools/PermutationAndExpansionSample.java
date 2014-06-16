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
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * 
 * @author ZHANGTIANHONG
 */
public class PermutationAndExpansionSample {

	/**
	 * @param args
	 *            the command line arguments
	 */
	public static void main(String[] args) {

		String[] chems = new String[] { "CHEM1" };
		String[] rnas = new String[] { "RNA1", "RNA2" };
		Map<String, String[]> map = new HashMap<String, String[]>();
		map.put("PEG3", chems);
		map.put("R(A)P", rnas);

		try {

			// System.out.println("RNA12".split("\\d")[0]);
			Map<String, List<String[]>> permuatedMap = new TreeMap<String, List<String[]>>();
			Set<String> keyset = map.keySet();
			for (String key : keyset) {
				String[] value = map.get(key);
				List<String[]> l = new ArrayList<String[]>();
				PermutationAndExpansion.permutate(l, value, value.length);
				permuatedMap.put(key, l);
			}

			List<List<String[]>> expandedList = new ArrayList<List<String[]>>();
			for (String key : keyset) {
				List<String[]> value = permuatedMap.get(key);
				PermutationAndExpansion.expand(expandedList, value);
			}

			List<List<String>> linearList = PermutationAndExpansion
					.linearize(expandedList);
			for (List<String> sl : linearList) {
				for (int j = 0; j < sl.size(); j++) {
					System.out.print(sl.get(j) + "$");
				}
				System.out.println();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void explodeMap(Map map) throws IOException,
			ClassNotFoundException {
		Set keyset = map.keySet();
		for (Object key : keyset) {
			Object value = map.get(key);
			if (value instanceof List) {
				List<String> l = (List<String>) value;
				if (l.size() > 1) {
					Map m = list2map(l);
					map.put(key, m);
				} else if (l.size() == 1) {
					map.put(key, l.get(0));
				} else {
					System.out.println("list is null or empty");
				}
			} else if (value instanceof Map) {
				Map m = (Map) value;
				explodeMap((Map) value);
			} else {
				System.out.println("map is fully exploded");
			}
		}
	}

	public static TreeMap<String, List<String>> list2map(List<String> list)
			throws IOException, ClassNotFoundException {
		TreeMap<String, List<String>> map = new TreeMap<String, List<String>>();
		for (String str : list) {
			List<String> l = DeepCopy.copy(list);
			l.remove(str);
			map.put(str, l);
		}
		return map;
	}

	public static void perm1(String s) {
		perm1("", s);
	}

	private static void perm1(String prefix, String s) {
		int N = s.length();
		if (N == 0) {
			System.out.println(prefix);
		} else {
			for (int i = 0; i < N; i++) {
				perm1(prefix + s.charAt(i),
						s.substring(0, i) + s.substring(i + 1, N));
			}
		}

	}
}
