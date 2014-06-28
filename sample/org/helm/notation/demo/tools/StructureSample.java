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
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import org.helm.notation.model.MoleculeInfo;

/**
 * This is an example on how to use StructureParser
 * 
 * @author ZHANGTIANHONG
 */
public class StructureSample {

	public static void main(String[] args) {

		try {

			String ethyl = "CC[*] |$;;_R1$|";
			String benzene = "Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|";

			Molecule ethylMol = StructureParser.getMolecule(ethyl);
			Molecule benzeneMol = StructureParser.getMolecule(benzene);
			// benzeneMol.dearomatize();
			System.out.println(benzeneMol.toFormat("smiles"));

			MolAtom ethylA = StructureParser.getRgroupAtom(ethylMol, 1);
			MolAtom benzeneA = StructureParser.getRgroupAtom(benzeneMol, 1);
			;

			StructureParser.merge(benzeneMol, benzeneA, ethylMol, ethylA);
			System.out.println(benzeneMol.toFormat("smiles:u"));

			String testSmiles = "[H]OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)[nH]1cnc2c1nc(N)[nH]c2=O)[nH]1ccc(=O)[nH]c1=O)[nH]1cnc2c(N)ncnc12";
			boolean isValid = StructureParser.validateSmiles(testSmiles);

			// SMILES validation
			String validSmiles = "CC(C)[C@H](N[*])C([*])=O |r,$;;;;;_R1;;_R2;$|";
			String invalidSmiles = "Cc7ccccc7Cc1(c7ccccc7)nnnc2cc(C(=O))c21";

			if (StructureParser.validateSmiles(validSmiles)) {
				System.out.println("Valid SMILES is valid");
			} else {
				System.out.println("Valid SMILES is invalid");
			}
			if (StructureParser.validateSmiles(invalidSmiles)) {
				System.out.println("Invalid SMILES is valid");
			} else {
				System.out.println("invalid SMILES is invalid");
			}

			MoleculeInfo mi = StructureParser.getMoleculeInfo(validSmiles);
			System.out.println("MF: " + mi.getMolecularFormula());
			System.out.println("MW: " + mi.getMolecularWeight());
			System.out.println("Exact Mass: " + mi.getExactMass());

			String aromatic = "[H]O[C@@H]1[C@@H](COP(S)(=O)O[C@@H]2[C@@H](COP(O)(=O)O[C@@H]3[C@@H](COP(O)(=O)O[C@@H]4[C@@H](COP(O)(=O)O[C@@H]5[C@@H](COP(O)(=O)O[C@@H]6[C@@H](COP(O)(=O)O[C@@H]7[C@@H](COP(O)(=O)O[C@@H]8[C@@H](COP(O)(=O)O[C@@H]9[C@@H](COP(O)(=O)O[C@@H]%10[C@@H](COP(O)(=O)O[C@@H]%11[C@@H](COP(O)(=O)O[C@@H]%12[C@@H](COP(O)(=O)O[C@@H]%13[C@@H](COP(O)(=O)O[C@@H]%14[C@@H](COP(O)(=O)O[C@@H]%15[C@@H](COP(O)(=O)O[C@@H]%16[C@@H](COP(O)(=O)O[C@@H]%17[C@@H](COP(O)(=O)O[C@@H]%18[C@@H](COP(O)(=O)O[C@@H]%19[C@@H](COP(O)(=O)O[C@@H]%20[C@@H](COP(O)(=O)O[C@@H]%21[C@@H](COP(S)(=O)OCCCCC#C[H])O[C@@H]([C@@H]%21O)n%21cnc%22c%21nc(N)[nH]c%22=O)O[C@@H]([C@@H]%20O)n%20ccc(=O)[nH]c%20=O)O[C@@H]([C@@H]%19O)n%19ccc(N)nc%19=O)O[C@@H]([C@@H]%18O)n%18cnc%19c(N)ncnc%18%19)O[C@@H]([C@@H]%17O)n%17ccc(=O)[nH]c%17=O)O[C@@H]([C@@H]%16O)n%16ccc(N)nc%16=O)O[C@@H]([C@@H]%15O)n%15cnc%16c(N)ncnc%15%16)O[C@@H]([C@@H]%14O)n%14ccc(N)nc%14=O)O[C@@H]([C@@H]%13O)n%13cnc%14c(N)ncnc%13%14)O[C@@H]([C@@H]%12O)n%12ccc(N)nc%12=O)O[C@@H]([C@@H]%11O)n%11ccc(=O)[nH]c%11=O)O[C@@H]([C@@H]%10O)n%10cnc%11c%10nc(N)[nH]c%11=O)O[C@@H]([C@@H]9O)n9cnc%10c(N)ncnc9%10)O[C@@H]([C@@H]8O)n8cnc9c(N)ncnc89)O[C@@H]([C@@H]7O)n7ccc(=O)[nH]c7=O)O[C@@H]([C@@H]6O)n6cnc7c(N)ncnc67)O[C@@H]([C@@H]5O)n5ccc(N)nc5=O)O[C@@H]([C@@H]4O)n4ccc(N)nc4=O)O[C@@H]([C@@H]3O)n3cnc4c(N)ncnc34)O[C@@H]([C@@H]2O)n2cnc3c(N)ncnc23)O[C@@H]([C@@H]1O)n1ccc(=O)[nH]c1=O";
			Molecule mol = StructureParser.getMolecule(aromatic);
			mol.dearomatize();
			String dearomatic = mol.toFormat("smiles");
			System.out.println(dearomatic);

			System.out.println("Canonicalization .................");
			String inSmi = null;
			String canSmi = null;
			String extSmi = null;

			inSmi = "[H]N[C@@H](CCCCN)C(=O)N[C@@H](CCCNC(N)=N)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(O)=O";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			inSmi = "[H]NCCCC[C@H](N[H])C(=O)N[C@@H](CCCNC(N)=N)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(O)=O";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			inSmi = "[*]N[C@@H](CCCCN[*])C([*])=O |r,$_R1;;;;;;;;_R3;;_R2;$|";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			// contains explicit hydrogen
			inSmi = "[H]NCCCC[C@H](N[*])C([*])=O |r,$;;;;;;;;_R1;;_R2;$|";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			// contains explicit hydrogen and huckle ring
			inSmi = "[H]NC(CCC[C@H](N[*])C([*])=O)C1=CC=CC=C1 |r,$;;;;;;;;_R1;;_R2;;;;;;;$,c:14,16,t:12|";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			inSmi = "[*]NCC([*])=O |$_R1;;;;_R2;$|";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			inSmi = "[*]NCC([*])=O";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			inSmi = "[*]C(=O)CN[*]";
			System.out.println("Input Smiles:\t\t" + inSmi);
			extSmi = StructureParser.getUniqueExtendedSMILES(inSmi);
			System.out.println("Extended Smiles:\t" + extSmi);
			canSmi = StructureParser.getUniqueSmiles(inSmi);
			System.out.println("Cannonical Smiles:\t" + canSmi);

			// arg-ser
			inSmi = "N[C@H](CCCNC(N)=N)C(=O)N[C@@H](CO)C=O";
			System.out.println("Input Smiles:\t\t" + inSmi);
			String pdb = StructureParser.SMILES2ChemAxonPDB(inSmi);
			System.out.println("ChemAxon PDB:\t" + pdb);

		} catch (Exception e) {
			e.printStackTrace();
		}
		System.exit(0);
	}
}
