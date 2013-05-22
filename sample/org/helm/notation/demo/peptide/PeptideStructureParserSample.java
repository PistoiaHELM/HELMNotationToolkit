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
package org.helm.notation.demo.peptide;

import org.helm.notation.StructureException;
import org.helm.notation.peptide.PeptideStructureParser;
import java.io.IOException;

/**
 *
 * @author ZHANGTIANHONG
 */
public class PeptideStructureParserSample {

    public static void main(String[] args) {

        try {
            PeptideStructureParser parser = PeptideStructureParser.getInstance();
            parser.initAminoAcidLists();

            //linear peptide
            String smiles = "[H]N[C@@H](C)C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(O)=O";
            String notation = "PEPTIDE1{A.G.G.A}$$$$";
            smiles2notation(parser, smiles, notation);

            //terminal ac
            smiles = "[H]SC[C@H](NC(=O)CNC(=O)[C@H](C)NC(C)=O)C(=O)N[C@@H](C)C(O)=O";
            notation = "PEPTIDE1{[ac].A.G.C.A}$$$$";
            smiles2notation(parser, smiles, notation);

            //backbone-backbone amide cycle
            smiles = "C[C@@H]1NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](C)NC1=O";
            notation = "PEPTIDE1{A.A.E.D}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$";
            smiles2notation(parser, smiles, notation);

            //backbone amine - branch carbonul amide cycle
            smiles = "C[C@@H]1NC(=O)C[C@H](NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](C)NC1=O)C(O)=O";
            notation = "PEPTIDE1{A.A.E.D}$PEPTIDE1,PEPTIDE1,4:R3-1:R1$$$";
            smiles2notation(parser, smiles, notation);


            //backbone-branch amide cycle
            smiles = "C[C@H](NC(=O)[C@@H]1CC(=O)N[C@@H](C)C(=O)N[C@@H](CCC(O)=O)C(=O)N1)C(O)=O";
            notation = "PEPTIDE1{A.E.D.A}$PEPTIDE1,PEPTIDE1,1:R1-3:R3$$$";
            smiles2notation(parser, smiles, notation);

            //branch-branch amide cycle
            smiles = "[H]N[C@@H](C)C(=O)N[C@H]1CCCCNC(=O)C[C@H](NC(=O)CNC(=O)[C@H](CC(C)C)NC1=O)C(=O)N[C@@H](C)C(O)=O";
            notation = "PEPTIDE1{A.K.L.G.D.A}$PEPTIDE1,PEPTIDE1,5:R3-2:R3$$$";
            smiles2notation(parser, smiles, notation);

            //branch peptides
            smiles = "[H]N[C@@H](C)C(=O)N[C@@H](CS[H])C(=O)N[C@@H](CC(=O)N[C@@H](CCSC)C(=O)N[C@@H](CC(N)=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(N)=O)C(O)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)NCC(O)=O";
            notation = "PEPTIDE1{A.C.D.E.F.G}|PEPTIDE2{M.N.P.Q}$PEPTIDE1,PEPTIDE2,3:R3-1:R1$$$";
            smiles2notation(parser, smiles, notation);

            //disulfide cycle
            smiles = "[H]N[C@@H](C)C(=O)N[C@H]1CSSC[C@H](NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(O)=O)NC1=O)C(=O)NCC(O)=O";
            notation = "PEPTIDE1{A.C.D.E.C.G}$PEPTIDE1,PEPTIDE1,5:R3-2:R3$$$";
            smiles2notation(parser, smiles, notation);

            //disulfide branch
            smiles = "[H]N[C@@H](C)C(=O)N[C@@H](CSSC[C@H](NC(=O)[C@@H](N[H])CCC(O)=O)C(=O)NCC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O";
            notation = "PEPTIDE1{A.C.D}|PEPTIDE2{E.C.G}$PEPTIDE2,PEPTIDE1,2:R3-2:R3$$$";
            smiles2notation(parser, smiles, notation);

            //modified aa
            smiles = "[H]NCCCC[C@H](NC(=O)[C@@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@H](CC1=CNC=N1)NC(=O)CNC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CS[H])NC(=O)[C@H](C)N[H])[C@@H](C)CC)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CC(N)=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](CO)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC1=CNC2=C1C=CC=C2)C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(=O)N[C@@H](C[SeH])C(=O)NCC(O)=O";
            notation = "PEPTIDE1{A.C.D.E.F.G.H.I.[dK].K.L.M.N.P.Q.R.S.T.V.W.Y.[seC].G}$$$$";
            smiles2notation(parser, smiles, notation);


            //mixture            
            smiles = "[H]N[C@@H](C)C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(O)=O.[H]N[C@@H](C)C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(O)=O";
            notation = "PEPTIDE1{A.G.G.A}|PEPTIDE1{A.G.G.A}$$$$";
            smiles2notation(parser, smiles, notation);

            //PF-05297218
            smiles = "CCc1cc(OC)ccc1-c1ccc(C[C@H](NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@](C)(Cc2c(F)cccc2F)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)C(C)(C)NC(=O)[C@@H](N)Cc2cnc[nH]2)[C@@H](C)O)[C@@H](C)O)C(=O)N[C@@H](Cc2ccc(cc2)-c2ccccc2C)C(N)=O)cc1";
            smiles2notation(parser, smiles, notation);

        } catch (Exception e) {
            e.printStackTrace();
        }
        System.exit(0);
    }

    private static void smiles2notation(PeptideStructureParser parser, String smiles, String notation) {
        try {
            System.out.println("Input    smiles:\t" + smiles);
            System.out.println("Input  notation:\t" + notation);
            String result = parser.smiles2notation(smiles);
            System.out.println("Output notation:\t" + result);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
