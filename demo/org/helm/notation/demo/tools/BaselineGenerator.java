package org.helm.notation.demo.tools;

import org.helm.notation.model.MoleculeInfo;
import org.helm.notation.model.Monomer;
import org.helm.notation.tools.ComplexNotationParser;
import org.helm.notation.tools.NucleotideSequenceParser;
import org.helm.notation.tools.SimpleNotationParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class BaselineGenerator {

    private static class CaseResult {
        private String id;
        private String operation;
        private String input;
        private String output;
        private String error;
    }

    public static void main(String[] args) throws Exception {
        String outPath = args.length > 0 ? args[0] : "migration/fixtures/java_baseline.json";
        List<CaseResult> results = new ArrayList<CaseResult>();

        addSimpleSmilesCase(results,
                "simple_peptide_smiles_1",
                "A.G.G.G.C.C.K.K.K.K",
                Monomer.PEPTIDE_POLYMER_TYPE);

        addSimpleSmilesCase(results,
                "simple_rna_smiles_1",
                "R(A)P.R(U)P.R(G)",
                Monomer.NUCLIEC_ACID_POLYMER_TYPE);

        addSequenceCase(results,
                "sequence_to_notation_1",
                "AUGCUU");

        String complexNotation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
        addComplexCanonicalCase(results, "complex_canonical_1", complexNotation);
        addComplexSmilesCase(results, "complex_smiles_1", complexNotation);
        addComplexMoleculeInfoCase(results, "complex_molinfo_1", complexNotation);

        String cyclicPeptide = "PEPTIDE1{K.A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$";
        addComplexCanonicalCase(results, "complex_canonical_2", cyclicPeptide);

        writeJson(outPath, results);
        System.out.println("Wrote baseline fixtures to " + outPath + " with " + results.size() + " entries.");
    }

    private static void addSimpleSmilesCase(List<CaseResult> results, String id, String notation, String polymerType) {
        CaseResult result = new CaseResult();
        result.id = id;
        result.operation = "simple_polymer_smiles";
        result.input = notation + "||" + polymerType;
        try {
            result.output = SimpleNotationParser.getSimplePolymerSMILES(notation, polymerType);
        } catch (Exception ex) {
            result.error = ex.getClass().getSimpleName() + ": " + ex.getMessage();
        }
        results.add(result);
    }

    private static void addSequenceCase(List<CaseResult> results, String id, String sequence) {
        CaseResult result = new CaseResult();
        result.id = id;
        result.operation = "sequence_to_notation";
        result.input = sequence;
        try {
            result.output = NucleotideSequenceParser.getNotation(sequence);
        } catch (Exception ex) {
            result.error = ex.getClass().getSimpleName() + ": " + ex.getMessage();
        }
        results.add(result);
    }

    private static void addComplexCanonicalCase(List<CaseResult> results, String id, String notation) {
        CaseResult result = new CaseResult();
        result.id = id;
        result.operation = "complex_canonical_notation";
        result.input = notation;
        try {
            result.output = ComplexNotationParser.getCanonicalNotation(notation);
        } catch (Exception ex) {
            result.error = ex.getClass().getSimpleName() + ": " + ex.getMessage();
        }
        results.add(result);
    }

    private static void addComplexSmilesCase(List<CaseResult> results, String id, String notation) {
        CaseResult result = new CaseResult();
        result.id = id;
        result.operation = "complex_polymer_smiles";
        result.input = notation;
        try {
            result.output = ComplexNotationParser.getComplexPolymerSMILES(notation);
        } catch (Exception ex) {
            result.error = ex.getClass().getSimpleName() + ": " + ex.getMessage();
        }
        results.add(result);
    }

    private static void addComplexMoleculeInfoCase(List<CaseResult> results, String id, String notation) {
        CaseResult result = new CaseResult();
        result.id = id;
        result.operation = "complex_molecule_info";
        result.input = notation;
        try {
            MoleculeInfo info = ComplexNotationParser.getMoleculeInfo(notation);
            result.output = "MW=" + info.getMolecularWeight()
                    + ";EXACT=" + info.getExactMass()
                    + ";FORMULA=" + info.getMolecularFormula();
        } catch (Exception ex) {
            result.error = ex.getClass().getSimpleName() + ": " + ex.getMessage();
        }
        results.add(result);
    }

    private static void writeJson(String outPath, List<CaseResult> results) throws Exception {
        File file = new File(outPath);
        File parent = file.getParentFile();
        if (parent != null && !parent.exists()) {
            parent.mkdirs();
        }

        BufferedWriter writer = new BufferedWriter(new FileWriter(file));
        writer.write("{\n");
        writer.write("  \"source\": \"java\",\n");
        writer.write("  \"generator\": \"org.helm.notation.demo.tools.BaselineGenerator\",\n");
        writer.write("  \"cases\": [\n");

        for (int i = 0; i < results.size(); i++) {
            CaseResult result = results.get(i);
            writer.write("    {\n");
            writer.write("      \"id\": \"" + escape(result.id) + "\",\n");
            writer.write("      \"operation\": \"" + escape(result.operation) + "\",\n");
            writer.write("      \"input\": \"" + escape(result.input) + "\",\n");
            writer.write("      \"output\": " + toJsonNullableString(result.output) + ",\n");
            writer.write("      \"error\": " + toJsonNullableString(result.error) + "\n");
            writer.write("    }");
            if (i < results.size() - 1) {
                writer.write(",");
            }
            writer.write("\n");
        }

        writer.write("  ]\n");
        writer.write("}\n");
        writer.close();
    }

    private static String toJsonNullableString(String value) {
        if (value == null) {
            return "null";
        }
        return "\"" + escape(value) + "\"";
    }

    private static String escape(String input) {
        if (input == null) {
            return "";
        }
        return input
                .replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("\n", "\\n")
                .replace("\r", "\\r")
                .replace("\t", "\\t");
    }
}
