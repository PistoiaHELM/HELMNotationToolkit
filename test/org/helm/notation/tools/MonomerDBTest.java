package org.helm.notation.tools;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.helm.notation.MonomerException;
import org.helm.notation.MonomerFactory;
import org.helm.notation.StructureException;
import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import org.jdom.JDOMException;
import org.junit.Test;

public class MonomerDBTest {

	@Test
    public void initializeMonomerDB() throws IOException, MonomerException, JDOMException, StructureException{
		
        //show all monomers
        Map monomerDB = MonomerFactory.getInstance().getMonomerDB();
        //showMonomerDB(monomerDB);

        //show all attachment
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


}

public void showMonomerDB(Map monomerDB) throws StructureException, IOException {
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



