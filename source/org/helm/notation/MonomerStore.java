package org.helm.notation;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.helm.notation.model.Monomer;
import org.helm.notation.tools.DeepCopy;

public class MonomerStore {
	private Map<String, Map<String, Monomer>> monomerDB;
    private Map<String, Monomer> smilesMonomerDB;
    
    public MonomerStore( Map<String, Map<String, Monomer>> monomerDB, Map<String, Monomer> smilesMonomerDB) {
    	this.monomerDB = monomerDB;
    	this.smilesMonomerDB = smilesMonomerDB;
    }
    
    public MonomerStore() {
    	monomerDB = new HashMap<String, Map<String, Monomer>>();
    	smilesMonomerDB = new HashMap<String, Monomer>();
    }

    public Map<String, Map<String, Monomer>> getMonomerDB() {
    	return monomerDB;
    }
    
    public  Map<String, Monomer> getSmilesMonomerDB() {
    	return smilesMonomerDB;
    }

    public void addMonomer( Monomer monomer) throws IOException, MonomerException {
	 Map<String, Monomer> monomerMap = monomerDB.get(monomer.getPolymerType());
	 String polymerType = monomer.getPolymerType();
 	 String alternateId = monomer.getAlternateId();
 	 String smilesString = monomer.getCanSMILES();
 	 
 	 boolean hasSmilesString =  (smilesString != null && smilesString.length() > 0);
 	 
	 if (null == monomerMap) {
		 monomerMap = new HashMap<String, Monomer>();         
         monomerDB.put(polymerType, monomerMap);
     } 
     
 	 
	 Monomer copyMonomer = DeepCopy.copy(monomer);

	 
	 boolean alreadyInMonomerMap = monomerMap.containsKey(alternateId);
	 if (!alreadyInMonomerMap){ 
		 monomerMap.put(alternateId, copyMonomer);
	 }
	 
	 
	 boolean alreadyInSMILESMap = hasSmilesString && (smilesMonomerDB.containsKey(smilesString));

	 if (!alreadyInSMILESMap) {
    	 smilesMonomerDB.put(smilesString, monomer);
     }
  /*   if (!alreadyInMonomerMap
    		 && !alreadyInSMILESMap) {
         monomerMap.put(alternateId, copyMonomer);
         if ( hasSmilesString) {
        	 smilesMonomerDB.put(monomer.getCanSMILES(), monomer);
         }
     }
     else {
    	 throw new MonomerException( "Monomer already exists in store " + monomer.getAlternateId());
     }*/
   
    }
    
    public Monomer getMonomer( String polymerType, String alternateId) {
    	return monomerDB.get( polymerType).get(alternateId);
    }
    
    public Monomer getMonomer( String smiles) {
    	return smilesMonomerDB.get( smiles);
    }
    
    public Map<String, Monomer> getMonomers( String polymerType) {
    	return monomerDB.get( polymerType);
    }
    
    public synchronized void addNewMonomer(Monomer monomer) throws IOException, MonomerException {
        monomer.setNewMonomer(true);
        addMonomer(monomer);
    }

}

