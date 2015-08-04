package org.helm.notation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.helm.notation.model.Monomer;
import org.helm.notation.tools.DeepCopy;
import org.helm.notation.tools.StructureParser;

/**
 * This class represents a store for monomers. It is mainly used to seperate the
 * monomers coming from a single (XHELM) file from the monomers within the local
 * database.
 * 
 * @author maisel
 * 
 */
public class MonomerStore {
	private Map<String, Map<String, Monomer>> monomerDB;
	private Map<String, Monomer> smilesMonomerDB;

	/**
	 * Constructor with Monomer- and SmilesDB
	 * 
	 * @param monomerDB
	 * @param smilesMonomerDB
	 */
	public MonomerStore(Map<String, Map<String, Monomer>> monomerDB,
			Map<String, Monomer> smilesMonomerDB) {
		this.monomerDB = monomerDB;
		this.smilesMonomerDB = smilesMonomerDB;
	}

	/**
	 * 
	 * Constructs empty MonomerStore
	 * 
	 */
	public MonomerStore() {
		monomerDB = new HashMap<String, Map<String, Monomer>>();
		smilesMonomerDB = new HashMap<String, Monomer>();
	}

	/**
	 * returns MonomerDB
	 * 
	 * @return MonomerDB as Map<String, Map<String, Monomer>>
	 */
	public Map<String, Map<String, Monomer>> getMonomerDB() {
		return monomerDB;
	}

	/**
	 * returns SmilesMonomerDB
	 * 
	 * @return SmilesMonomerDB as Map<String, Monomer>
	 */
	public Map<String, Monomer> getSmilesMonomerDB() {
		return smilesMonomerDB;
	}

	/**
	 * Adds a monomer to the store
	 * 
	 * @param monomer
	 * @throws IOException
	 * @throws MonomerException
	 */
	public void addMonomer(Monomer monomer) throws IOException,
			MonomerException {
		addMonomer(monomer, false);
	}

	/**
	 * Adds a monomer to the store and optionally sets the dbChanged flag
	 * 
	 * @param monomer
	 * @param dbChanged
	 * @throws IOException
	 * @throws MonomerException
	 */
	public void addMonomer(Monomer monomer, boolean dbChanged)
			throws IOException, MonomerException {
		Map<String, Monomer> monomerMap = monomerDB.get(monomer
				.getPolymerType());
		String polymerType = monomer.getPolymerType();
		String alternateId = monomer.getAlternateId();
		String smilesString = monomer.getCanSMILES();

		try {
			smilesString = StructureParser
					.getUniqueExtendedSMILES(smilesString);
		} catch (Exception e) {
			smilesString = monomer.getCanSMILES();
		}

		boolean hasSmilesString = (smilesString != null && smilesString
				.length() > 0);

		if (null == monomerMap) {
			monomerMap = new HashMap<String, Monomer>();
			monomerDB.put(polymerType, monomerMap);
		}

		Monomer copyMonomer = DeepCopy.copy(monomer);
		
		// ensure the canonical SMILES is indexed in the monomer store
		if (hasSmilesString) {
			copyMonomer.setCanSMILES(smilesString);
		}

		boolean alreadyAdded = false;
		alreadyAdded = monomerMap.containsKey(alternateId);

		if (!alreadyAdded) {
			monomerMap.put(alternateId, copyMonomer);

			boolean alreadyInSMILESMap = hasSmilesString
					&& (smilesMonomerDB.containsKey(smilesString));

			if (!alreadyInSMILESMap) {
				smilesMonomerDB.put(smilesString, copyMonomer);
			}
		}

		if (dbChanged) {
			MonomerFactory.setDBChanged(true);
		}
	}

	/**
	 * Checks if a specific monomer exists in the store
	 * 
	 * @param polymerType
	 * @param alternateId
	 * @return true if monomer exists, false if not
	 */
	public boolean hasMonomer(String polymerType, String alternateId) {
		return ((monomerDB.get(polymerType) != null) && getMonomer(polymerType,
				alternateId) != null);
	}

	/**
	 * Returns the monomer specified by polymerType and alternatId
	 * 
	 * @param polymerType
	 * @param alternateId
	 * @return the matching monomer
	 */
	public Monomer getMonomer(String polymerType, String alternateId) {
		return monomerDB.get(polymerType).get(alternateId);
	}

	/**
	 * Returns the monomer by smiles string
	 * 
	 * @param smiles
	 * @return the matching monomer
	 */
	public Monomer getMonomer(String smiles) {
		return smilesMonomerDB.get(smiles);
	}

	/**
	 * Returns all monomers by polymerType
	 * 
	 * @param polymerType
	 * @return All monomers with polymerType
	 */
	public Map<String, Monomer> getMonomers(String polymerType) {
		return monomerDB.get(polymerType);
	}

	/**
	 * Adds a monomer to the store and makes it a temporary new monomer
	 * 
	 * @param monomer
	 * @throws IOException
	 * @throws MonomerException
	 */
	public synchronized void addNewMonomer(Monomer monomer) throws IOException,
			MonomerException {
		monomer.setNewMonomer(true);
		addMonomer(monomer, true);
	}

	/**
	 * Checks for the empty store
	 * 
	 * @return true if the store is empty, false if not
	 */
	public boolean isMonomerStoreEmpty() {
		return (this.monomerDB == null || this.monomerDB.values() == null || this.monomerDB
				.values().size() == 0);
	}

	/**
	 * Clears the MonomerStore
	 */
	public synchronized void clearMonomers() {
		this.monomerDB.clear();
		this.smilesMonomerDB.clear();
	}

	public String toString() {
		String str = "";
		for (Map<String, Monomer> val : this.monomerDB.values()) {
			for (Monomer mon : val.values()) {
				str += mon.getAlternateId() + "(" + mon.getPolymerType()
						+ "); ";
			}
			str += System.getProperty("line.separator");
		}

		return str;
	}

	/**
	 * Returns the polymer type set
	 * 
	 * @return the polymer type set as Set<String>
	 */
	public Set<String> getPolymerTypeSet() {
		return monomerDB.keySet();
	}

	/**
	 * This method returns all monomers of the store as list sorted by polymer
	 * type
	 * 
	 * @return all monomers of store as List<Monomer>
	 */
	public List<Monomer> getAllMonomersList() {
		List<Monomer> monomers = new ArrayList<Monomer>();
		for (String polymerType : getPolymerTypeSet()) {
			Map<String, Monomer> map = getMonomers(polymerType);
			monomers.addAll(map.values());
		}
		return monomers;

	}
}
