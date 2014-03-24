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
package org.helm.notation;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.helm.notation.model.Attachment;
import org.helm.notation.model.Monomer;
import org.helm.notation.model.MonomerCache;
import org.helm.notation.tools.DeepCopy;
import org.helm.notation.tools.MonomerParser;
import org.jdom.Attribute;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

/**
 * This is a factory class to build monomer database from
 * MonomerDBGZEnconded.xml document
 * 
 * @author zhangtianhong
 */
public class MonomerFactory {

	public static final String NOTATION_DIRECTORY = NotationConstant.NOTATION_DIRECTORY;
	public static final String MONOMER_CACHE_FILE_NAME = "MonomerCache.ser";
	public static final String MONOMER_CACHE_FILE_PATH = NOTATION_DIRECTORY
			+ System.getProperty("file.separator") + MONOMER_CACHE_FILE_NAME;
	public static final String MONOMER_DB_FILE_NAME = "MonomerDBGZEncoded.xml";
	public static final String MONOMER_DB_FILE_PATH = NOTATION_DIRECTORY
			+ System.getProperty("file.separator") + MONOMER_DB_FILE_NAME;
	public static final String MONOMER_DB_XML_RESOURCE = "resources/MonomerDBGZEncoded.xml";
	public static final String MONOMER_DB_SCHEMA_RESOURCE = "resources/MonomerDBSchema.xsd";
	public static final String XML_SCHEMA_VALIDATION_FEATURE = "http://apache.org/xml/features/validation/schema";
	public static final String EXTERNAL_SCHEMA_LOCATION_KEY = "http://apache.org/xml/properties/schema/external-schemaLocation";
	public static final String DEFAULT_NAME_SPACE = "lmr";
	public static final String POLYMER_LIST_ELEMENT = "PolymerList";
	public static final String POLYMER_ELEMENT = "Polymer";
	public static final String POLYMER_TYPE_ATTRIBUTE = "polymerType";
	public static final String ATTACHMENT_LIST_ELEMENT = "AttachmentList";
	private static MonomerFactory instance;
	/**
	 * First key is polymer Type, such as "RNA" Second key is monomer ID, such
	 * as "A"
	 */
	private static Map<String, Map<String, Monomer>> monomerDB; // key is
																// monomer
																// SMILES, value
																// is Monomer
	private static Map<String, Monomer> smilesMonomerDB; // key is
															// AttachementID,
															// value is
															// Attachment
	private static Map<String, Attachment> attachmentDB;
	// private static Map<String, Map<String, Monomer>> externalMonomerDB;
	private static SAXBuilder builder;
	private static Logger logger = Logger.getLogger(MonomerFactory.class
			.toString());

	private static boolean dbChanged = true;

	/**
	 * retruns the monomer database
	 * 
	 * @return Map as Map<String, Map<String, Monomer>>
	 */
	public synchronized Map<String, Map<String, Monomer>> getMonomerDB() {
		return monomerDB;
	}

	/*
	 * public synchronized Map<String, Map<String, Monomer>>
	 * getExternalMonomerDB() { if ( externalMonomerDB == null) {
	 * externalMonomerDB = new HashMap<String, Map<String, Monomer>>(); } return
	 * externalMonomerDB; } public synchronized void
	 * setExternalMonomerDB(Map<String, Map<String, Monomer>> map) {
	 * externalMonomerDB=map; } public void clearExternalMonomerDB(){
	 * externalMonomerDB=null; }
	 */

	protected MonomerStore monomerStore;

	public synchronized MonomerStore getMonomerStore() {
		if (monomerStore == null) {
			monomerStore = new MonomerStore(monomerDB, smilesMonomerDB);
		}
		return monomerStore;
	}

	public synchronized Map<String, Attachment> getAttachmentDB() {
		return attachmentDB;
	}

	public synchronized Map<String, Monomer> getSmilesMonomerDB() {
		return smilesMonomerDB;
	}

	public synchronized List<String> getPolymerTypes() {
		List<String> l = new ArrayList<String>();
		l.addAll(monomerDB.keySet());
		Collections.sort(l);
		return l;
	}

	public synchronized List<String> getMonomerTypes() {
		List<String> monomerTypeList = new ArrayList<String>();
		Object[] col = monomerDB.values().toArray();
		for (int i = 0; i < col.length; i++) {
			Map<String, Monomer> map = (Map<String, Monomer>) col[i];
			Monomer[] monomers = map.values().toArray(new Monomer[0]);
			for (int j = 0; j < monomers.length; j++) {
				if (!monomerTypeList.contains(monomers[j].getMonomerType())) {
					monomerTypeList.add(monomers[j].getMonomerType());
				}
			}
		}
		Collections.sort(monomerTypeList);
		return monomerTypeList;
	}

	public synchronized Map<String, List<String>> getAttachmentLabelIDs() {
		Map<String, List<String>> labelMap = new HashMap<String, List<String>>();

		// group attachments based on R value (label)
		Set<String> idSet = attachmentDB.keySet();
		for (String id : idSet) {
			String label = attachmentDB.get(id).getLabel();
			List<String> ids = labelMap.get(label);
			if (null == ids || ids.isEmpty()) {
				ids = new ArrayList<String>();
				ids.add(id);
				labelMap.put(label, ids);
			} else {
				ids.add(id);
			}
		}

		// Sort ids in each R group
		Set<String> labelSet = labelMap.keySet();
		for (String label : labelSet) {
			List<String> ids = labelMap.get(label);
			Collections.sort(ids);
		}

		return labelMap;
	}

	public static void setupBuilder() {
		URL schema = MonomerFactory.class
				.getResource(MONOMER_DB_SCHEMA_RESOURCE);
		builder = new SAXBuilder(true); // checks both well-formedness and
										// validity
		builder.setFeature(XML_SCHEMA_VALIDATION_FEATURE, true);
		builder.setProperty(EXTERNAL_SCHEMA_LOCATION_KEY, DEFAULT_NAME_SPACE
				+ " " + schema.toString());
	}

	private MonomerFactory() {
	}

	/**
	 * Initialize MonomerCache and returns the singlton Factory class
	 * 
	 * @return MonomerFactory
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public static MonomerFactory getInstance() throws MonomerException,
			IOException, JDOMException {
		if (null == instance) {
			initializeMonomerCache();
			instance = new MonomerFactory();
		}
		return instance;
	}

	public static void setDBChanged(boolean isChanged) {
	    dbChanged = isChanged;	
	}
	
	/**
	 * Returns whether one of the stored databases has changed, for example by
	 * adding or removing monomers.
	 * 
	 * @return true when database has changed else false
	 */
	public static boolean hasDBChanged() {
		return dbChanged;
	}
	public static void resetDBChanged() {
		dbChanged = false;
	}

	private static void serializeMonomerCache(MonomerCache monomerCache,
			String fileName) throws IOException {
		FileOutputStream fos = new FileOutputStream(fileName);
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(monomerCache);
		oos.close();
		fos.close();
	}

	private static MonomerCache deserializeMonomerCache(String fileName)
			throws IOException, MonomerException {
		FileInputStream fis = new FileInputStream(fileName);
		ObjectInputStream ois = new ObjectInputStream(fis);
		MonomerCache cache = null;
		try {
			cache = (MonomerCache) ois.readObject();
		} catch (ClassNotFoundException cnfe) {
			throw new MonomerException(
					"Unable to deserialize monomer cache from file");
		} finally {
			ois.close();
		}
		fis.close();
		return cache;
	}

	/**
	 * To add new monomer into monomerCache
	 * 
	 * @param monomer
	 */
	public synchronized void addNewMonomer(Monomer monomer) throws IOException,
			MonomerException {
		monomer.setNewMonomer(true);
		addMonomer(monomerDB, smilesMonomerDB, monomer);

		dbChanged = true;
	}

	private void addMonomer(Map<String, Map<String, Monomer>> monomerDB,
			Map<String, Monomer> smilesMonomerDB, Monomer monomer)
			throws IOException, MonomerException {
		Map<String, Monomer> monomerMap = monomerDB.get(monomer
				.getPolymerType());

		if (null == monomerMap) {
			Map<String, Monomer> map = new HashMap<String, Monomer>();
			Monomer copyMonomer = DeepCopy.copy(monomer);
			map.put(monomer.getAlternateId(), copyMonomer);
			monomerDB.put(monomer.getPolymerType(), map);
		} else {
			if (!monomerMap.containsKey(monomer.getAlternateId())) {
				monomerMap.put(monomer.getAlternateId(), monomer);
			}
		}

		if (monomer.getCanSMILES() != null
				&& monomer.getCanSMILES().length() > 0) {
			if (!smilesMonomerDB.containsKey(monomer.getCanSMILES())) {
				smilesMonomerDB.put(monomer.getCanSMILES(), monomer);
			}
		}

		dbChanged = true;
	}

	/*
	 * public void addExternalMonomer( Monomer monomer) throws IOException,
	 * MonomerException { String polymerType = monomer.getPolymerType(); String
	 * alternateId = monomer.getAlternateId();
	 * 
	 * Map<String, Monomer> monomerMap =
	 * getExternalMonomerDB().get(polymerType);
	 * 
	 * Monomer copyMonomer = DeepCopy.copy(monomer);
	 * 
	 * if (null == monomerMap) { monomerMap = new HashMap<String, Monomer>();
	 * externalMonomerDB.put(polymerType, monomerMap); }
	 * 
	 * if (!monomerMap.containsKey(alternateId)) { monomerMap.put(alternateId,
	 * copyMonomer); } //TODO if (monomer.getCanSMILES() != null &&
	 * monomer.getCanSMILES().length() > 0) { if
	 * (!smilesMonomerDB.containsKey(monomer.getCanSMILES())) {
	 * smilesMonomerDB.put(monomer.getCanSMILES(), monomer); } } }
	 */

	/**
	 * Build an MonomerCache object with monomerDBXML String
	 * 
	 * @param monomerDBXML
	 * @return MonomerCache
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	public MonomerCache buildMonomerCacheFromXML(String monomerDBXML)
			throws MonomerException, IOException, JDOMException,
			StructureException {
		ByteArrayInputStream bais = new ByteArrayInputStream(
				monomerDBXML.getBytes());
		return buildMonomerCacheFromXML(bais);
	}

	/**
	 * merge remote monomerCache with local monomerCache, will throw exception
	 * if conflicts found. Client needs to resolve conflicts prior to calling
	 * merge
	 * 
	 * @param remoteMonomerCache
	 * @throws java.io.IOException
	 * @throws org.helm.notation.MonomerException
	 */
	public synchronized void merge(MonomerCache remoteMonomerCache)
			throws IOException, MonomerException {
		Map<Monomer, Monomer> conflicts = getConflictedMonomerMap(remoteMonomerCache);
		if (conflicts.size() > 0) {
			throw new MonomerException(
					"Local new monomer and remote monomer database conflict found");
		} else {
			Map<String, Map<String, Monomer>> monoDB = remoteMonomerCache
					.getMonomerDB();

			Set<String> polymerTypeSet = monoDB.keySet();
			for (Iterator i = polymerTypeSet.iterator(); i.hasNext();) {
				String polymerType = (String) i.next();
				Map<String, Monomer> map = monoDB.get(polymerType);
				Set<String> monomerSet = map.keySet();
				for (Iterator it = monomerSet.iterator(); it.hasNext();) {
					String id = (String) it.next();
					Monomer m = map.get(id);
					addMonomer(monomerDB, smilesMonomerDB, m);
				}
			}
		}

		dbChanged = true;
	}

	/**
	 * replace local cache with remote one completely, may cause loss of data
	 * 
	 * @param remoteMonomerCache
	 * @throws java.io.IOException
	 * @throws org.helm.notation.MonomerException
	 */
	public synchronized void setMonomerCache(MonomerCache remoteMonomerCache)
			throws IOException, MonomerException {
		monomerDB = remoteMonomerCache.getMonomerDB();
		attachmentDB = remoteMonomerCache.getAttachmentDB();
		smilesMonomerDB = remoteMonomerCache.getSmilesMonomerDB();

		dbChanged = true;
	}

	/**
	 * 
	 * @param remoteMonomerCache
	 * @return localMonomer and remoteMonomer mismatch, key is local, value is
	 *         remote
	 * @throws java.io.IOException
	 * @throws org.helm.notation.MonomerException
	 */
	public synchronized Map<Monomer, Monomer> getConflictedMonomerMap(
			MonomerCache remoteMonomerCache) throws IOException,
			MonomerException {
		Map<String, Map<String, Monomer>> remoteMonomerDB = remoteMonomerCache
				.getMonomerDB();
		Map<String, Monomer> remoteSmilesDB = remoteMonomerCache
				.getSmilesMonomerDB();

		Map<Monomer, Monomer> map = new HashMap<Monomer, Monomer>();
		List<Monomer> newMonomers = getNewMonomers(monomerDB);
		if (newMonomers.size() > 0) {

			for (int i = 0; i < newMonomers.size(); i++) {
				Monomer local = newMonomers.get(i);
				if (remoteMonomerDB.containsKey(local.getPolymerType())) {
					Map<String, Monomer> monomers = remoteMonomerDB.get(local
							.getPolymerType());

					if (monomers.containsKey(local.getAlternateId())) {
						Monomer remote = monomers.get(local.getAlternateId());
						if (local.getCanSMILES().equals(remote.getCanSMILES())) {
							logger.log(Level.INFO, "Perfect Match");
						} else {
							map.put(local, remote);
						}
					} else {
						if (remoteSmilesDB.containsKey(local.getCanSMILES())) {
							Monomer remote = remoteSmilesDB.get(local
									.getCanSMILES());
							map.put(local, remote);
						} else {
							logger.log(Level.INFO, "Really New");
						}
					}
				} else {
					logger.log(Level.INFO, "New Polymer Type");
				}
			}
		}
		return map;
	}

	private List<Monomer> getNewMonomers(
			Map<String, Map<String, Monomer>> monomerDB) {
		List<Monomer> l = new ArrayList<Monomer>();
		Set typeSet = monomerDB.keySet();
		for (Iterator i = typeSet.iterator(); i.hasNext();) {
			String polymerType = (String) i.next();
			Map<String, Monomer> map = (Map<String, Monomer>) monomerDB
					.get(polymerType);
			Object[] monomers = map.values().toArray();
			for (int j = 0; j < monomers.length; j++) {
				Monomer m = (Monomer) monomers[j];
				if (m.isNewMonomer()) {
					l.add(m);
				}
			}
		}
		return l;
	}

	private static MonomerCache buildMonomerCacheFromXML(
			InputStream monomerDBInputStream) throws MonomerException,
			IOException, JDOMException {
		if (null == builder) {
			setupBuilder();
		}
		Document doc = builder.build(monomerDBInputStream);
		Element root = doc.getRootElement();

		Element polymerList = root.getChild(POLYMER_LIST_ELEMENT,
				root.getNamespace());
		Element attachmentList = root.getChild(ATTACHMENT_LIST_ELEMENT,
				root.getNamespace());

		Map<String, Map<String, Monomer>> newMonomerDB = buildMonomerDB(polymerList);
		Map<String, Attachment> newAttachmentDB = buildAttachmentDB(attachmentList);
		Map<String, Monomer> newSmilesMonomerDB = buildSmilesMonomerDB(newMonomerDB);

		MonomerCache cache = new MonomerCache();
		cache.setMonomerDB(newMonomerDB);
		cache.setAttachmentDB(newAttachmentDB);
		cache.setSmilesMonomerDB(newSmilesMonomerDB);

		return cache;
	}

	private static String buildMonomerDbXMLFromCache(MonomerCache cache)
			throws MonomerException {
		XMLOutputter outputer = new XMLOutputter(Format.getCompactFormat());

		StringBuilder sb = new StringBuilder();
		sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?><MONOMER_DB xmlns=\"lmr\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">");

		Map<String, Map<String, Monomer>> mDB = cache.getMonomerDB();
		Element polymerListElement = new Element(POLYMER_LIST_ELEMENT);
		Set<String> polymerTypeSet = mDB.keySet();
		for (Iterator i = polymerTypeSet.iterator(); i.hasNext();) {
			String polymerType = (String) i.next();
			Element polymerElement = new Element(POLYMER_ELEMENT);
			Attribute att = new Attribute(POLYMER_TYPE_ATTRIBUTE, polymerType);
			polymerElement.setAttribute(att);
			polymerListElement.getChildren().add(polymerElement);

			Map<String, Monomer> monomerMap = mDB.get(polymerType);
			Set<String> monomerSet = monomerMap.keySet();

			for (Iterator it = monomerSet.iterator(); it.hasNext();) {
				String monomerID = (String) it.next();
				Monomer m = monomerMap.get(monomerID);
				Element monomerElement = MonomerParser.getMonomerElement(m);
				polymerElement.getChildren().add(monomerElement);
			}
		}
		String polymerListString = outputer.outputString(polymerListElement);
		sb.append(polymerListString);

		Map<String, Attachment> aDB = cache.getAttachmentDB();
		Element attachmentListElement = new Element(ATTACHMENT_LIST_ELEMENT);
		Set<String> attachmentSet = aDB.keySet();
		for (Iterator itr = attachmentSet.iterator(); itr.hasNext();) {
			String attachmentID = (String) itr.next();
			Attachment attachment = aDB.get(attachmentID);
			Element attachmentElement = MonomerParser
					.getAttachementElement(attachment);
			attachmentListElement.getChildren().add(attachmentElement);
		}
		String attachmentListString = outputer
				.outputString(attachmentListElement);
		sb.append(attachmentListString);

		sb.append("</MONOMER_DB>");

		return sb.toString();
	}

	/**
	 * This method is called during startup, use serialized version if exists,
	 * otherwise use XML version (First from local, then from jar)
	 * 
	 * @throws org.helm.notation.MonomerException
	 * @throws java.io.IOException
	 * @throws org.jdom.JDOMException
	 */
	private static void initializeMonomerCache() throws MonomerException,
			IOException, JDOMException {
		MonomerCache cache = null;
		InputStream in = null;

		File cacheFile = new File(MONOMER_CACHE_FILE_PATH);
		if (cacheFile.exists()) {
			try {
				cache = deserializeMonomerCache(MONOMER_CACHE_FILE_PATH);
				validate(cache.getMonomerDB());
				logger.log(Level.INFO, MONOMER_CACHE_FILE_PATH
						+ " is used for monomer cache initialization");
			} catch (Exception e) {
				logger.log(Level.INFO,
						"Unable to use local monomer cache file: "
								+ MONOMER_CACHE_FILE_NAME);
				cacheFile.delete();
				logger.log(Level.INFO, "Deleted local monomer cache file: "
						+ MONOMER_CACHE_FILE_NAME);
			}
		}

		File localMonomerDBFile = new File(MONOMER_DB_FILE_PATH);
		if (null == cache && localMonomerDBFile.exists()) {
			try {
				in = new FileInputStream(MONOMER_DB_FILE_PATH);
				cache = buildMonomerCacheFromXML(in);
				validate(cache.getMonomerDB());
				logger.log(Level.INFO, MONOMER_DB_FILE_PATH
						+ " is used for monomer cache initialization");
			} catch (Exception e) {
				logger.log(Level.INFO, "Unable to use local monomer DB file: "
						+ MONOMER_DB_FILE_NAME);
				localMonomerDBFile.delete();
				logger.log(Level.INFO, "Deleted local monomer DB file: "
						+ MONOMER_DB_FILE_NAME);
			}
		}

		if (null == cache) {
			in = MonomerFactory.class
					.getResourceAsStream(MONOMER_DB_XML_RESOURCE);
			cache = buildMonomerCacheFromXML(in);
			validate(cache.getMonomerDB());
			logger.log(Level.INFO, MONOMER_DB_XML_RESOURCE
					+ " is used for monomer cache initialization");
		}

		monomerDB = cache.getMonomerDB();
		attachmentDB = cache.getAttachmentDB();
		smilesMonomerDB = cache.getSmilesMonomerDB();

		dbChanged = true;
	}

	/**
	 * save monomerCache to disk file
	 * 
	 * @throws java.io.IOException
	 */
	public void saveMonomerCache() throws IOException, MonomerException {
		File f = new File(NOTATION_DIRECTORY);
		if (!f.exists()) {
			f.mkdir();
		}
		MonomerCache cache = new MonomerCache();
		cache.setMonomerDB(getMonomerDB());
		cache.setAttachmentDB(getAttachmentDB());
		cache.setSmilesMonomerDB(getSmilesMonomerDB());
		serializeMonomerCache(cache, MONOMER_CACHE_FILE_PATH);

		String monomerDbXML = buildMonomerDbXMLFromCache(cache);
		FileOutputStream fos = new FileOutputStream(MONOMER_DB_FILE_PATH);
		fos.write(monomerDbXML.getBytes());
	}

	private static Map<String, Map<String, Monomer>> buildMonomerDB(
			Element polymerList) throws MonomerException, IOException,
			JDOMException {
		Map<String, Map<String, Monomer>> map = new HashMap<String, Map<String, Monomer>>();
		List poplymers = polymerList.getChildren();

		Iterator i = poplymers.iterator();
		while (i.hasNext()) {
			Element polymer = (Element) i.next();
			Attribute polymerType = polymer
					.getAttribute(POLYMER_TYPE_ATTRIBUTE);

			Map idMonomerMap = new HashMap<String, Monomer>();

			List monomers = polymer.getChildren();
			Iterator it = monomers.iterator();
			while (it.hasNext()) {
				Element monomer = (Element) it.next();
				Monomer m = MonomerParser.getMonomer(monomer);
				if (MonomerParser.validateMonomer(m)) {
					idMonomerMap.put(m.getAlternateId(), m);
				}
			}
			map.put(polymerType.getValue(), idMonomerMap);
		}

		return map;
	}

	private static Map<String, Attachment> buildAttachmentDB(
			Element attachmentList) throws MonomerException, IOException,
			JDOMException {
		Map<String, Attachment> map = new HashMap<String, Attachment>();

		List attachments = attachmentList.getChildren();
		Iterator i = attachments.iterator();
		while (i.hasNext()) {
			Element attachment = (Element) i.next();
			Attachment att = MonomerParser.getAttachment(attachment);

			if (MonomerParser.validateAttachement(att)) {
				map.put(att.getAlternateId(), att);
			}

		}

		return map;
	}

	private static Map<String, Monomer> buildSmilesMonomerDB(
			Map<String, Map<String, Monomer>> monomerDB) {
		Map<String, Monomer> map = new HashMap<String, Monomer>();
		Set<String> polymerSet = monomerDB.keySet();
		for (Iterator i = polymerSet.iterator(); i.hasNext();) {
			String polymer = (String) i.next();
			Map<String, Monomer> monomerMap = monomerDB.get(polymer);
			Set<String> monomerSet = monomerMap.keySet();
			for (Iterator it = monomerSet.iterator(); it.hasNext();) {
				String monomerID = (String) it.next();
				Monomer monomer = monomerMap.get(monomerID);
				String smiles = monomer.getCanSMILES();
				map.put(smiles, monomer);
			}
		}
		return map;
	}

	private static boolean validate(Map<String, Map<String, Monomer>> monomerDB)
			throws MonomerException, IOException {
		Set<String> polymers = monomerDB.keySet();
		for (String polymer : polymers) {
			Map<String, Monomer> monomerMap = monomerDB.get(polymer);
			Set<String> monomers = monomerMap.keySet();
			for (String monomer : monomers) {
				Monomer m = monomerMap.get(monomer);
				MonomerParser.validateMonomer(m);
			}
		}
		return true;
	}

	public static void finalizeMonomerCache() {

		monomerDB = null;
		attachmentDB = null;
		smilesMonomerDB = null;
		dbChanged = true;
		// externalMonomerDB=null;
		instance = null;
	}

}
