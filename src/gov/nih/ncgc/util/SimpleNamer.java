package gov.nih.ncgc.util;

import java.util.HashMap;
import java.util.Map;

import chemaxon.marvin.util.MolExportException;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

public class SimpleNamer {
	private static Map<String,String> knownNames=new HashMap<String,String>();
	static{
		knownNames.put(getHash("*-Cl"), "chloro");
		knownNames.put(getHash("*-Br"), "bromo");
		knownNames.put(getHash("*-I"), "iodo");
		knownNames.put(getHash("*-F"), "flouro");
		knownNames.put(getHash("*-c1ccccc1"), "phenyl");
		knownNames.put(getHash("*-C"), "methyl");
		knownNames.put(getHash("*-CC"), "ethyl");
		knownNames.put(getHash("*-CCC"), "propyl");
		knownNames.put(getHash("*-CCCC"), "n-butyl");
		knownNames.put(getHash("*-C(C)C"), "isopropyl");
		knownNames.put(getHash("*-CC(C)C"), "isobutyl");
		knownNames.put(getHash("*-C(C)(C)C"), "tertbutyl");
		knownNames.put(getHash("*-OC"), "methoxy");
		knownNames.put(getHash("*-[#6]"), "methyl");
		knownNames.put(getHash("*-[N+]([O-])=O"), "nitro");
		knownNames.put(getHash("*-[N+](=O)=O"), "nitro");
		knownNames.put(getHash("*-[N]=[N+]=[N-]"), "azide");
		knownNames.put(getHash("*-C(=O)C"), "acetyl");
		//pleasing rearrangements
		knownNames.put(getHash("*-O"), "OH");
		knownNames.put(getHash("*-CO"), "COH");
		knownNames.put(getHash("*-C(=O)O"), "COOH");
		knownNames.put(getHash("*-SC"), "SCH3");
		knownNames.put(getHash("*-N"), "NH2");
		
	}
	public static String getName(Molecule m){
		//String form=
		String n=knownNames.get(MolStandardizer.hashKey(m));
		
		if (n == null){
			try{
				String n2 = m.exportToFormat("smiles:q");
				n=knownNames.get(getHash(n2));
				if(n==null){
					if(m.getAtomCount()<8){
						n=m.getFormula();
					}else{
						n=n2;
					}
				}
			}catch(Exception e){
				try {
					String n2 = m.exportToFormat("cxsmiles:q").split(" ")[0];
					n=knownNames.get(getHash(n2));
					if(n==null){
						if(m.getAtomCount()<8){
							n=m.getFormula();
						}else{
							n=n2;
						}
					}
				} catch (MolExportException e1) {
					e1.printStackTrace();
				}
			}
		}
		return n;
	}
	private static Molecule getMol(String smi){
		MolHandler mh=new MolHandler();
		try{
		mh.setMolecule(smi);
		return mh.getMolecule();
		}catch(Exception e){
			return null;
		}
	}
	private static String getHash(String smi){
		return MolStandardizer.hashKey(getMol(smi));
	}
}
