// $Id$

package gov.nih.ncgc.util;

import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;


public class SaltIdentifier {
    private static Logger logger = Logger.getLogger
        (SaltIdentifier.class.getName());

    private static final String SALT_DATA = "/resources/SaltData.smi";
    private static volatile SaltIdentifier INSTANCE;

    private Map<String, String> salts = 
        new ConcurrentHashMap<String, String>();

    private SaltIdentifier () throws Exception {
        //System.err.println(SaltIdentifier.class.getResource(SALT_DATA));
        MolImporter molimp = new MolImporter 
            (SaltIdentifier.class.getResourceAsStream(SALT_DATA));

        SMIRKS[] smirks = new SMIRKS[MolStandardizer.NormalizedRules.length];
        for (int i = 0; i < smirks.length; ++i) {
            String rule = MolStandardizer.NormalizedRules[i];
            smirks[i] = new SMIRKS (rule);
        }

        SMIRKS[] neutralize = 
            new SMIRKS[MolStandardizer.NeutralizedRules.length];
        for (int i = 0; i < neutralize.length; ++i) {
            String rule = MolStandardizer.NeutralizedRules[i];
            neutralize[i] = new SMIRKS (rule);
        }

        TautomerGenerator tg = new TautomerGenerator ();
        for (Molecule mol = new Molecule (); molimp.read(mol); ) {
            String name = mol.getName();
            //neutralize (mol);

            for (SMIRKS tx : smirks) {
                tx.transform(mol);
            }

            tg.setMolecule(mol);
            tg.run();

            mol = tg.getCanonicalTautomer();
            for (SMIRKS tx : neutralize) {
                tx.transform(mol);
            }
            
            String[] keys = MolStandardizer.hashKeyArray(mol);
            // ignore stereo....
            salts.put(keys[0]+keys[1], name);
        }
        molimp.close();
        logger.info("Salt data initialized!");
    }

    void neutralize (Molecule m) {
        for (MolAtom a : m.getAtomArray()) {
            if (a.getAtno() == 8) {
                a.setCharge(0);
            }
        }
    }

    protected static Molecule clean (Molecule m) {
        Molecule c = m.cloneMolecule();
        for (MolAtom a : c.getAtomArray()) {
            a.setAtomMap(0);
        }
        return c;
    }

    // assume that molecule has been properly standardized
    public String identifySaltOrSolvent (Molecule m) {
        String[] keys = MolStandardizer.hashKeyArray(m);
        return salts.get(keys[0]+keys[1]);
    }

    public boolean isSaltOrSolvent (Molecule m) {
        String[] keys = MolStandardizer.hashKeyArray(m);
        return salts.containsKey(keys[0]+keys[1]);
    }

    public static SaltIdentifier getInstance () {
        if (INSTANCE == null) {
            synchronized (SaltIdentifier.class) {
                if (INSTANCE == null) {
                    try {
                        INSTANCE = new SaltIdentifier ();
                    }
                    catch (Exception ex) {
                        logger.log(Level.SEVERE, "Can't initialize class "
                                   +SaltIdentifier.class.getName(), ex);
                    }
                }
            }
        }
        return INSTANCE;
    }

    public static void main (String[] argv) throws Exception {
        if (argv.length == 0) {
            System.err.println("SaltIdentifier FILES...");
            System.exit(1);
        }

        SaltIdentifier id = SaltIdentifier.getInstance();
        for (String file : argv) {
            MolImporter molimp = new MolImporter (file);
            for (Molecule m = new Molecule (); molimp.read(m); ) {
                System.out.println(m.getName() + "\t"+id.isSaltOrSolvent(m));
            }
            molimp.close();
        }
    }
}
