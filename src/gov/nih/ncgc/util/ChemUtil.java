// $Id: ChemUtil.java 3906 2010-01-12 19:44:54Z nguyenda $

package gov.nih.ncgc.util;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.*;
import chemaxon.formats.*;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.*;


public class ChemUtil {
    static private Logger logger = Logger.getLogger(ChemUtil.class.getName());
    static boolean DEBUG = false;
    static {
        try {
            DEBUG = Boolean.getBoolean("chemutil.debug");
        }
        catch (Exception ex) {
        }
    }

    static public final int BITCOUNT[] = {
        0,1,1,2,1,2,2,3,
        1,2,2,3,2,3,3,4,
        1,2,2,3,2,3,3,4,
        2,3,3,4,3,4,4,5,
        1,2,2,3,2,3,3,4,
        2,3,3,4,3,4,4,5,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        1,2,2,3,2,3,3,4,
        2,3,3,4,3,4,4,5,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        3,4,4,5,4,5,5,6,
        4,5,5,6,5,6,6,7,
        1,2,2,3,2,3,3,4,
        2,3,3,4,3,4,4,5,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        3,4,4,5,4,5,5,6,
        4,5,5,6,5,6,6,7,
        2,3,3,4,3,4,4,5,
        3,4,4,5,4,5,5,6,
        3,4,4,5,4,5,5,6,
        4,5,5,6,5,6,6,7,
        3,4,4,5,4,5,5,6,
        4,5,5,6,5,6,6,7,
        4,5,5,6,5,6,6,7,
        5,6,6,7,6,7,7,8
    };


    /**
     * Any change to calcMolScore must update this value.  The 
     * MolStandardizer class is very sensitive to how calMolScore
     * behaves!
     */
    public static final int MOLSCORE_VERSION = 0x1;


    static class IntPair implements Comparable<IntPair> {
        public int ival1, ival2;
        public IntPair (int ival1, int ival2) {
            this.ival1 = ival1;
            this.ival2 = ival2;
        }
        public int compareTo (IntPair ip) {
            return ival2 - ip.ival2;
        }
    }

    static public int countBits (int b32) {
        int x0 = (b32 >> 24) & 0xff;
        int x1 = (b32 >> 16) & 0xff;
        int x2 = (b32 >>  8) & 0xff;
        int x3 = (b32 & 0xff);
        return BITCOUNT[x0] + BITCOUNT[x1] + BITCOUNT[x2] + BITCOUNT[x3];
    }

    static public int countBits (byte b8) {
        return BITCOUNT[b8 & 0xff];
    }

    static private int compareFingerprints (byte[] target, byte[] query) {
        /*
          if (target.length != query.length) {
          throw new IllegalArgumentException 
          ("Fingerprints are not of the same sizes!");
          }
        */

        int Nt = 0, Nq = 0;
        for (int i = 0; i < target.length; ++i) {
            Nt += BITCOUNT[(int)(target[i] & ~query[i] & 0xff)];
            Nq += BITCOUNT[(int)(~target[i] & query[i] & 0xff)];
        }
        return Nt - Nq;
    }

    // check whether fp1 contains fp2
    protected static boolean contains (byte[] fp1, byte[] fp2) {
        assert (fp1.length == fp2.length);
        for (int i = 0; i < fp1.length; ++i) {
            if ((byte)(fp1[i] & fp2[i]) != fp2[i])
                return false;
        }
        return true;
    }

    protected static boolean contains (int[] fp1, int[] fp2) {
        assert (fp1.length == fp2.length);
        for (int i = 0; i < fp1.length; ++i) {
            if ((fp1[i] & fp2[i]) != fp2[i])
                return false;
        }
        return true;
    }

    // make sure target & query have been properly processed (e.g., aromatize)
    public boolean contains (Molecule query) {

        // assuming query doesn't contain salts... do the obvious
        if (query.getAtomCount() > msearch.getTarget().getAtomCount())
            return false;

        boolean contains = false;
        int queryFp[] = new MolHandler(query)
            .generateFingerprintInInts(32, 1, 8);

        //  int res = compareFingerprints (targetFp, queryFp);
        //System.out.print(" res=" + res + " ");
        //if (res >= 0) {
        if (contains (targetFp, queryFp)) {
            msearch.setQuery(query);
            try { contains = msearch.isMatching(); }
            catch (SearchException ex) { ex.printStackTrace(); }
        }
        return contains;
    }

    private int[] targetFp; 
    private MolSearch msearch;
    public ChemUtil (Molecule mol) {
        // 1024 bit fingerprint
        targetFp = new MolHandler(mol).generateFingerprintInInts(32, 1, 8);

        msearch = new MolSearch ();
        msearch.setTarget(mol);
    }

    public static int complexity (String smarts) {
        Molecule mol = toMolecule (smarts);
        if (mol != null) return complexity (mol);
        return -1;
    }

    public static int complexity (Molecule mol) {
        int index = 0;
    
        if (mol != null) {
            int[][] sssr = mol.getSSSR();
            for (int i = 0; i < sssr.length; ++i) {
                index += sssr[i].length * 6;
            }
        
            for (int i = 0; i < mol.getAtomCount(); ++i) {
                MolAtom a = mol.getAtom(i);
                int nb = a.getBondCount();
                switch (nb) {
                case 4: index += 24; break;
                case 3: index += 12; break;
                case 2: index += 6; break;
                case 1: index += 3; break;
                }
        
                index += a.getAtno() == 6 ? 3 : 6;
            }
        }
        return index;
    }

    public static int[] topologyInvariant (Molecule mol) {
        Molecule m = mol.cloneMolecule();
        m.hydrogenize(false);
        m.expandSgroups();

        for (MolAtom a : m.getAtomArray()) {
            a.setAtno(6);
            a.setRadical(0);
            a.setCharge(0);
            a.setFlags(0);
        }
        for (MolBond b : m.getBondArray()) {
            b.setFlags(0);
            b.setType(1);
        }

        int[] map = new int[m.getAtomCount()];
        m.getGrinv(map);

        return map;
    }

    public static int[] graphInvariantOrder (Molecule mol) {
        int[] gi = new int[mol.getAtomCount()];
        mol.getGrinv(gi);
        IntPair[] pairs = new IntPair[gi.length];
        for (int i = 0; i < pairs.length; ++i) {
            pairs[i] = new IntPair (i, gi[i]);
        }
        Arrays.sort(pairs);
        for (int i = 0; i < pairs.length; ++i) {
            gi[i] = pairs[i].ival1;
        }
        pairs = null;
        return gi;
    }

    public static int[] graphInvariantOrder (Molecule mol, final int[] rank) {
        if (rank.length != mol.getAtomCount()) {
            throw new IllegalArgumentException
                ("Input rank doesn't match number of atoms in molecule");
        }

        int[] gi = new int[mol.getAtomCount()];
        mol.getGrinv(gi, Molecule.GRINV_NOHYDROGEN);

        IntPair[] pairs = new IntPair[gi.length];
        for (int i = 0; i < pairs.length; ++i) {
            pairs[i] = new IntPair (i, gi[i]);
        }

        Arrays.sort(pairs, new Comparator<IntPair>() {
                public int compare (IntPair ip1, IntPair ip2) {
                    int d = ip1.ival2 - ip2.ival2;
                    if (d == 0) {
                        d = rank[ip1.ival1] - rank[ip2.ival1];
                    }
                    return d;
                }
            });

        for (int i = 0; i < pairs.length; ++i) {
            //System.out.println(i + ": " + pairs[i].ival1 + " " + pairs[i].ival2 + " " + rank[pairs[i].ival1]);
            gi[i] = pairs[i].ival1;
        }
        pairs = null;
        return gi;
    }


    public static int[] graphInvariantOrder (Molecule mol, Molecule ref) {
        if (mol.getAtomCount() != ref.getAtomCount()) {
            throw new IllegalArgumentException 
                ("Input and ref molecules don't match");
        }
        int[] rank = new int[ref.getAtomCount()];
        for (int i = 0; i < rank.length; ++i) {
            rank[i] = ref.getAtom(i).getAtno();
        }
        return graphInvariantOrder (mol, rank);
    }


    /*
     * Unlike SSSR, this function considers fused and spiro rings
     *  as part of the same ring system.
     */
    public static int[][] getRingSystems (Molecule m) {
        return getRingSystems (m.getSSSR());
    }

    public static int[][] getRingSystems (int[][] sssr) {

        BitSet[] rings = new BitSet[sssr.length];
        for (int i = 0; i < sssr.length; ++i) {
            BitSet ri = new BitSet (sssr[i].length);
            for (int k = 0; k < sssr[i].length; ++k) {
                ri.set(sssr[i][k]);
            }
            rings[i] = ri;
        }

        UnionFind eqv = new UnionFind (sssr.length);
        for (int i = 0; i < rings.length; ++i) {
            for (int j = i+1; j < rings.length; ++j) {
                // see if rings ri & rj share common atoms
                if (rings[i].intersects(rings[j])) {
                    eqv.union(i, j);
                }
            }
        }
        rings = null;

        // equivalence classes
        int[][] eqc = eqv.getComponents();
        int[][] ringSystems = new int[eqc.length][];
        BitSet bs = new BitSet ();
        for (int i = 0; i < eqc.length; ++i) {
            bs.clear();
            for (int j = 0; j < eqc[i].length; ++j) {
                int[] r = sssr[eqc[i][j]];
                for (int n = 0; n < r.length; ++n)
                    bs.set(r[n]);
            }
            int[] ring = new int[bs.cardinality()];
            for (int j = bs.nextSetBit(0), k = 0; 
                 j >= 0; j = bs.nextSetBit(j+1)) 
                ring[k++] = j;
            ringSystems[i] = ring;
        }
        sssr = null;

        Arrays.sort(ringSystems, new Comparator<int[]>() {
                public int compare (int[] r1, int[] r2) {
                    int d = r2.length - r1.length;
                    for (int i = 0; d==0 && i < r1.length; ++i) {
                        d = r1[i] - r2[i];
                    }
                    return d;
                }
            });
    
        return ringSystems;
    }


    // the scoring function below is a modified version of the one 
    // described by
    //   Oellien et al. j. chem inf model 2006, 46, 2342-54
    public static int calcMolScore (Molecule m) {
        int score = 0;

        Molecule mol = m.cloneMolecule();
        //mol.aromatize(Molecule.AROM_BASIC);
        mol.aromatize();

        int[][] sssr = mol.getSSSR();
        int[] rsizes = mol.getSmallestRingSizeForIdx();

        /* don't use this method... very expensive!
           int[][][] rings = tau.getAromaticAndAliphaticRings
           (Molecule.AROM_BASIC, true, false, 0, 0);
           score += rings[0].length *100;
        */

        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond b = bonds[i];
            int a1 = b.getAtom1().getAtno();
            int a2 = b.getAtom2().getAtno();

            if (b.getType() == 2) {
                // each double bond between carbon and heteroatoms +1
                if ((a1 == 6 && a2 != 6) || (a1 != 6 && a2 == 6)) {
                    score += 1;
                }
                // discourage carbon = carbon
                else if (a1 == 6 && a2 == 6 && !mol.isRingBond(i)) {
                    score -= 55;
                }
                // each double bond between oxygen and nitrogen +2
                else if ((a1 == 7 && a2 == 8) || (a1 == 8 && a2 == 7)) {
                    score += 2;
                }

                if (b.getAtom1().isTerminalAtom()
                    || b.getAtom2().isTerminalAtom()) {
                    score += 1;
                }

                if (!mol.isRingBond(i)) {
                    int s = b.calcStereo2();
                    switch (s) {
                    case MolBond.CIS:
                    case MolBond.TRANS:
                    case MolBond.CTUNSPEC:
                        // penalize extra e/z bonds
                        score -= 55;
                        break;
                    }
                }

                // discourage (non-terminal) =N coming off rings
                if (((a1 == 7 && rsizes[mol.indexOf(b.getAtom2())] > 0) 
                     || (a2 == 7 && rsizes[mol.indexOf(b.getAtom1())] > 0))
                    /*&& ((a1 == 7 && !b.getAtom1().isTerminalAtom())
                      || (a2 == 7 && !b.getAtom2().isTerminalAtom()))*/) {
                    score -= 4;
                }

                // each C=N[OH] fragment +4
                if ((a1 == 6 && a2 == 7) || (a1 == 7 && a2 == 6)) {
                    MolAtom xa = a1 == 6 && a2 == 7 
                        ? b.getAtom2() : b.getAtom1();
                    int nb = xa.getBondCount();
                    if (nb == 2) {
                        for (int j = 0; j < nb; ++j) {
                            MolBond xb = xa.getBond(j);
                            if (xb.getType() == 1 
                                && xb.getOtherAtom(xa).getAtno() == 8) {
                                score += 4;
                            }
                        }
                    }
                }

                // each double bond between carbon and oxygen +2
                if ((a1 == 6 && a2 == 8) || (a1 == 8 && a2 == 6)) {
                    score += 4;
                    if (rsizes[mol.indexOf(b.getAtom1())] > 0
                        || rsizes[mol.indexOf(b.getAtom2())] > 0) {
                        // favor =O coming off of a ring
                        score += 48; 
                    }
                }

                if (a1 == 8 || a2 == 8) {
                    score += 2;
                }

                if (mol.isRingBond(i)) {
                    score += 1;
                }
            }
            else if (b.getType() == 1) {
                // each P-H, S-H, Se-H, and Te-H bond -1
                if (b.getAtom1().isTerminalAtom()
                    && (a1 == 15 || a1 == 16 || a1 == 34 || a1 == 52)) {
                    score -= 1;
                }
                else if (b.getAtom2().isTerminalAtom()
                         && (a2 == 15 || a2 == 16 || a2 == 34 || a2 == 52)) {
                    score -= 1;
                }
            }
        }
        bonds = null;



        if (DEBUG) {
            logger.info("score 1: "+score);
        }

        //mol.hydrogenize(false);

        // each methyl group (applying a penalty to structures with 
        //   terminal double bonds) +1
        MolAtom[] atoms = mol.getAtomArray(); 
        int chargeCount = 0;
        for (MolAtom a : atoms) {
            // there shouldn't be any explicit H, but what the hell...
            if ((a.getImplicitHcount() + a.getExplicitHcount()) == 3) {
                score += 1;
            }

            // not in the paper
            int c = 2*a.getCharge(); // scale the charge values
            if (c != 0) {
                switch (a.getAtno()) {
                case 7: score += c; break;
                case 8: score -= c; break;
                default: // discourage charged ions
                    score -= Math.abs(c); 
                    break;
                }
                ++chargeCount;
            }
        }
        // penalize molecule with too many charged atoms
        score -= chargeCount*10;

        if (DEBUG) {
            logger.info("score 2: "+score);
        }


        for (int i = 0; i < sssr.length; ++i) {
            int naro = 0;
            for (int j = 0; j < sssr[i].length; ++j) {
                if (atoms[sssr[i][j]].hasAromaticBond()) {
                    ++naro;
                }
            }

            if (naro == sssr[i].length) {
                //each aromatic ring system +100
                score += 100;
            }
        }
        atoms = null;

        if (DEBUG) {
            logger.info("score 3: "+score);
        }

        /*
         * we have to do this hack to get the alternating single-double
         * bond to be in a consistent layout for an aromatic ring
         */
        try {
            MolHandler mh = new MolHandler (mol.toFormat("smiles:q"));
            m = mh.getMolecule();
            m.dearomatize();
            for (MolAtom a : m.getAtomArray()) {
                if (a.hasAromaticBond()) {
                    throw new Exception ("bailing out");
                }
            }
            mol = m;
            sssr = mol.getSSSR();
            rsizes = mol.getSmallestRingSizeForIdx();
        }
        catch (Exception ex) {
            mol.dearomatize();
        }
    
        // prefer double-bond shared between rings
        Map<Integer, Integer> ringCount = new HashMap<Integer, Integer>();
    
        for (int i = 0; i < sssr.length; ++i) {
            for (int j = 0; j < sssr[i].length; ++j) {
                Integer c = ringCount.get(sssr[i][j]);
                ringCount.put(sssr[i][j], c == null ? 1 : (c+1));
            }
        }
    
        for (MolBond b : mol.getBondArray()) {
            int a1 = mol.indexOf(b.getAtom1());
            int a2 = mol.indexOf(b.getAtom2());
        
            if (b.getType() == 2 && rsizes[a1] > 0 && rsizes[a2] > 0) {
                if (ringCount.get(a1) > 1 && ringCount.get(a2) > 1) {
                    score += 1;
                }
            }
        }
        ringCount.clear();
        sssr = null;

        if (DEBUG) {
            logger.info("score 4: "+score);
        }

        return score;
    }

    static public String generateSkeleton (Molecule m) {
        List<MolAtom> remove = new ArrayList<MolAtom>();
        for (MolAtom a : m.getAtomArray()) {
            if (a.isTerminalAtom()) {
                // remove all terminal atoms
                remove.add(a);
            }
            else {
                a.setRadical(0);
                a.setCharge(0);
                a.setFlags(0);
                a.setAtno(6);
            }
        }
        for (MolBond b : m.getBondArray()) {
            b.setFlags(0);
            b.setType(1);
        }

        for (MolAtom a : remove) {
            m.removeNode(a);
        }

        return m.toFormat("cxsmarts");
    }

    public static void resetEZ (Molecule mol, int flags) {
        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond b = bonds[i];
            if (b.getType() == 2) {
                MolAtom a1 = b.getCTAtom1();
                MolAtom a4 = b.getCTAtom4();
                if (a1 != null && a4 != null) {
                    b.setStereo2Flags(a1, a4, flags);
                }
            }
        }
    }

    public static void resetEZ (Molecule mol) {
        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond b = bonds[i];
            int type = b.getType();
            if (type == 2) {
                b.setFlags(0, MolBond.CTUMASK);
            }
            else if (type == 1) {
                int flags = b.getFlags() & MolBond.STEREO1_MASK;
                if (flags  == (MolBond.UP | MolBond.DOWN)) {
                    // WAVY flag... 
                    b.setFlags(0, MolBond.STEREO1_MASK);
                }
            }
        }
    }

    public static Molecule toMolecule (String smarts) {
        try { return new MolHandler (smarts).getMolecule(); }
        catch (Exception ex) { /*ex.printStackTrace();*/ }
        return null;
    }

    public static Molecule createMolecule 
        (Molecule query, Molecule target, int[] matches) {
        return createMolecule (query, target, matches, 0);
    }

    /**
     * Create a molecule based on the matches array.  Only substructure
     * defined by matches[i] >= 0 will be created.  The mapdir specifies
     * the atom mapping: none mapdir = 0, query mapdir < 0, 
     * and target mapdir > 0
     */
    public static Molecule createMolecule 
        (Molecule query, Molecule target, int[] matches, int mapdir) {
        return createMolecule (query, target, matches, 
                               new DefaultBondComparator (), mapdir);
    }

    public static Molecule createMolecule 
        (Molecule query, Molecule target, int[] matches, 
         BondComparator bc, int mapdir) {

        Molecule mol = new Molecule ();
        MolBond[] bonds = query.getBondArray();
        MolAtom[] atoms = new MolAtom[matches.length];
        int[][] ttab = target.getBtab();
        Map<MolAtom, Integer> qmap = new HashMap<MolAtom, Integer>();
        Map<MolAtom, Integer> tmap = new HashMap<MolAtom, Integer>();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond bond = bonds[i];
            int u1 = query.indexOf(bond.getAtom1());
            int u2 = query.indexOf(bond.getAtom2());
            int v1 = matches[u1], v2 = matches[u2];
            if (v1 >= 0 && v2 >= 0) {
                // this bond is in the result structure
                MolAtom a1 = atoms[u1];
                if (a1 == null) {
                    atoms[u1] = a1 = 
                        new MolAtom (bond.getAtom1().getAtno());
                    //(MolAtom)bond.getAtom1().clone();
                    /*
                      a1.setImplicitHcount(0);
                      a1.setRadical(0);
                      a1.setCharge(0);
                    */
                    if (mapdir < 0) {
                        a1.setAtomMap(u1+1);
                    }
                    else if (mapdir > 0) {
                        a1.setAtomMap(v1+1);
                    }
                    else {
                        a1.setAtomMap(0);
                    }
                    qmap.put(a1, u1);
                    tmap.put(a1, v1);
                    mol.add(a1);
                }
                MolAtom a2 = atoms[u2];
                if (a2 == null) {
                    atoms[u2] = a2 = 
                        new MolAtom (bond.getAtom2().getAtno());
                    //(MolAtom)bond.getAtom2().clone();
                    /*
                      a2.setImplicitHcount(0);
                      a2.setRadical(0);
                      a2.setCharge(0);
                    */
                    if (mapdir < 0) {
                        a2.setAtomMap(u2+1);
                    }
                    else if (mapdir > 0) {
                        a2.setAtomMap(v2+1);
                    }
                    else {
                        a2.setAtomMap(0);
                    }
                    qmap.put(a2, u2);
                    tmap.put(a2, v2);
                    mol.add(a2);
                }

                if (ttab[v1][v2] >= 0 
                    && bc.match(bond, target.getBond(ttab[v1][v2]))) {
                    MolBond b = new MolBond (a1, a2);
                    //bond.cloneBond(a1, a2);

                    b.setFlags(bond.getType(), MolBond.TYPE_MASK);
                    // this forces the valences of a1 & a2 to be checked...
                    //b.setType(bond.getType());
            
                    mol.add(b);
                }
            }
        }
        atoms = null;
        bonds = null;

        mol.valenceCheck();
        if (mol.getAtomCount() > 0) {
            try {
                mol.clean(2, null);
                mol.aromatize();

                /*
                  if (!mol.isQuery()) {
                  String smiles = mol.toFormat("smarts:q0-H");
                  MolImporter.importMol(smiles, mol);
                  }
                */
                resetEZ (mol);

                atoms = mol.getAtomArray();
                int[] qm = new int[atoms.length];
                int[] tm = new int[atoms.length];
                for (int i = 0; i < atoms.length; ++i) {
                    qm[i] = qmap.get(atoms[i]);
                    tm[i] = tmap.get(atoms[i]);
                }
                mol.setPropertyObject("QMAP", qm);
                mol.setPropertyObject("TMAP", tm);

                /*
                  for (MolAtom atom : mol.getAtomArray()) {
                  int atno = atom.getAtno();
                  if (atno == 7) {
                  int nb = atom.getBondCount();
                  int valence = 0;
                  for (int i = 0; i < nb; ++i) {
                  valence += atom.getBond(i).getType();
                  }
                  if (nb == 3 && (valence == 4 
                  || atom.hasAromaticBond())) {
                  atom.setCharge(1);
                  }
                  }
                  }
                */
            }
            catch (Exception ex) {
                //ex.printStackTrace();
                logger.log(Level.SEVERE, "Can't read molecule", ex);
            }
        }

        return mol;
    }

    public static int[] neighborDensity (int[] rho, MoleculeGraph g) {
        int acount = g.getAtomCount();
        if (rho == null) {
            rho = new int[acount];
        }
        else if (rho.length < acount) {
            throw new IllegalArgumentException
                ("Bad array dimension: "+rho.length);
        }

        for (int i = 0; i < acount; ++i) {
            int nb = g.getNeighborCount(i);

            IntPair[] pairs = new IntPair[nb];
            for (int j = 0; j < nb; ++j) {
                int kb = g.getNeighbor(i, j);
                pairs[j] = new IntPair (j, g.getNeighborCount(kb));
            }
            Arrays.sort(pairs); // sort in non-increasing order

            int r = 0;
            for (int j = 0; j < nb; ++j) {
                if (pairs[j].ival2 > r) {
                    ++r;
                }
                else {
                    break;
                }
            }
            rho[i] = r;
        }
        return rho;
    }

    public static String encode (Molecule mol) {
        try {
            String str = mol.toFormat("mrv");
            return Base64.encodeBytes(str.getBytes("utf-8"));
        }
        catch (Exception ex) {
            logger.log(Level.SEVERE, "Can't encode molecule base64", ex);
        }
        return null;
    }

    public static Molecule decode (String str) {
        try {
            byte[] bytes = Base64.decode(str);
            MolHandler mh = new MolHandler (new String (bytes));
            return mh.getMolecule();
        }
        catch (Exception ex) {
            logger.log(Level.SEVERE, "Can't decode string "+str, ex);
        }
        return null;
    }

    // determines wether a contains b
    public static boolean contains (Molecule a, Molecule b) {
        return new ChemUtil(a).contains(b);
    }

    public static void main (String argv[]) throws Exception {
        if (argv.length == 0) {
            System.out.println("ChemUtil FILE");
            System.exit(1);
        }

        MolImporter molimp = new MolImporter (argv[0]);
        Vector<Molecule> mols = new Vector<Molecule>();
        for (Molecule mol; (mol = molimp.read()) != null; ) {
            mol.aromatize();
            String name = mol.getName();
            if (name == null || name.length() == 0) {
                for (int i = 0; i < mol.getPropertyCount(); ++i) {
                    name = mol.getProperty(mol.getPropertyKey(i));
                    if (name.length() > 0) {
                        mol.setName(name);
                        break;
                    }
                }
            }
            mols.add(mol);
        
            int[][] rs = getRingSystems (mol);
            System.out.println(name + " has "+rs.length + " ring system(s)!");
            for (int i = 0; i < rs.length; ++i) {
                System.out.print(" " + (i+1) + ":");
                for (int j = 0; j < rs[i].length; ++j) {
                    System.out.print(" "+(rs[i][j]+1));
                }
                System.out.println();
            }
        }
        molimp.close();
    }
}
