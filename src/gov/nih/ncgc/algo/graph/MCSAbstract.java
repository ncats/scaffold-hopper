// $Id: MCSAbstract.java 3551 2009-11-12 23:05:13Z nguyenda $

package gov.nih.ncgc.algo.graph;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Collections;
import java.util.ArrayList;
import java.util.List;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolFormatException;
import chemaxon.util.MolHandler;

import gov.nih.ncgc.util.UnionFind;
import gov.nih.ncgc.util.ElementData;
import gov.nih.ncgc.util.ChemUtil;
import gov.nih.ncgc.util.AtomComparator;
import gov.nih.ncgc.util.BondComparator;
import gov.nih.ncgc.util.DefaultAtomComparator;
import gov.nih.ncgc.util.DefaultBondComparator;


public abstract class MCSAbstract implements Serializable {
    static final int MIN_MATCH_SIZE = 2;

    protected AtomComparator atomComparator = new DefaultAtomComparator ();
    protected BondComparator bondComparator = new DefaultBondComparator ();

    protected Molecule query, target;
    protected int[][] qCtab, tCtab;
    protected int[][] qBtab, tBtab;
    protected MolAtom[] qAtoms, tAtoms;
    protected int minsize = MIN_MATCH_SIZE;
    
    public MCSAbstract () {
    }

    // copy constructor
    public MCSAbstract (MCSAbstract mcs) {
        this.query = mcs.query;
        this.target = mcs.target;
        this.qCtab = mcs.qCtab;
        this.tCtab = mcs.tCtab;
        this.qBtab = mcs.qBtab;
        this.tBtab = mcs.tBtab;
        this.atomComparator = mcs.atomComparator;
        this.bondComparator = mcs.bondComparator;
    }

    public void setAtomComparator (AtomComparator comparator) {
        atomComparator = comparator;
    }
    public AtomComparator getAtomComparator () { return atomComparator; }
    public void setBondComparator (BondComparator comparator) {
        bondComparator = comparator;
    }
    public BondComparator getBondComparator () { return bondComparator; }

    public void setQuery (Molecule query) {
        this.query = query;
        qCtab = query.getCtab();
        qBtab = query.getBtab();
        qAtoms = query.getAtomArray();
    }
    public void setQuery (String query) {
        try {
            MolHandler mh = new MolHandler (query);
            setQuery (mh.getMolecule());
        }
        catch (Exception ex) {
            throw new IllegalArgumentException (ex);
        }
    }
    public Molecule getQuery () { return query; }
    public void setTarget (Molecule target) {
        this.target = target; 
        tCtab = target.getCtab();
        tBtab = target.getBtab();
        tAtoms = target.getAtomArray();
    }
    public void setTarget (String target) {
        try {
            MolHandler mh = new MolHandler (target);
            setTarget (mh.getMolecule());
        }
        catch (Exception ex) {
            throw new IllegalArgumentException (ex);
        }
    }
    public Molecule getTarget () { return target; }
    public void setMolecules (Molecule query, Molecule target) {
        setQuery (query);
        setTarget (target);
    }
    public void setMolecules (String query, String target) {
        setQuery (query);
        setTarget (target);
    }

    public abstract boolean search ();
    public abstract int[] getResult ();
    public abstract int getResultSize ();

    public Molecule getResultAsMolecule () {
        return getResultAsMolecule (true);
    }

    public Molecule getResultAsMolecule (boolean atomMapping) {
        int[] result = getResult ();
        if (result == null || result.length == 0) {
            return null;
        }
        return ChemUtil.createMolecule
            (query, target, result, atomMapping ? -1 : 0);
    }

    public String getResultAsSmiles (boolean unique, boolean atomMap) {
        Molecule mol = getResultAsMolecule (atomMap);
        if (mol != null) {
            try {
                for (MolAtom a : mol.getAtomArray()) {
                    if (a.isQuery()) {
                        return mol.toFormat("cxsmarts:u");
                    }
                }
                return mol.toFormat("smiles" + (unique?":a_basu0":""));
            }
            catch (Exception ex) {
                System.err.println("** query="+query.getName() + " target="
                                   + target.getName());
                ex.printStackTrace();
            }
        }
        return null;
    }

    public String getResultAsSmiles () { 
        return getResultAsSmiles (true, false);  
    }

    protected int[][] getComponents 
        (Molecule core, int[] fmap, int[] rmap, boolean largestOnly) {
        UnionFind disjoint = new UnionFind (qAtoms.length);

        MolBond[] qbonds = query.getBondArray();
        for (int i = 0; i < qbonds.length; ++i) {
            MolBond qb = qbonds[i];
            int a1 = query.indexOf(qb.getAtom1());
            int a2 = query.indexOf(qb.getAtom2());
            int i1 = fmap[a1], i2 = fmap[a2];
            if (i1 >= 0 && i2 >= 0 && tBtab[i1][i2] >= 0) {
                disjoint.union(a1, a2);
            }
        }

        List<BitSet> qcomps = new ArrayList<BitSet>();
        for (int[] comp : disjoint.getComponents()) {
            // figure out which components to keep
            if (comp.length < minsize) { // don't bother
            }
            else {
                BitSet bs = new BitSet ();
                for (int j = 0; j < comp.length; ++j) {
                    int v = comp[j];
                    if (fmap[v] >= 0) {
                        bs.set(v);
                    }
                }
                if (bs.cardinality() > 0) {
                    qcomps.add(bs);
                }
            }
        }

        if (qcomps.isEmpty()) {
            return null;
        }

        disjoint = new UnionFind (tAtoms.length);
        MolBond[] tbonds = target.getBondArray();
        for (int i = 0; i < tbonds.length; ++i) {
            MolBond tb = tbonds[i];
            int a1 = target.indexOf(tb.getAtom1());
            int a2 = target.indexOf(tb.getAtom2());
            int i1 = rmap[a1], i2 = rmap[a2];
            if (i1 >= 0 && i2 >= 0 && qBtab[i1][i2] >= 0) {
                disjoint.union(a1, a2);
            }
        }

        List<BitSet> tcomps = new ArrayList<BitSet>();
        for (int[] comp : disjoint.getComponents()) {
            if (comp.length < minsize) {
            }
            else {
                BitSet bs = new BitSet ();
                for (int j = 0; j < comp.length; ++j) {
                    int v = comp[j];
                    if (rmap[v] >= 0) {
                        bs.set(v);
                    }
                }
                if (bs.cardinality() > 0) {
                    tcomps.add(bs);
                }
            }
        }

        if (tcomps.isEmpty()) {
            return null;
        }

        List<int[]> disconnected = new ArrayList<int[]>();
        List<Molecule> fragments = new ArrayList<Molecule>();
        for (BitSet qs : qcomps) {
            for (BitSet ts : tcomps) {
                BitSet cs = new BitSet ();
                for (int b = ts.nextSetBit(0); 
                     b >= 0; b = ts.nextSetBit(b+1)) {
                    cs.set(rmap[b]);
                }
                // now cs contains the common fragment between 
                //   query and target
                cs.and(qs);
                if (cs.cardinality() >= minsize) {
                    if (core != null) {
                        fragments.add(extractFragment (cs, fmap));
                    }
            
                    int[] c = new int[cs.cardinality()];
                    for (int i = cs.nextSetBit(0), j = 0; 
                         i >= 0; i = cs.nextSetBit(i+1)) {
                        c[j++] = i;
                    }
                    disconnected.add(c);
                }
            }
        }

        Collections.sort(disconnected, new Comparator<int[]>() {
                public int compare (int[] a, int[] b) {
                    return b.length - a.length;
                }
            });

        Collections.sort(fragments, new Comparator<Molecule>() {
                public int compare (Molecule a, Molecule b) {
                    return b.getAtomCount() - a.getAtomCount();
                }
            });

        if (core != null) {
            core.clear();
            for (Molecule frag : fragments) {
                //System.out.println(" frag " + frag.toFormat("smiles:u0-H"));
                core.fuse(frag);
                if (largestOnly) {
                    break; // these are sorted in desc order
                }
            }
        }

        /*
          System.out.println("** " + disconnected.size() + " component(s)");
          for (int[] f : disconnected) {
          System.out.print(f.length + ":");
          for (int i = 0; i < f.length; ++i) {
          System.out.print(" " + (f[i]+1));
          }
          System.out.println();
          }
        */

        return disconnected.toArray(new int[0][]);
    }

    protected Molecule extractFragment 
        (Molecule frag, BitSet atomset, BitSet bondset, int[] map) {
        if (frag == null) {
            frag = new Molecule ();
        }
        else {
            frag.clear();
        }

        MolAtom[] atoms = new MolAtom[qAtoms.length];
        MolBond[] bonds = query.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            if (!bondset.get(i)) {
                continue;
            }

            MolBond bond = bonds[i];
            int i1 = query.indexOf(bond.getAtom1());
            int i2 = query.indexOf(bond.getAtom2());
            if (atomset.get(i1) && atomset.get(i2) 
                && map[i1] >= 0 && map[i2] >= 0 
                && tBtab[map[i1]][map[i2]] >= 0) {
                MolAtom a1 = atoms[i1];
                if (a1 == null) {
                    atoms[i1] = a1 = (MolAtom)bond.getAtom1().clone();
                    a1.setImplicitHcount(0);
                    a1.setRadical(0);
                    a1.setCharge(0);
                    //a1.setAtomMap(i1+1);
                    /*
                      a1.setSetSeq(i1);
                      a1.setExtraLabelSetSeq(map[i1]);
                    */
                    frag.add(a1);
                }
                MolAtom a2 = atoms[i2];
                if (a2 == null) {
                    atoms[i2] = a2 = (MolAtom)bond.getAtom2().clone();
                    a2.setImplicitHcount(0);
                    a2.setRadical(0);
                    a2.setCharge(0);
                    //a2.setAtomMap(i2+1);
                    /*
                      a2.setSetSeq(i2);
                      a2.setExtraLabelSetSeq(map[i2]);
                    */
                    frag.add(a2);
                }
                MolBond b = bond.cloneBond(a1, a2);

                b.setFlags(0);
                // this forces the valences of a1 & a2 to be checked...
                b.setType(bond.getType());

                frag.add(b);
            }
        }

        frag.aromatize();
        frag.dearomatize();
        for (MolAtom atom : frag.getAtomArray()) {
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
        frag.aromatize();

        return frag;
    }

    protected Molecule extractFragment (BitSet atoms, int[] map) {
        BitSet bondset = new BitSet (query.getBondCount());
        bondset.flip(0, query.getBondCount()-1);
        return extractFragment (null, atoms, bondset, map);
    }
    protected Molecule extractFragment 
        (BitSet atoms, BitSet bonds, int[] map) {
        return extractFragment (null, atoms, bonds, map);
    }
}
