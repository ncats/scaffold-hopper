// $Id: DefaultAtomComparator.java 3566 2009-11-15 06:21:22Z nguyenda $

package gov.nih.ncgc.util;

import java.util.BitSet;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;

public class DefaultAtomComparator implements AtomComparator {
    boolean literal = false; // treat query atoms as literal?

    public DefaultAtomComparator () {
        this (false);
    }

    public DefaultAtomComparator (boolean literal) {
        this.literal = literal;
    }

    public boolean isLiteral () { return literal; }
    public void setLiteral (boolean literal) { this.literal = literal; }

    public boolean match (MolAtom a, MolAtom b) {
        int atnoA = a.getAtno(), atnoB = b.getAtno();

        if (!literal) {
            if (atnoA == MolAtom.ANY || atnoB == MolAtom.ANY) {
                return true;
            }

            // treat pseudo atom same as ANY
            if (atnoA == MolAtom.PSEUDO || atnoB == MolAtom.PSEUDO) {
                return true;
            }
        
            if (atnoA == MolAtom.RGROUP && atnoB == MolAtom.RGROUP) {
                return a.getRgroup() == b.getRgroup();
            }
            else if (atnoA == MolAtom.RGROUP || atnoB == MolAtom.RGROUP) {
                return true;
            }
        }

        BitSet aset = new BitSet ();
        BitSet bset = new BitSet ();

        if (!literal && atnoA == MolAtom.LIST) {
            int[] list = a.getList();
            for (int i = 0; i < list.length; ++i) {
                aset.set(list[i]);
            }
        }
        else {
            aset.set(atnoA);
        }

        if (!literal && atnoB == MolAtom.LIST) {
            int[] list = b.getList();
            for (int i = 0; i < list.length; ++i) {
                bset.set(list[i]);
            }
        }
        else {
            bset.set(atnoB);
        }

        /*      
                return (aset.intersects(bset) 
                && (a.isTerminalAtom() || b.isTerminalAtom() ||
                a.getHybridizationState() 
                == b.getHybridizationState()));
        */
        return (aset.intersects(bset)  
                && a.hasAromaticBond() == b.hasAromaticBond()
                );
    }
}
