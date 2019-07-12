// $Id: DefaultBondComparator.java 2441 2008-12-01 22:18:58Z nguyenda $

package gov.nih.ncgc.util;

import chemaxon.struc.MolBond;

public class DefaultBondComparator implements BondComparator {
    public boolean match (MolBond a, MolBond b) {
	int atype = a.getType(), btype = b.getType();
	return atype == btype
	    || atype == MolBond.ANY
	    || btype == MolBond.ANY;
    }
}
