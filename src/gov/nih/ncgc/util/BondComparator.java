// $Id: BondComparator.java 2440 2008-12-01 22:18:46Z nguyenda $
package gov.nih.ncgc.util;

import chemaxon.struc.MolBond;

public interface BondComparator {
    public boolean match (MolBond a, MolBond b);
}
