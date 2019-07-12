// $Id: AtomComparator.java 2439 2008-12-01 22:18:31Z nguyenda $

package gov.nih.ncgc.util;

import chemaxon.struc.MolAtom;

public interface AtomComparator {
    public boolean match (MolAtom a, MolAtom b);
}
