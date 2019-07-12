// $Id: RGroupNode.java 2427 2008-11-04 21:30:33Z nguyenda $
package gov.nih.ncgc.rgroup;

import java.util.*;
import java.io.*;

import javax.swing.tree.*;

import chemaxon.struc.Molecule;
import chemaxon.struc.RgMolecule;
import chemaxon.struc.MolBond;
import chemaxon.struc.MolAtom;
import chemaxon.util.MolHandler;
import chemaxon.sss.screen.HashCode;
import chemaxon.sss.search.MolSearch;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.RGroupDecomposition;

import gov.nih.ncgc.algo.graph.MCSMaxClique;
import gov.nih.ncgc.util.ChemUtil;

public class RGroupNode extends DefaultMutableTreeNode 
    implements Comparable<RGroupNode> {

    static final int FP_SIZE = 32; // numInts
    static final int FP_DEPTH = 2;
    static final int FP_EDGES = 7;

    static final int MAX_DEPTH = 64;

    private int[] fp; // fingerprint of fragment
    private HashCode hash = new HashCode ();
    private MCSMaxClique mcs = new MCSMaxClique ();
    private MolSearch msearch = new MolSearch ();
    private MolHandler mh = new MolHandler ();

    private RgMolecule rgMol = new RgMolecule ();
    private int[] rmap = null;
    private Molecule[][] ligands = null; // ligand table

    private static int debug = Integer.getInteger("rgroup.debug", 0);

    private static Object scaffoldLock = new Object ();
    private static int scaffoldCounter = 0;
    
    public RGroupNode () {
    }
    
    public RGroupNode (Molecule fragment) {
	setFragment (fragment);
    }

    public int[] getFp () { return fp; }
    public Molecule getFragment () { 
	return (Molecule)getUserObject (); 
    }
    public String getFragmentAsSmiles () {
	Molecule frag = getFragment ();
	if (frag != null) {
	    return frag.toFormat("smiles:u0-H");
	}
	return null;
    }

    public void buildRgroup () {
	Molecule core = getFragment ();
	if (isLeaf () || core == null) {
	    return;
	}

	rgMol.setRoot(core.cloneMolecule());
	RGroupDecomposition.addRGroups(rgMol);

	RGroupDecomposition rgd = new RGroupDecomposition();
	rgd.setAttachmentType(RGroupDecomposition.ATTACHMENT_MAP);
	rgd.setQuery(rgMol);

	Map<Integer, Set<String>> rgroups 
	    = new TreeMap<Integer, Set<String>>();

	rmap = rgd.getQueryRMap();
	ligands = new Molecule[getChildCount ()][];
	for (RGroupNode child = (RGroupNode)getFirstChild ();
	     child != null; child = (RGroupNode)child.getNextSibling()) {
	    rgd.setTarget(child.getFragment());

	    try {
		int[] hit = rgd.findFirst();
		if (hit != null) {
		    Molecule[] ligs = rgd.findLigands(hit);
		    Molecule[] copy = new Molecule[ligs.length];
		    for (int i = 0; i < ligs.length; ++i) {
			if (rmap[i] != 0) {
			    Set<String> group = rgroups.get(rmap[i]);
			    if (group == null) {
				rgroups.put
				    (rmap[i], group = new HashSet<String>());
			    }
			    MolAtom[] a = ligs[i].getAtomArray();
			    if (a.length == 1 && a[0].getAtno() == 1) {
				// ignore 
			    }
			    else {
				group.add(ligs[i].toFormat("smiles:u0"));
				copy[i] = ligs[i].cloneMolecule();
			    }
			}
		    }
		    // save the ligands
		    ligands[getIndex (child)] = copy;
		}
		else {
		    System.err.println
			("** warning: No R-group found for node " + child);
		}
	    }
	    catch (chemaxon.sss.search.SearchException ex) {
		ex.printStackTrace();
	    }
	}
	
	// prune empty groups
	MolAtom[] atoms = rgMol.getAtomArray();
	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom a = atoms[i];
	    if (a.getAtno() == MolAtom.RGROUP) {
		Set<String> group = rgroups.get(a.getRgroup());
		if (group == null || group.isEmpty()) {
		    rgMol.removeNode(a);
		}
	    }
	}
	// make all 
	rgMol.hydrogenize(false);
	rgroups.clear();
    }

    public Molecule getRgroup () { return rgMol; }
    public Molecule[][] getLigands () { return ligands; }
    public Molecule[] getLigands (int child) { return ligands[child]; }
    public Molecule[] getLigands (RGroupNode child) {
	return getLigands (getIndex (child));
    }

    public void setFragment (Molecule fragment) {
	if (fragment != null) {
	    mh.setMolecule(fragment);
	    fp = mh.generateFingerprintInInts(FP_SIZE, FP_DEPTH, FP_EDGES);
	    mcs.setQuery(fragment);
	    msearch.setQuery(fragment);
	}
	setUserObject (fragment);
    }

    public String getName () { 
	Molecule m = getFragment ();
	if (m != null) {
	    return m.getName();
	}
	return null;
    }

    public int hashCode () {
	Molecule fragment = getFragment ();
	if (fragment != null) {
	    return hash.getHashCode(fragment, fp);
	}
	return -1;
    }

    protected void add0 (RGroupNode node) {
	if (isLeaf () || (getLevel()+1) >= MAX_DEPTH) {
	    add (node);
	}
	else {
	    double max = 0.;
	    RGroupNode best = null;
	    for (RGroupNode child = (RGroupNode) getFirstChild();
		 child != null; 
		 child = (RGroupNode)child.getNextSibling()) {
		double sim = child.tanimoto(node);
		if (sim > max) {
		    max = sim;
		    best = child;
		}
	    }

	    if (max < 1.) {
		best.add1(node);
	    }
	    else {
		best.add(node);
	    }
	}
    }

    protected void add1 (RGroupNode node) {
	int dir = compareTo (node);
	if (debug > 0) {
	    System.out.println("** " + getName () + " vs " 
			       + node.getName() + " => " + dir);
	}

	if (dir < 0 || hashCode() == node.hashCode()) {
	    add0 (node);
	}
	else if (dir > 0) {
	    node.add0(this);
	}
	else {
	    // create a parent node with the current and 
	    //   new fragments as children
	    mcs.setTarget(node.getFragment());

	    RGroupNode parent = (RGroupNode) getParent ();

	    mcs.setSeed(parent != null ? parent.getFragment() : null);
	    mcs.search();
		
	    Molecule core = mcs.getResultAsMolecule(false);
	    if (core != null) {
		//core.aromatize(Molecule.AROM_BASIC);
		core.aromatize();
		core.calcHybridization();
		synchronized (scaffoldLock) {
		    ++scaffoldCounter;
		}
		core.setName("Scaffold" + scaffoldCounter);
	    }
	    else if (parent != null) {
		System.err.println("** fatal error: failed to find MCS "
				   +"with mincore = " 
				   + parent.getFragmentAsSmiles());
	    }

	    if (debug > 0) {
		System.out.println("creating new parent: "
				   + (core != null ? 
				      (core.getName() + " " 
				       + core.toFormat("smiles:u0-a")) 
				      : "null"));
		System.out.println
		    ("  + " + getName () + " " + getFragmentAsSmiles());
		System.out.println
		    ("  + " + node.getName() + " " 
		     + node.getFragmentAsSmiles());
	    }

	    RGroupNode newParent = new RGroupNode (core);
	    newParent.add(this);
	    newParent.add(node);
	    if (parent != null) {
		parent.add(newParent);
	    }
	}
    }

    public void add (Molecule mol) {
	Molecule frag = getFragment ();
	if (frag == null) {
	    setFragment (mol);
	}
	else {
	    RGroupNode node = new RGroupNode (mol);
	    add1 (node);
	}
    }


    protected static void pruneChildren (RGroupNode node) {
	if (node == null || node.isLeaf()) {
	    return;
	}

	RGroupNode parent = (RGroupNode) node.getParent();
	if (debug > 0) {
	    System.out.println("## pruning node " + node.getName());
	}
	
	if (parent != null) {
	    String smiles = parent.getFragmentAsSmiles();
	    if (smiles.equals(node.getFragmentAsSmiles())) {
		// add all of node's children to its parent
		if (debug > 0) {
		    System.out.println("** merging " + parent.getName() 
				       + " and " + node.getName());
		}
		
		node.removeFromParent();
		for (RGroupNode c = (RGroupNode)node.getFirstChild();
		     c != null; ) {
		    RGroupNode n = (RGroupNode)c.getNextSibling();
		    parent.add(c);
		    c = n;
		}
		node = parent;
	    }
	}

	for (RGroupNode child = (RGroupNode)node.getFirstChild();
	     child != null; ) {
	    pruneChildren (child);
	    child = (RGroupNode)child.getNextSibling();
	}
    }

    public void prune () {
	pruneChildren (this);
    }

    public String toString () { return getName (); }

    public int compareTo (RGroupNode n) {
	int c = compare (getFp (), n.getFp());
	if (c != 0) {
	    msearch.setTarget(n.getFragment());
	    try {
		if (msearch.isMatching()) {
		    c = -1;
		}
		else {
		    n.msearch.setTarget(getFragment ());
		    c = n.msearch.isMatching() ? 1 : 0;
		}
	    }
	    catch (chemaxon.sss.search.SearchException ex) {
		ex.printStackTrace();
	    }
	}
	return c;
    }

    // return 0 if fp1 & fp2 are not completely within the other
    // -1 if fp1 is a substructure of fp2
    // +1 if fp2 is a substructure of fp1
    public static int compare (int[] fp1, int[] fp2) {
	if (fp1.length != fp2.length) {
	    throw new IllegalArgumentException
		("Can't compare fingerprints of different sizes");
	}
	int fwd = 0, rev = 0;
	for (int i = 0; i < fp1.length; ++i) {
	    int mask = fp1[i] & fp2[i];
	    if (mask == fp1[i]) {
		++fwd;
	    }
	    if (mask == fp2[i]) {
		++rev;
	    }
	}
	if (fwd == fp1.length) return -1;
	if (rev == fp1.length) return +1;
	return 0;
    }

    public double tanimoto (RGroupNode node) {
	return tanimoto (getFp (), node.getFp());
    }

    public static double tanimoto (int[] fp1, int[] fp2) {
	if (fp1.length != fp2.length) {
	    throw new IllegalArgumentException 
		("Can't compute tanimoto for bit arrays of different sizes");
	}
	int Nc = 0, Nb = 0, Na = 0;
	for (int i = 0; i < fp1.length; ++i) {
	    Nc += ChemUtil.countBits(fp1[i] & fp2[i]);
	    Na += ChemUtil.countBits(fp1[i]);
	    Nb += ChemUtil.countBits(fp2[i]);
	}
	return (double)Nc/(Na + Nb - Nc);
    }

    public String toTreeML () {
	StringBuffer s = new StringBuffer ();
	
	String pad = "";
	for (int i = 0; i <= getLevel (); ++i) {
	    pad += " ";
	}

	s.append(pad);
	if (isLeaf ()) {
	    s.append("<leaf>\n");
	    s.append(pad + " <attribute name=\"structure\" value=\"" 
		     + getFragmentAsSmiles () + "\"/>\n");
	    s.append(pad + " <attribute name=\"name\" value=\"" 
		      + getName() + "\"/>\n");
	    s.append(pad + "</leaf>");
	}
	else {
	    s.append("<branch>\n");
	    s.append(pad + " <attribute name=\"name\" value=\"" 
		     + getName () + "\"/>\n");
	    s.append(pad + " <attribute name=\"structure\" value=\"" 
		     + rgMol.toFormat("cxsmarts") + "\"/>\n");
	    for (RGroupNode child = (RGroupNode)getFirstChild();
		 child != null; child = (RGroupNode)child.getNextSibling()) {
		s.append(child.toTreeML() + "\n");
	    }
	    s.append(pad + "</branch>");
	}
	return s.toString();
    }

    public static void printHierarchy (PrintStream ps, RGroupNode node) {
	Molecule m = node.getFragment();
	int depth = node.getLevel();
	if (m != null) {
	    for (int i = 0; i < depth; ++i) {
		ps.print(" ");
	    }

	    ps.print(depth + ". " + m.getName() 
			     + " " + m.toFormat("smiles:u0"));
	    if (!node.isLeaf()) {
		ps.println(" R-group: " + node.rgMol.toFormat("cxsmarts"));
	    }
	    else {
		ps.println();
		RGroupNode parent = (RGroupNode)node.getParent();
		Molecule[] rgs = parent.getLigands(node);
		if (rgs != null) {
		    for (int r = 0; r < rgs.length; ++r) {
			if (rgs[r] != null) {
			    for (int i = 0; i < depth; ++i) {
				ps.print(" ");
			    }
			    ps.println("    R" + parent.rmap[r] + ". " 
				       + rgs[r].toFormat("smiles:u0"));
			}
		    }
		}
	    }
	}

	if (!node.isLeaf()) {
	    for (RGroupNode child = (RGroupNode)node.getFirstChild();
		 child != null; child = (RGroupNode)child.getNextSibling()) {
		printHierarchy (ps, child);
	    }
	}
    }

    public static void main (String[] argv) throws Exception {
	MolImporter mi = null;
	if (argv.length == 0) {
	    System.err.println("** Reading from stdin");
	    mi = new MolImporter (System.in);
	}
	else {
	    mi = new MolImporter (argv[0]);
	}

	RGroupNode node = new RGroupNode ();
	Vector<RGroupNode> molv = new Vector<RGroupNode>();
	for (Molecule mol; (mol = mi.read()) != null; ) {
	    mol.aromatize(Molecule.AROM_BASIC);
	    mol.calcHybridization();
	    mol.hydrogenize(false);

	    System.out.println("..adding " + mol.getName());
	    node.add(mol);

	    RGroupNode root = (RGroupNode)node.getRoot();
	    if (root != node) {
		node = root;
		System.out.print("## new root: " + root.getName());
		if (root.getFragment() != null) {
		    System.out.print(" " + root.getFragment().toFormat
				     ("smiles:u0-a"));
		}
		System.out.println();
	    }

	    /*
	    System.out.println("** " + mol.getName());
	    RGroupNode n = new RGroupNode (mol);
	    for (RGroupNode m : molv) {
		System.out.println("  " + String.format("%1$.3f ", 
							n.tanimoto(m))
				   + m.getName());
	    }
	    molv.add(n);
	    */
	}
	mi.close();

	printHierarchy (System.out, node);

	System.out.println("** pruned hierarchy **");
	pruneChildren (node);
	printHierarchy (System.out, node);
    }
}
