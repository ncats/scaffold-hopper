// $Id: RGroupTreeModel.java 2434 2008-11-18 15:45:12Z nguyenda $

package gov.nih.ncgc.rgroup;

import java.util.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.tree.TreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.DefaultTreeModel;

import chemaxon.struc.Molecule;
import chemaxon.struc.RgMolecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.MolSearch;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.RGroupDecomposition;

public class RGroupTreeModel extends DefaultTreeModel {
    private RGroupNode rootNode = new RGroupNode ();

    public RGroupTreeModel () {
	super (null);
    }

    public void add (Molecule mol) {
	//mol.aromatize(Molecule.AROM_BASIC);
	mol.aromatize();
	mol.calcHybridization();
	mol.hydrogenize(false);

	rootNode.add(mol);
	RGroupNode newRoot = (RGroupNode)rootNode.getRoot();
	if (newRoot != rootNode) {
	    // update the root node
	    setRoot (rootNode = newRoot);
	}
    }

    public void buildRgroup () {
	rootNode.prune();

	// now walk the tree and build R-group for each non-leaf node
	buildRgroup (rootNode);
    }

    protected static void buildRgroup (RGroupNode node) {
	if (!node.isLeaf()) {
	    node.buildRgroup();

	    for (RGroupNode child = (RGroupNode)node.getFirstChild();
		 child != null; ) {
		buildRgroup (child);
		child = (RGroupNode) child.getNextSibling();
	    }
	}
    }

    public String toString () {
	ByteArrayOutputStream bos = new ByteArrayOutputStream ();
	RGroupNode.printHierarchy(new PrintStream (bos), rootNode);
	return new String (bos.toByteArray());
    }

    public void printTreeML (PrintStream ps) {
	ps.println("<?xml version=\"1.0\"?>");
	ps.println("<tree>");
	ps.println(" <declarations>");
	ps.println("   <attributeDecl name=\"structure\" type=\"String\"/>");
	ps.println("   <attributeDecl name=\"name\" type=\"String\"/>");
	ps.println("  </declarations>");
	ps.println(rootNode.toTreeML());
	ps.println("</tree>");
    }

    public static void main (String[] argv) throws Exception {
	MolImporter mi = null;
	String out = null;
	if (argv.length == 0) {
	    System.err.println("** Reading from stdin");
	    mi = new MolImporter (System.in);
	}
	else {
	    mi = new MolImporter (argv[0]);
	    if (argv.length > 1) {
		out = argv[1];
	    }
	}

	RGroupTreeModel model = new RGroupTreeModel ();

	MolHandler mh = new MolHandler ("Cc1cnc2ccccc2c1");
	mh.aromatize();
	Molecule core = mh.getMolecule();
	core.setName("Core");
	//model.add(core.cloneMolecule());

	RgMolecule rgmol = new RgMolecule ();
	rgmol.setRoot(core);
	RGroupDecomposition.addRGroups(rgmol);

	RGroupDecomposition rgd = new RGroupDecomposition();
	rgd.setAttachmentType(RGroupDecomposition.ATTACHMENT_MAP);
	rgd.setQuery(rgmol);

	MolAtom[] atoms = rgmol.getAtomArray();
	int[] gi = new int[atoms.length];
	rgmol.getGrinv(gi);
	for (int i = 0; i < atoms.length; ++i) {
	    int atno = atoms[i].getAtno();
	    if (atno == MolAtom.RGROUP) {
		System.out.println
		    (gi[i] + " R" + atoms[i].getRgroup() + " at atom " + (i+1));
	    }
	}

	// rmap: the i-th element is the r-group index of the i-th 
	// query atom, 0 for scaffold atoms
	int[] rmap = rgd.getQueryRMap();
	Map<Integer,Set<String>> rgroups 
	    = new TreeMap<Integer, Set<String>>();

	System.out.println("** Rgroup: " + rgmol.toFormat("cxsmarts"));

	for (Molecule mol; (mol = mi.read()) != null; ) {
	    model.add(mol);

	    rgd.setTarget(mol);
	    int[] hit  = rgd.findFirst();
	    if (hit != null) {
		System.out.println("## " + mol.getName());
		Molecule[] ligands = rgd.findLigands(hit);
		for (int i = 0; i < ligands.length; ++i) {
		    if (rmap[i] != 0) {
			System.out.println
			    ("** R" + rmap[i] + ": " 
			     + ligands[i].toFormat("smiles:u0"));
			Set<String> rg = rgroups.get(rmap[i]);
			if (rg == null) {
			    rgroups.put(rmap[i], rg = new HashSet<String>());
			}
			rg.add(ligands[i].toFormat("smiles:u0"));
		    }
		}
		System.out.println();
	    }
	}

	//System.out.println("** RGroup hierarchy");
	//System.out.println(model);

	model.buildRgroup();

	if (out != null) {
	    System.out.println("** writing tree to file " +out);
	    FileOutputStream fos = new FileOutputStream (out);
	    model.printTreeML(new PrintStream (fos));
	    fos.close();
	}
	else {
	    System.out.println(model);
	}

	/*
	System.out.println("** pruned RGroup hierarchy");
	System.out.println(model);

	for (Map.Entry<Integer, Set<String>> e : rgroups.entrySet()) {
	    Set<String> rg = e.getValue();
	    if (rg.size() < 2) { // ignore

	    }
	    else {
		System.out.println("R" + e.getKey() + ": " + rg);
	    }
	}
	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom a = atoms[i];
	    
	    if (a.getAtno() == MolAtom.RGROUP) {
		Set<String> rg = rgroups.get(a.getRgroup());
		if (rg != null && rg.size() < 2) {
		    // remove this R-group
		    System.out.println("** removing R" + a.getRgroup());
		    rgmol.removeNode(a);
		}
	    }
	}
	System.out.println
	    ("## pruned Rgroup: " + rgmol.toFormat("cxsmarts:u"));
	*/
    }
}
