// $Id: RGroupUtil.java 3549 2009-11-12 23:02:59Z nguyenda $
package gov.nih.ncgc.rgroup;

import nu.xom.Attribute;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Serializer;

import java.io.*;
import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;

public class RGroupUtil {
    static final String LIGAND_FORMAT = "cxsmiles";

    // not meant to be instantiate
    private RGroupUtil () { 
    }

    public static Element getScaffoldsAsXML 
	(RGroupGenerator rgroup, String format) {
        Element root = new Element("scaffoldList");
        // print each scaffold and its members
        for (int i = 0; i < rgroup.getScaffoldCount(); ++i) {
            RGroupTable rtab = rgroup.getRGroupTable(i);

            Element scaffold = new Element("scaffold");
            Attribute attr = new Attribute
		("core", rtab.getScaffold().toFormat(format));
            scaffold.addAttribute(attr);

	    // first column is ID
	    // second column is Structure
	    // ... rest R's
            for (int r = 0; r < rtab.getRowCount(); ++r) {
                Molecule m = (Molecule) rtab.getValueAt(r, 1);

                Element member = new Element("member");
                Attribute molid = new Attribute("id", m.getName());
                member.addAttribute(molid);

                // now the rest of the columns are r-group ligands
                for (int c = 2; c < rtab.getColumnCount(); ++c) {
                    Molecule ligand = (Molecule) rtab.getValueAt(r, c);
                    if (ligand != null) {
                        Element relem = new Element("ligand");
                        Attribute ligandattr = new Attribute("id", "R"+(c-1));
                        relem.addAttribute(ligandattr);
                        relem.appendChild(ligand.toFormat(format));
                        member.appendChild(relem);
                    }
                }
                scaffold.appendChild(member);
            }
            root.appendChild(scaffold);
        }
        return root;
    }

    public static void exportXML (OutputStream os, RGroupGenerator rgroup) 
	throws IOException {
        // dump the document
        Document doc = new Document 
	    (getScaffoldsAsXML (rgroup, LIGAND_FORMAT));
	Serializer serializer = new Serializer (os, "ISO-8859-1");
	serializer.setIndent(2);
	serializer.setMaxLength(512);
	serializer.write(doc);
    }

    public static void main (String[] argv) throws Exception {
        RGroupGenerator rgroup = new RGroupGenerator ();
        for (int i = 0; i < argv.length; ++i) {
            MolImporter molimp = new MolImporter(argv[i]);
            for (Molecule mol; (mol = molimp.read()) != null;) {
                String name = mol.getName();
                if (name == null || name.equals("")) {
                    for (int j = 0; j < mol.getPropertyCount(); ++j) {
                        name = mol.getProperty(mol.getPropertyKey(j));
                        if (name != null && !name.equals("")) {
                            break;
                        }
                    }
                    mol.setName(name);
                }
                rgroup.add(mol);
            }
            molimp.close();
        }
        rgroup.run();
	exportXML (System.out, rgroup);
    }
}
