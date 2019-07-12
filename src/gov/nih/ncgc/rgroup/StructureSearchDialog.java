
package gov.nih.ncgc.rgroup;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.*;
import org.jdesktop.swingworker.SwingWorker;

import chemaxon.marvin.beans.MViewPane;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.*;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolFormatException;
import chemaxon.util.MolHandler;
import gov.nih.ncgc.search.SearchParams;


/**
 * subclass must implement the search () method!
 */
public abstract class StructureSearchDialog extends JDialog 
    implements ActionListener, PropertyChangeListener {

    private static final Logger logger = Logger.getLogger
	(StructureSearchDialog.class.getName());
    
    private MSketchPane msketch;
    private JRadioButton subRb, simRb, exactRb, supRb;
    private JTextField tanimotoField;
    private JButton searchBtn, cancelBtn;
    protected SearchParams query;

    public StructureSearchDialog (Frame owner) {
	super (owner, false);
	initDialog ();
	setTitle ("Structure Search");
    }

    protected void initDialog () {
	JPanel pane = new JPanel (new BorderLayout (0, 5));

	final JLabel tanLabel = new JLabel ("Tanimoto");
	tanLabel.setEnabled(false);
	ActionListener action = new ActionListener () {
		public void actionPerformed (ActionEvent e) {
		    Object source = e.getSource();
		    if (source == subRb
			|| source == exactRb
			|| source == supRb) {
			tanimotoField.setEnabled(false);
			tanLabel.setEnabled(false);
		    }
		    else if (source == simRb) {
			tanimotoField.setEnabled(true);
			tanLabel.setEnabled(true);
		    }
		}
	    };

	msketch = new MSketchPane ();
	msketch.setCloseEnabled(false);
	msketch.addPropertyChangeListener(this);

        JMenuBar menubar = msketch.getJMenuBar();
        ArrayList<JMenu> remove = new ArrayList<JMenu>();
        for (int i = 0; i < menubar.getMenuCount(); ++i) {
            JMenu menu = menubar.getMenu(i);
            if (menu.getText().equals("Help") 
                || menu.getText().equals("Tools")) {
                remove.add(menu);
            }
        }
        for (JMenu m : remove) {
            menubar.remove(m);
        }

	pane.add(msketch);

	Box box = Box.createHorizontalBox();
	ButtonGroup bg = new ButtonGroup ();
	box.add(subRb = new JRadioButton ("Substructure"));
	subRb.setToolTipText
	    ("Search for references with the specified substructure");
	subRb.addActionListener(action);
	subRb.setSelected(true);
	bg.add(subRb);
	box.add(Box.createHorizontalStrut(10));
        box.add(Box.createHorizontalGlue());

	box.add(supRb = new JRadioButton ("Superstructure"));
	supRb.setToolTipText
	    ("Search for all references containg compounds that "
             +"are superstructure of the query");
	supRb.addActionListener(action);
	bg.add(supRb);
	box.add(Box.createHorizontalStrut(10));
        box.add(Box.createHorizontalGlue());

	box.add(exactRb = new JRadioButton ("Exact"));
	exactRb.setToolTipText
            ("Search for references containg exact compound");
	exactRb.addActionListener(action);
	bg.add(exactRb);
	box.add(Box.createHorizontalStrut(10));
        box.add(Box.createHorizontalGlue());

	box.add(simRb = new JRadioButton ("Similarity"));
	simRb.setToolTipText
	    ("Search for references with similar "
             +"compounds (in the Tanimoto sense)");
	simRb.addActionListener(action);
	bg.add(simRb);

	box.add(Box.createHorizontalStrut(5));
	box.add(tanimotoField = new JTextField (5));
	tanimotoField.setEnabled(false);
	tanimotoField.setText("0.85");
	tanimotoField.setToolTipText("Specify minimum Tanimoto cutoff");
	pane.add(box, BorderLayout.SOUTH);

	JPanel control = new JPanel (new BorderLayout ());
	control.setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
	control.add(pane);

	JPanel bpane = new JPanel (new GridLayout (1, 2, 2, 0));
	JButton btn;
	bpane.add(btn = new JButton ("Search"));
	btn.addActionListener(this);
	btn.setEnabled(false);
	searchBtn = btn;
	bpane.add(btn = new JButton ("Close"));
	btn.addActionListener(this);

	JPanel cpane = new JPanel ();
	cpane.add(bpane);

	control.add(cpane, BorderLayout.SOUTH);
	getContentPane().add(control);

	pack ();
    }

    abstract protected void search (Molecule query, SearchParams params);


    public void actionPerformed (ActionEvent e) {
	String cmd = e.getActionCommand();

	if (cmd.equalsIgnoreCase("close")) {
	    setVisible (false);
	    return;
	}

	SearchParams params = null;
        Molecule mol = msketch.getMol();
        mol.aromatize();
	if (subRb.isSelected()) {
            params = SearchParams.substructure();
	}
        else if (exactRb.isSelected()) {
            params = SearchParams.exact();
        }
        else if (supRb.isSelected()) {
            params = SearchParams.superstructure();
        }
	else {
	    try {
		for (MolAtom a : mol.getAtomArray()) {
		    int atno = a.getAtno();
		    if (a.isQuery() 
			|| atno == MolAtom.LIST 
			|| atno == MolAtom.ANY) {
			throw new Exception ("Query molecule");
		    }
		}

		for (MolBond b : mol.getBondArray()) {
		    if (b.isQuery() || b.getType() == MolBond.ANY) {
			throw new Exception ("Query molecule");
		    }
		}
	    }
	    catch (Exception ex) {
		JOptionPane.showMessageDialog
		    (this, "Markush structure not supported "
		     +"for the selected search type!", "Error", 
		     JOptionPane.ERROR_MESSAGE);
		return;
	    }

	    if (simRb.isSelected()) {
		String sim = tanimotoField.getText();
		try {
                    params = SearchParams.similarity(Double.parseDouble(sim));
		}
		catch (NumberFormatException ex) {
		    JOptionPane.showMessageDialog
			(this, "Invalid Tanimoto value: "+sim, 
			 "Error", JOptionPane.ERROR_MESSAGE);
		    return;
		}
	    }
	}
        setVisible (false);

        // now execute the search
        search (mol, params);
    }

    public void propertyChange (PropertyChangeEvent e) {
	String name = e.getPropertyName();
	//logger.info(name+": old="+e.getOldValue()+" new="+e.getNewValue());
        /*
	if ("doc0".equals(name)) {
	    statusField.setText(null);
	    searchBtn.setEnabled(mview.getM(0).getAtomCount() > 0);
	}
        */
        if ("mol".equals(name)) {
	    searchBtn.setEnabled(msketch.getMol().getAtomCount() > 0); 
        }
    }

    public static void main (String[] argv) throws Exception {
	StructureSearchDialog sd = new StructureSearchDialog (null) {
                protected void search (Molecule query, SearchParams params) {
                    System.out.println(params+" "+query.toFormat("cxsmarts"));
                }
            };
	sd.setVisible(true);
    }
}
