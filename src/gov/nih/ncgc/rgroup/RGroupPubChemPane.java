package gov.nih.ncgc.rgroup;

import java.beans.*;
import java.io.*;
import java.net.*;
import javax.net.ssl.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.locks.ReentrantLock;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import javax.swing.text.AttributeSet;
import javax.swing.text.html.HTML;

import org.jdesktop.jxlayer.JXLayer;
import org.jdesktop.jxlayer.plaf.LayerUI;
import org.jdesktop.jxlayer.plaf.AbstractLayerUI;
import org.jdesktop.swingworker.SwingWorker;
import org.jdesktop.swingx.JXTaskPane;
import org.jdesktop.swingx.JXTaskPaneContainer;
import org.jdesktop.swingx.JXErrorPane;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.formats.MolImporter;

import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.node.ArrayNode;
import org.codehaus.jackson.JsonNode;

import gov.nih.ncgc.viz.LabelRenderer;
import gov.nih.ncgc.util.DummySSLSocketFactory;

public class RGroupPubChemPane extends JPanel {
    static final Logger logger = 
        Logger.getLogger(RGroupPubChemPane.class.getName());

    /**
     * PubChem search server
     */
    static final String PCS_URL = "https://tripod.nih.gov/pcs/search?";
    static final String BARD_URL = "https://bard.nih.gov/api/v17.3/compounds/";


    class BARDLoader extends SwingWorker<Void, Void> {
        PubChemTableModel model;
        BARDLoader (PubChemTableModel model) {
            this.model = model;
        }

        @Override
        protected Void doInBackground () {
            setMessage ("Loading compound info...");
            model.instruments();
            return null;
        }

        @Override
        protected void done () {
            clearMessage ();
            try {
                get ();
            }
            catch (Exception ex) {
                JXErrorPane.showDialog(ex);
            }
        }
    }

    class PCSLoader extends SwingWorker<Throwable, Molecule> {
        Molecule query;
        PubChemTableModel model = (PubChemTableModel)pcTable.getModel();

        PCSLoader (Molecule query) {
            this.query = query;
        }
        
        @Override
        protected Throwable doInBackground () {
            setMessage ("Searching PubChem...");
            model.clear();
            try {
                URL url  = new URL (PCS_URL+"s="
                                    +URLEncoder.encode(getQuery (query), "utf-8")
                                    +"&mode=text&format=smiles");
                HttpsURLConnection con =
                    (HttpsURLConnection)url.openConnection();
                con.setSSLSocketFactory(new DummySSLSocketFactory ());
                
                logger.info("Loading "+url);
                MolImporter mi = new MolImporter (con.getInputStream());

                ArrayList<Molecule> queue = new ArrayList<Molecule>();
                for (Molecule mol; (mol = mi.read()) != null; ) {
                    if (Thread.currentThread().isInterrupted()) {
                        logger.warning(Thread.currentThread()+" interrupted!");
                        return null;
                    }

                    for (MolAtom a : mol.getAtomArray()) {
                        int map = a.getAtomMap();
                        if (map > 0) {
                            // highlighting
                            a.setSetSeq(3);
                            a.setAtomMap(0);
                        }
                    }

                    queue.add(mol);
                }

                Collections.sort
                    (queue, new Comparator<Molecule>() {
                        public int compare (Molecule m1, Molecule m2) {
                            Long l1 = Long.parseLong(m1.getName());
                            Long l2 = Long.parseLong(m2.getName());
                            if (l1 < l2) return -1;
                            if (l1 > l2) return 1;
                            return 0;
                        }
                    });

                publish (queue.toArray(new Molecule[0]));
            }
            catch (Exception ex) {
                if (ex.getMessage().indexOf("empty") > 0)
                    ; // nothing is return, so not an error
                else
                    return ex;
            }
            return null;
        }

        @Override
        protected void process (Molecule... mols) {
            for (Molecule m : mols)
                model.add(m);
        }

        @Override
        protected void done () {
            clearMessage ();
            try {
                Throwable t = get ();
                if (t != null) {
                    JXErrorPane.showDialog(t);
                }
                else {
                    model.fireTableDataChanged();
                    firePubChemPropertyChange ();
                    future = threadPool.submit(new BARDLoader (model));
                }
            }
            catch (Exception ex) {
                JXErrorPane.showDialog(ex);
            }
        }
    }

    class PubChemTableModel extends AbstractTableModel {
        String[] columns = new String[]{
            "Structure",
            "CID",
            "Name",
            "Class",
            "Assays",
            "Actives",
            "MolWt",
            "TPSA",
            "XLOGP",
            "Complexity",
            "Rotatable",
            "Donor",
            "Acceptor",
            "IUPAC"
        };
        ArrayList<Molecule> mols = new ArrayList<Molecule>();
        AtomicBoolean instrumented = new AtomicBoolean (false);
        final ReentrantLock lock = new ReentrantLock ();

        PubChemTableModel () {
        }

        public void clear () { 
            instrumented.set(false);
            mols.clear(); 
            fireTableStructureChanged ();
        }
        public void add (Molecule m) {
            instrumented.set(false);
            mols.add(m);
        }
        public Collection<Molecule> getMols () { return mols; }

        public int getRowCount () { return mols.size(); }
        public int getColumnCount () { return columns.length; }
        public Class getColumnClass (int col) {
            if (col == 0) return Molecule.class;
            if (col == 1) return Long.class;

            String p = columns[col];
            if (p.equalsIgnoreCase("xlogp") 
                || p.equalsIgnoreCase("molwt")
                || p.equalsIgnoreCase("tpsa"))
                return Double.class;
            else if (p.equalsIgnoreCase("complexity")
                     || p.equalsIgnoreCase("rotatable")
                     || p.equalsIgnoreCase("donor")
                     || p.equalsIgnoreCase("acceptor")
                     || p.equalsIgnoreCase("assays")
                     || p.equalsIgnoreCase("actives"))
                return Integer.class;
            return String.class;
        }

        @Override
        public String getColumnName (int col) { return columns[col]; }

        public void instruments () {
            if (instrumented.get())
                return;

            lock.lock();
            try {
                ObjectMapper mapper = new ObjectMapper ();
                for (Molecule mol : mols) {
                    try {
                        if (Thread.currentThread().isInterrupted())
                            return;

                        URL u = new URL (BARD_URL+mol.getName());
                        HttpsURLConnection con =
                            (HttpsURLConnection)u.openConnection();
                        con.setSSLSocketFactory(new DummySSLSocketFactory ());
                        
                        //logger.info("Fetching "+u+"...");

                        setMessage ("Fetching "+mol.getName()+"...");
                        JsonNode node = mapper.readTree(con.getInputStream());
                        if (node != null) {
                            if (node.isArray()) {
                                node = node.get(0);
                            }
                            mol.setProperty("Name", asText (node, "name"));
                            mol.setProperty
                                ("IUPAC", asText (node, "iupacName"));
                            mol.setProperty("MolWt", asText (node, "mwt"));
                            mol.setProperty("TPSA", asText (node, "tpsa"));
                            mol.setProperty("XLOGP", asText (node, "xlogp"));
                            mol.setProperty("Complexity", 
                                            asText (node, "complexity"));
                            mol.setProperty("Rotatable", 
                                            asText (node, "rotatable"));
                            mol.setProperty("Donor", 
                                            asText (node, "hbondDonor"));
                            mol.setProperty("Acceptor",
                                            asText (node, "hbondAcceptor"));
                            mol.setProperty("Class", 
                                            asText (node, "compoundClass"));
                            mol.setProperty("Assays", 
                                            asText (node, "numAssay"));
                            mol.setProperty("Actives",
                                            asText (node, "numActiveAssay"));
                        }

                        fireTableDataChanged ();
                    }
                    catch (Exception ex) {
                        logger.warning("Can't retrieve compound info for "
                                       +mol.getName());
                    }
                }
                
                instrumented.set(true);
            }
            finally {
                lock.unlock();
            }
        }

        public Object getValueAt (int row, int col) {
            Molecule m = mols.get(row);
            if (m != null) {
                if (col == 0) return m;
                if (col == 1) return Long.parseLong(m.getName());

                //if (lock.tryLock()) {
                    try {
                        String p = columns[col];
                        String value = m.getProperty(p);
                        if (value == null) {
                        }
                        else if (p.equalsIgnoreCase("xlogp") 
                                 || p.equalsIgnoreCase("molwt")
                                 || p.equalsIgnoreCase("tpsa"))
                            return Double.parseDouble(value);
                        else if (p.equalsIgnoreCase("complexity")
                                 || p.equalsIgnoreCase("rotatable")
                                 || p.equalsIgnoreCase("donor")
                                 || p.equalsIgnoreCase("acceptor")
                                 || p.equalsIgnoreCase("assays")
                                 || p.equalsIgnoreCase("actives"))
                            return Integer.parseInt(value);
                        return value;
                    }
                    finally {
                        //lock.unlock();
                    }
                    //}
            }
            
            return null;
        }
    }

    volatile String message;
    volatile Future future; // current running task
    JTable pcTable;
    PropertyChangeSupport pcs = new PropertyChangeSupport (this);
    ExecutorService threadPool = Executors.newSingleThreadExecutor();

    public RGroupPubChemPane () {
        super (new BorderLayout ());

        pcTable = new RTable (new PubChemTableModel ());

        JPanel panel = new JPanel (new BorderLayout ());
        panel.add(new JScrollPane (pcTable));

        final LabelRenderer labeler = new LabelRenderer ();
        labeler.setLabelFilled(false);
        JXLayer layer = new JXLayer (panel, new AbstractLayerUI () {
                @Override
                public long getLayerEventMask() { 
                    /*
                     * this prevents the layer from intercepting the input 
                     * masks such as keyboard and mouse events.  we have to
                     * do this because turning on the masks will mess around 
                     * with the AWT's default toolkit event listener that 
                     * prevents it to run within webstart sandbox.
                     */
                    return 0; 
                }

                @Override
                protected void paintLayer (Graphics2D g2, JXLayer layer) {
                    super.paintLayer(g2, layer);
                    String mesg = getMessage ();
                    if (mesg != null) {
                        Rectangle r = layer.getBounds();
                        labeler.renderLabel(g2, r.height/2, 
                                            5, 5, r.width, r.height, 
                                            mesg);
                    }
                }
            });
        add (layer);
    }

    public void setRowHeight (int height) {
        pcTable.setRowHeight(height);
    }

    public void setMessage (String message) {
        this.message = message;
        repaint ();
    }
    public void clearMessage () { 
        setMessage (null); 
    }
    public String getMessage () { return message; }

    public void search (Molecule scaffold) {
        load (new PCSLoader (scaffold));
    }

    public int getCount () {
        return pcTable.getRowCount();
    }

    protected void load (PCSLoader loader) {
        if (future != null) {
            future.cancel(true);
        }
        future = threadPool.submit(loader);
    }

    protected void firePubChemPropertyChange () {
        pcs.firePropertyChange(new PropertyChangeEvent 
                               (this, "pubchem", null, 
                                ((PubChemTableModel)pcTable
                                 .getModel()).getMols().toArray
                                (new Molecule[0])));
    }

    public void addPropertyChangeListener (PropertyChangeListener l) {
        pcs.addPropertyChangeListener(l);
    }
    public void removePropertyChangeListener (PropertyChangeListener l) {
        pcs.removePropertyChangeListener(l);
    }

    static String getQuery (Molecule mol) {
        String smarts = mol.toFormat("smarts");
        StringBuilder query = new StringBuilder ();
        for (int i = 0; i < smarts.length(); ++i) {
            char ch = smarts.charAt(i);
            if (ch == '/' || ch == '\\') {
                // ignore e/z
            }
            else {
                query.append(ch);
            }
        }
        return query.toString();
    }

    static String asText (JsonNode node, String name) {
        JsonNode n = node.get(name);
        return n.isNull() ? null : n.asText();
    }
}
