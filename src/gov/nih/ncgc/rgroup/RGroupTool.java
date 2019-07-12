// $Id: RGroupTool.java 3864 2009-12-18 22:09:45Z nguyenda $
package gov.nih.ncgc.rgroup;

import chemaxon.formats.MolImporter;
import chemaxon.marvin.beans.MViewPane;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import com.jgoodies.looks.plastic.Plastic3DLookAndFeel;
import com.jidesoft.swing.CheckBoxList;
import com.jidesoft.swing.JideTabbedPane;
import gov.nih.ncgc.descriptor.TopologicalIndices;
import gov.nih.ncgc.model.DataSeq;
import gov.nih.ncgc.search.SearchParams;
import gov.nih.ncgc.util.MolFpFactory;
import gov.nih.ncgc.util.MolStandardizer;
import gov.nih.ncgc.util.DummySSLSocketFactory;

import org.codehaus.jackson.JsonNode;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.node.ArrayNode;
import org.codehaus.jackson.node.ObjectNode;
import org.jdesktop.swingworker.SwingWorker;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;

import javax.imageio.ImageIO;
import javax.jnlp.FileContents;
import javax.jnlp.FileOpenService;
import javax.jnlp.FileSaveService;
import javax.jnlp.ServiceManager;
import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.plaf.UIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.net.*;
import javax.net.ssl.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class RGroupTool extends JFrame 
    implements PropertyChangeListener, 
               ListSelectionListener,
               ActionListener {

    static final Logger logger = Logger.getLogger(RGroupTool.class.getName());
    static final Border EMPTY = BorderFactory.createEmptyBorder(1,1,1,1);
    static boolean UnsignedWebstart = false;

    static final ImageIcon BACK_ICON =
        new ImageIcon (RGroupTool.class.getResource
                       ("/resources/arrow-180-medium.png"));
    static final ImageIcon FORWARD_ICON = 
        new ImageIcon (RGroupTool.class.getResource
                       ("/resources/arrow-000-medium.png"));
    static final ImageIcon OPEN_ICON = 
        new ImageIcon (RGroupTool.class.getResource
                       ("/resources/folder-open.png"));
    static final ImageIcon COMPOUND_ICON = 
        new ImageIcon (RGroupTool.class.getResource
                       ("/resources/pccompound_small.png"));

    static final Icon PLUS_ICON = new ImageIcon 
        (RGroupTool.class.getResource("/resources/icons/add_obj.gif"));

    /**
     * content tab
     */
    static final int REFERENCE_TAB = 0;
    static final int COMPOUND_TAB = 1;
    static final int RGROUP_TAB = 2;
    static final int PUBCHEM_TAB = 3;

    /**
     * navigation tab
     */
    static final int SCAFFOLD_TAB = 0;
    static final int NETWORK_TAB = 1;
    static final int SINGLETON_TAB = 2;

    JFileChooser chooser;
    JTable scaffoldTab;
    JTable instanceTab;
    JLabel scaffoldLabel;
    JTable rgroupTab;
    JTable singletonTab;
    JTabbedPane contentTab;
    JProgressBar progress;
    MViewPane mview;
    JTextField statusField;
    JPopupMenu popup;
    JButton extBtn;
    RGroupNetwork network;
    RGroupDocsPane docsPane;
    RGroupPubChemPane pubchemPane;
    JideTabbedPane navtab;
    JButton navBackBtn, navForwardBtn;
    JTextField searchField;
    StructureSearchDialog _strucSearchDialog;

    CheckBoxList propCheckBoxList;
    JDialog propDialog;
    Map<String, Class> props = new TreeMap<String, Class>();
    String activityDataFileName = null;

    RGroupWebResources webservice;

    /**
     * history stack
     */
    LinkedList<RGroupGenerator> history = new LinkedList<RGroupGenerator>();
    LinkedList<RGroupGenerator> forward = new LinkedList<RGroupGenerator>();
    
    //there's got to be a better way to do this.
    Map<Integer,Integer> columnPreviousWidth = new HashMap<Integer,Integer>();



    class TabControl extends JPanel implements UIResource, ActionListener {
        JButton addBtn;

        TabControl () {
            //setOpaque (true);
            setLayout (new GridLayout (1, 1));

            addBtn = new JButton (PLUS_ICON);
            addBtn.setMargin(new Insets (0, 0, 0, 0));
            addBtn.setRolloverEnabled(true);
            addBtn.setFocusPainted(false);
            addBtn.setBorderPainted(false);
            addBtn.setToolTipText("Search for scaffolds in new dataset");
            addBtn.addActionListener(this);
            addBtn.setBackground(getBackground());
            
            add (addBtn);
        }

        public void actionPerformed (ActionEvent e) {
            RGroupGenerator rgroup = getRGroup ();
            if (rgroup == null || rgroup.getScaffoldCount() == 0) {
                JOptionPane.showMessageDialog
                    (RGroupTool.this, "No reference scaffold(s) have "
                     +"been generated.", "Error",JOptionPane.ERROR_MESSAGE);
                return;
            }

            if (UnsignedWebstart) {
                try {
                    FileOpenService fos = 
                        (FileOpenService)ServiceManager.lookup
                        ("javax.jnlp.FileOpenService");
                    FileContents contents = fos.openFileDialog
                        (".", new String[]{"sdf","mol","smi","smiles",
                                           "cml","sd"});
                    if (!contents.canRead()) {
                        throw new Exception ("Can't open file");
                    }
                    new SearchRGroup (rgroup.getRGroupTables(), 
                                      contents).execute();
                }
                catch (Exception ex) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, "Your Java Webstart doesn't have "
                         +"proper permission to open files.", 
                         "Fatal Error", JOptionPane.ERROR_MESSAGE);
                }
            }
            else {
                chooser.setDialogTitle("Please select input file");
                if (JFileChooser.APPROVE_OPTION == 
                    chooser.showOpenDialog(RGroupTool.this)) {
                    File file = chooser.getSelectedFile();
                    new SearchRGroup (rgroup.getRGroupTables(), 
                                      file).execute();
                }
            }
        }
    }

    class ExtendScaffold extends SwingWorker<RGroupTable, String> {
        RGroupTable tab;
        ExtendScaffold (RGroupTable tab) {
            this.tab = tab;
        }

        @Override
        protected RGroupTable doInBackground () {
            publish ("Extending scaffold... "); 
            return getRGroup().extendScaffold(tab);
        }

        @Override
        protected void process (String... mesg) {
            for (String m : mesg) {
                statusField.setText(m);
            }
        }

        @Override
        protected void done () {
            try {
                RGroupTable tab = get ();
                if (tab != null) {
                    updateProperties ();

                    ((AbstractTableModel)scaffoldTab.getModel())
                    .fireTableDataChanged();

                    for (int i = 1; i < contentTab.getTabCount(); ++i) {
                        SearchRGroup srg = (SearchRGroup)
                            ((JComponent)contentTab.getComponentAt(i))
                            .getClientProperty("rgroup.search");
                        srg.addScaffold(tab);
                    }

                    JOptionPane.showMessageDialog
                        (RGroupTool.this, "New scaffold added!", "Message",
                         JOptionPane.INFORMATION_MESSAGE);
                }
                else {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, 
                         "Unable to extend this scaffold any further!",
                         "Message", JOptionPane.INFORMATION_MESSAGE);
                }
                statusField.setText(null);
            }
            catch (Exception ex) {
                JOptionPane.showMessageDialog
                    (RGroupTool.this, ex, "Exception", 
                     JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace();
            }
        }
    }

    class RGroupWorker extends SwingWorker<Throwable, String> {
        Object[] argv;
        
        RGroupWorker (String... argv) {
            this.argv = new Object[argv.length];
            for (int i = 0; i < argv.length; ++i) {
                this.argv[i] = argv[i];
            }
        }

        RGroupWorker (File... files) {
            this.argv = new Object[files.length];
            for (int i = 0; i < files.length; ++i) {
                argv[i] = files[i].getPath();
            }
        }

        RGroupWorker (URL... urls) {
            this.argv = new Object[urls.length];
            for (int i = 0; i < urls.length; ++i) {
                argv[i] = urls[i];
            }
        }

        RGroupWorker (FileContents... files) {
            this.argv = new Object[files.length];
            for (int i = 0; i < files.length; ++i) {
                argv[i] = files[i];
            }       
        }

        void checkProps (Molecule mol) {
            for (int i = 0; i < mol.getPropertyCount(); ++i){
                String key = mol.getPropertyKey(i);
                if (key == null) {
                    continue;
                }

                Object val = mol.getPropertyObject(key);
                Class clazz = props.get(key);
                if (clazz == null) {
                    try {
                        double x = Double.parseDouble(val.toString());
                        if (Math.rint(x) == x) {
                            props.put(key, Long.class);
                        }
                        else {
                            props.put(key, Double.class);
                        }
                    }
                    catch (NumberFormatException ex) {
                        props.put(key, String.class);
                    }
                }
                else if (clazz == Double.class) {
                    try {
                        Double.parseDouble(val.toString());
                    }
                    catch (NumberFormatException ex) {
                        props.put(key, String.class);
                    }
                }
                else if (clazz == Long.class) {
                    try {
                        double x = Double.parseDouble(val.toString());
                        if (Math.rint(x) != x) {
                            props.put(key, Double.class);
                        }
                    }
                    catch (NumberFormatException ex) {
                        props.put(key, String.class);
                    }
                }
                else if (clazz == String.class) {
                }
            }
        }

        InputStream getInputStream (Object input) throws Exception {
            if (input instanceof FileContents) {
                return ((FileContents)input).getInputStream();
            }
            else if (input instanceof File) {
                return new FileInputStream ((File)input);
            }
            else if (input instanceof URL) {
                URL url = (URL)input;
                URLConnection con = url.openConnection();
                if ("https".equals(url.getProtocol())) {
                    ((HttpsURLConnection)con)
                        .setSSLSocketFactory(new DummySSLSocketFactory ());
                }
                return con.getInputStream();
            }
            else { // assume String
                return new FileInputStream ((String)input);
            }
        }

        @Override
        protected Throwable doInBackground () {
            //props.clear();
            ((DefaultListModel)propCheckBoxList.getModel()).clear();

            try {
                publish ("Searching for scaffolds...");
                progress.setIndeterminate(true);
                int cnt = 0;
                for (Object file : argv) {
                    MolImporter mi = new MolImporter (getInputStream (file));
                    //mi.setGrabbingEnabled(true);
                    Molecule last = null;
                    try {
                        for (Molecule mol; (mol = mi.read()) != null; ) {
                            publish ("Searching for scaffolds..."
                                     +String.format("%1$5d %2$s", 
                                                    ++cnt, mol.getName()));
                            getRGroup().add(mol);
                            last = mol;
                            checkProps (mol);
                        }
                        mi.close();
                        getRGroup().setName(getName (file));
                    }
                    catch (Exception ex) {
                        logger.log(Level.WARNING, 
                                   "Can't parse molecule:  "
                                   +mi.getGrabbedMoleculeString(), ex);
                    }
                }
                /*
                  publish (rgroup.getFragmentCount() 
                  + " possible scaffold(s) generated...");
                */
                SwingUtilities.invokeLater(new Runnable () {
                        public void run () {
                            progress.setIndeterminate(false);
                        }
                    });

                getRGroup().run();

                publish ("Done!");

                }
            catch (Exception ex) {
                return ex;
            }
            return null;
        }

        @Override
        protected void process (String... mesg) {
            for (String m : mesg) {
                statusField.setText(m);
            }
        }

        String getName (Object input) throws IOException {
            if (input instanceof FileContents) {
                return ((FileContents)input).getName();
            }
            else if (input instanceof File) {
                return ((File)input).getName();
            }
            else if (input instanceof URL) {
                return ((URL)input).getFile();
            }
            return input.toString();
        }

        @Override
        protected void done() {
            try {
                Throwable t = get();
                if (t != null) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, t.getMessage(),
                         "Fatal Exception", JOptionPane.ERROR_MESSAGE);
                    t.printStackTrace();
                    statusField.setText(t.getMessage());
                    return;
                } 

                for (Map.Entry<String, Class> e : props.entrySet()) {
                    System.out.println
                        (e.getKey() + " => " + e.getValue());
                    ((DefaultListModel)propCheckBoxList
                     .getModel()).addElement(e.getKey());
                }
                propCheckBoxList.selectAll();
                
                statusField.setText(null);
                contentTab.setTitleAt(RGROUP_TAB, "R-group");
                updateRGroup ();

                progress.setValue(0);
                progress.setString(null);
                progress.setStringPainted(false);

            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    class SearchRGroup extends SwingWorker<Throwable, String> 
        implements DataSeq<Molecule> {

        RGroupTable[] rgroups;
        Object file;

        Vector<Molecule> molseq = new Vector<Molecule>();
        Map<RGroupTable, RGroupTable> rgmap = 
            new HashMap<RGroupTable, RGroupTable>();

        DefaultTableModel singletons = new DefaultTableModel 
            (new Object[]{"ID", "Structure"}, 0) {
                public Class getColumnClass (int col) {
                    if (col == 0) return String.class;
                    return Molecule.class;
                }
            };


        SearchRGroup (RGroupTable[] rgroups, Object file) {
            this.rgroups = rgroups;
            this.file = file;
        }

        /**
         * DataSeq interface
         */
        public int size () { return molseq.size(); }
        public Molecule get (int index) { return molseq.get(index); }

        public Object getFile () { return file; }
        public RGroupTable getRGroupTable (RGroupTable scaffold) {
            return rgmap.get(scaffold);
        }
        public TableModel getSingletons () { return singletons; }

        public String getName () { 
            if (UnsignedWebstart) {
                try {
                    return ((FileContents)file).getName() ;
                }
                catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
            else {
                return ((File)file).getName();
            }
            return null;
        }

        public RGroupTable addScaffold (RGroupTable scaffold) throws Exception {
            RGroupTable rg = RGroupSolver.solve
                (scaffold.getScaffold().cloneMolecule(), this, null);
            if (rg != null) {
                rgmap.put(scaffold, rg);
            }
            logger.info(getName ()+": scaffold " 
                        + scaffold.getScaffold().toFormat("cxsmarts:q") 
                        + "..." 
                        + (rg != null ? rg.getRowCount() : 0));
            return rg;
        }

        @Override
        protected Throwable doInBackground () {
            try {
                progress.setIndeterminate(true);
                publish ("Generating r-groups for "+getName()+"...");
                MolImporter mi = new MolImporter 
                    (UnsignedWebstart ?  ((FileContents)file).getInputStream()
                     : new FileInputStream ((File)file));
                Molecule last = null;
                try {
                    MolStandardizer molstd = new MolStandardizer ();
                    Set<Molecule> singles = new HashSet<Molecule>();

                    for (Molecule mol; (mol = mi.read()) != null; ) {
                        try {
                            molstd.standardize(mol);
                            publish ("Generating r-groups for "
                                     + getName ()+"..."+mol.getName());
                        }
                        catch (Exception ex) {
                            logger.log(Level.WARNING, 
                                       "Can't standardize "+mol.getName(), 
                                       ex);
                        }
                        molseq.add(mol);
                        singles.add(mol);
                        last = mol;
                    }
                    mi.close();

                    for (int i = 0; i < rgroups.length; ++i) {
                        RGroupTable rg = addScaffold (rgroups[i]);
                        if (rg != null) {
                            for (int j = 0; j < rg.rows.length; ++j) {
                                singles.remove(molseq.get(rg.rows[j]));
                            }
                        }
                        publish (getName()
                                 +": searching for scaffold " +(i+1) +"..."
                                 + (rg != null ? rg.getRowCount():0));
                    }

                    for (Molecule mol : singles) {
                        singletons.addRow(new Object[]{mol.getName(), mol});
                    }
                }
                catch (Exception ex) {
                    logger.log(Level.WARNING, 
                               "Can't parse molecule; last molecule "
                               +"parsed is " 
                               + (last != null
                                  ?last.getName() :"null"), ex);
                }
            }
            catch (Exception ex) {
                return ex;
            }
            return null;
        }

        @Override
        protected void process (String... mesg) {
            for (String m : mesg) {
                statusField.setText(m);
            }
        }

        @Override
        protected void done() {
            try {
                Throwable t = get();
                if (t != null) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, t.getMessage(),
                         "Fatal Exception", JOptionPane.ERROR_MESSAGE);
                    t.printStackTrace();
                    statusField.setText(t.getMessage());
                }
                else {
                    statusField.setText(null);
                    EventQueue.invokeLater(new Runnable () {
                            public void run () {
                                addSearchRGroup (SearchRGroup.this);
                            }
                        });
                }
                progress.setIndeterminate(false);
                progress.setValue(0);
                progress.setString(null);
                progress.setStringPainted(false);
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    class AddScaffold extends SwingWorker<RGroupTable, String> {
        Molecule scaffold;
        AddScaffold (Molecule scaffold) {
            this.scaffold = scaffold;
        }

        @Override
        protected RGroupTable doInBackground () {
            // don't prune
            RGroupTable rtab = getRGroup().generate(scaffold, false);
            if (rtab != null) {
                updateProperties ();
                // generate new scaffold for the rest of the dataset
                for (int i = 1; i < contentTab.getTabCount(); ++i) {
                    SearchRGroup srg = (SearchRGroup)
                        ((JComponent)contentTab.getComponentAt(i))
                        .getClientProperty("rgroup.search");
                    if (srg != null) {
                        try {
                            srg.addScaffold(rtab);
                        } catch (Exception e) {
                        }
                    }
                }
            }
            return rtab;
        }

        @Override
        protected void process (String... mesg) {
            for (String m : mesg) {
                statusField.setText(m);
            }
        }
        
        @Override
        protected void done () {
            try {
                RGroupTable tab = get ();
                if (tab != null) {
                    ((AbstractTableModel)scaffoldTab.getModel())
                    .fireTableDataChanged();
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, "New scaffold added!", "Message",
                         JOptionPane.INFORMATION_MESSAGE);
                }
                else {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, 
                         "No scaffold added because scaffold either\n"
                         +"already exists or is redundant!", 
                         "Message", JOptionPane.INFORMATION_MESSAGE);
                }
                statusField.setText(null);
            }
            catch (Exception ex) {
                JOptionPane.showMessageDialog
                    (RGroupTool.this, ex, "Exception", 
                     JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace();
            }
        }
    }



    class CRCCellRenderer implements TableCellRenderer {
        ChartPanel chart;
        XYPlot xyplot;
        JPanel empty;
        java.awt.Shape shape;
        
        public CRCCellRenderer () {
            chart = new ChartPanel 
                (ChartFactory.createScatterPlot
                 (null, null, null, null, PlotOrientation.HORIZONTAL, 
                  false, false, false));
            chart.getChart().setBackgroundPaint(null);
            chart.setOpaque(false);
            chart.setBackground (new Color (0, 0, 0, 0));
            chart.setDoubleBuffered(false);
            //chart.getChart().setBorderPaint(Color.white);
            chart.getChart().getPlot().setBackgroundAlpha(.0f);
            
            xyplot = chart.getChart().getXYPlot();
            xyplot.setRangeGridlinesVisible(false);
            xyplot.setOutlineVisible(false);
            //xyplot.setDomainCrosshairVisible(false);
            xyplot.setDomainGridlinesVisible(false);
            xyplot.setRangeZeroBaselineVisible(false);
            xyplot.setDomainZeroBaselineVisible(false);
            
            XYItemRenderer renderer = new  XYLineAndShapeRenderer ();
            /*
              for (int i = 0; i < CURVE_COLOR.length; ++i) {
              renderer.setSeriesPaint(i, CURVE_COLOR[i]);
              }
            */
            xyplot.setRenderer(renderer); 
            
            shape = new java.awt.geom.Ellipse2D.Double(0., 0., 4., 4.);
            
            empty = new JPanel ();
            empty.setBackground(Color.white);
        }
        
        public Component getTableCellRendererComponent 
            (JTable table, Object value, boolean selected, boolean focus,
             int row, int col) {
            if (value == null) {
                empty.setBackground(selected ? table.getSelectionBackground()
                                    : table.getBackground());
                return empty;
            }
            
            /*
              chart.getChart().setBackgroundPaint
              (selected ? table.getSelectionBackground() 
              : table.getBackground());
            */
            
            RGroupTable.CRC[] crc = (RGroupTable.CRC[])value;
            //xyplot.getDomainAxis().setVisible(table.getRowHeight() >= 100);
            int width = table.getColumnModel().getColumn(col)
                .getPreferredWidth();
            //xyplot.getRangeAxis().setVisible(width >= 100);
            
            DefaultXYDataset dataset = new DefaultXYDataset();
            XYLineAndShapeRenderer renderer = 
                (XYLineAndShapeRenderer)xyplot.getRenderer();
            
            for (int i = 0; i < crc.length; ++i) {
                double[][] xy = crc[i].xy();
            
                int series = dataset.getSeriesCount();
                renderer.setSeriesShapesVisible(series, true);
                renderer.setSeriesLinesVisible(series, false);
                renderer.setSeriesShape(series, shape);
                double[][] yx = new double[2][];
                yx[0] = xy[1];
                yx[1] = xy[0];
                dataset.addSeries("raw-"+i, yx);
            
                series = dataset.getSeriesCount();
                renderer.setSeriesShapesVisible(series, false);
                renderer.setSeriesLinesVisible(series, true);
                xy = crc[i].hillfit();
                if (xy != null) {
                    yx = new double[2][];
                    yx[0] = xy[1];
                    yx[1] = xy[0];
                    dataset.addSeries("fit-"+i, yx);
                }
                xyplot.setDataset(dataset);
            }
            
            return chart;
        }
    }


    class InstanceTableModel extends AbstractTableModel {
        Map<String, Class> types = new HashMap<String, Class>();
        ArrayList<String> columns = new ArrayList<String>();
        RGroupGenerator rgen;

        InstanceTableModel () {
        }

        InstanceTableModel (RGroupGenerator rgen, Map<String, Class> types) {
            setProperties (types);
            setRGroupGenerator (rgen);
        }

        public void setRGroupGenerator (RGroupGenerator rgen) {
            this.rgen = rgen;
            fireTableDataChanged ();
        }

        public void setProperties (Map<String, Class> types) {
            columns.clear();
            this.types.clear();
            this.types.put("Molecule", Molecule.class);
            columns.add("Molecule");
            for (Map.Entry<String, Class> e : types.entrySet()) {
                columns.add(e.getKey());
                this.types.put(e.getKey(), e.getValue());
            }
            fireTableStructureChanged ();
        }

        public int getColumnCount () { return columns.size(); }
        public int getRowCount () { return rgen != null ? rgen.size() : 0; }
        public Class getColumnClass (int col) {
            return types.get(columns.get(col));
        }
        public String getColumnName (int col) { return columns.get(col); }
        public Object getValueAt (int row, int col) {
            Molecule mol = rgen.get(row);
            if (mol != null) {
                if (col == 0) return mol;
                String prop = columns.get(col);
                Class cls = types.get(prop);
                String value = mol.getProperty(prop);
                if (value != null) {
                    if (cls == Long.class) {
                        return Float.parseFloat(value);
                    }
                    if (Double.class == cls) {
                        return Double.parseDouble(value);
                    }
                    return value;
                }
            }
            return null;
        }

        public boolean isCellEditable (int r, int c) { 
            return getColumnClass (c) == Molecule.class; 
        }

        public void setValueAt (int row, int col, Object value) {
            if (col != 0 || value == null) return;
            Molecule mol = rgen.get(row);
            mol.clear();
            ((Molecule)value).clonecopy(mol);
        }
    }

    class ScaffoldTableModel extends AbstractTableModel {
        Map<String, Class> types = new HashMap<String, Class>();
        ArrayList<String> columns = new ArrayList<String>();
        Map<String, Map<RGroupTable, ScaffoldData>> values = 
            new TreeMap<String, Map<RGroupTable, ScaffoldData>>();
        ArrayList<RGroupTable> rgs = new ArrayList<RGroupTable>();

        public ScaffoldTableModel () {
            types.put("Scaffold", Molecule.class);
            columns.add("Scaffold");
            types.put("Scaffold Score", Double.class);
            columns.add("Scaffold Score");
            types.put("Complexity", Integer.class);
            columns.add("Complexity");
            types.put("Count", Integer.class);
            columns.add("Count");
        }

        public ScaffoldTableModel (Map<String, Class> types) {
            this.types.put("Scaffold", Molecule.class);
            columns.add("Scaffold");
            this.types.put("Scaffold Score", Double.class);
            columns.add("Scaffold Score");
            this.types.put("Complexity", Integer.class);
            columns.add("Complexity");
            this.types.put("Count", Integer.class);
            columns.add("Count");
            addProperties (types);
        }

        public void addProperties (Map<String, Class> types) {
            for (Map.Entry<String, Class> e : types.entrySet()) {
                if (e.getValue() == Double.class) {
                    Class c = this.types.put(e.getKey(), ScaffoldData.class);
                    if (c == null) {
                        columns.add(e.getKey());
                    }
                }
            }
            /*
              for (String s : columns) {
              System.out.println(s + " => " + this.types.get(s));
              }
            */
            fireTableStructureChanged ();
        }

        public void removeProperties (Collection<String> props) {
            for (String s : props) {
                Class c = types.remove(s);
                if (c != null) {
                    columns.remove(s);
                }
            }
            fireTableStructureChanged ();
        }

        public void setScaffolds (RGroupTable... scaffolds) {
            rgs.clear();
            if (scaffolds != null && scaffolds.length > 0) {
                for (RGroupTable rg : scaffolds) {
                    rgs.add(rg);
                }
            }
            fireTableDataChanged ();
        }

        public void update () {
            fireTableDataChanged ();
        }
                
        public RGroupTable getScaffold (int r) {
            return rgs.isEmpty() ? getRGroup().getRGroupTable(r) : rgs.get(r);
        }

        public boolean isCellEditable (int r, int c) { 
            return getColumnClass (c) == Molecule.class; 
        }
        public Class getColumnClass (int c) { 
            return types.get(columns.get(c)); 
        }
        public String getColumnName (int c) { return columns.get(c); }
        public int getColumnCount () { return columns.size(); }
        public int getRowCount () { 
            return rgs.isEmpty() 
                ?  getRGroup().getScaffoldCount() : rgs.size();
        }
        public Object getValueAt (int row, int col) {
            RGroupTable rtab = getScaffold (row);
            if (col == 0) {
                try {
                    return rtab.getScaffold();
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
            if (col == 1) {
                return rtab.getScaffoldScore();
            }

            if (col == 2) {
                return rtab.getScaffoldComplexity();
            }
            
            if (col == 3) {
                return rtab.getRowCount();
            }

            String prop = columns.get(col);
            Map<RGroupTable, ScaffoldData> rows = values.get(prop);
            if (rows == null) {
                rows = new HashMap<RGroupTable, ScaffoldData>();
                values.put(prop, rows);
            }

            ScaffoldData data = rows.get(rtab);
            if (data == null) {
                rows.put(rtab, data = new ScaffoldData (prop, rtab));
            }

            return data;
        }
    }

    public RGroupTool () {
        this (null);
    }

    public RGroupTool (RGroupWebResources web) {
        if (web == null) {
            try {
                URL url = new URL ("https://tripod.nih.gov/chembl/version");
                HttpsURLConnection con =
                    (HttpsURLConnection)url.openConnection();
                con.setSSLSocketFactory(new DummySSLSocketFactory ());
                BufferedReader br = new BufferedReader 
                    (new InputStreamReader (con.getInputStream()));
                web = new RGroupWebResources
                    ("https://tripod.nih.gov", br.readLine());                
            }
            catch (Exception ex) {
                ex.printStackTrace();
                logger.warning("Can't retrieve chembl version; use default");
                web = new RGroupWebResources
                    ("https://tripod.nih.gov", "chembl-v23");
            }
            logger.info("## Using "+web.getBase());
        }
        webservice = web;
        initGUI ();
    }

    void  initGUI () {
        if (!UnsignedWebstart) {
            chooser = new JFileChooser (".");
            chooser.setMultiSelectionEnabled(true);
        }

        propCheckBoxList = new CheckBoxList ();
        propCheckBoxList.setCheckBoxEnabled(true);
        propCheckBoxList.setClickInCheckBoxOnly(true);
        propCheckBoxList.setModel(new DefaultListModel ());

        propDialog = new JDialog (this, true);
        propDialog.setTitle("Add/remove data columns");
        { JPanel pane = new JPanel (new BorderLayout (5, 0));
            pane.add(new JScrollPane (propCheckBoxList));
            JPanel bp = new JPanel (new GridLayout (1, 2));
            ((GridLayout)bp.getLayout()).setHgap(5);
            JButton btn;
            bp.add(btn = new JButton ("OK"));
            btn.addActionListener(new ActionListener () {
                    public void actionPerformed (ActionEvent e) {
                        updateProperties ();
                        propDialog.setVisible(false);
                    }
                });
            bp.add(btn = new JButton ("Cancel"));
            btn.addActionListener(new ActionListener () {
                    public void actionPerformed (ActionEvent e) {
                        propDialog.setVisible(false);
                    }
                });
            JPanel bpp = new JPanel ();
            bpp.setBorder(BorderFactory.createTitledBorder(""));
            bpp.add(bp);
            pane.add(bpp, BorderLayout.SOUTH);
            propDialog.getContentPane().add(pane);
            propDialog.pack();
            propDialog.setSize(400, 500);
        }

        popup = new JPopupMenu ();
        JMenuItem item = new JMenuItem ("Extend");
        item.setToolTipText("Extend (if at all possible) selected scaffold");
        item.addActionListener(this);
        popup.add(item);

        popup.add(item = new JMenuItem ("Dump"));
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    RGroupTable rtab = getSelectedScaffold ();
                    if (rtab != null)
                        dumpXml (System.out, rtab);
                }
            });

        setJMenuBar (createMenuBar ());

        JPanel pane = new JPanel (new BorderLayout (0, 0));
        pane.setBorder(EMPTY);

        JSplitPane split = new JSplitPane (JSplitPane.HORIZONTAL_SPLIT);
        split.setDividerSize(7);
        split.setResizeWeight(.35);
        split.setLeftComponent(createNavPane ());
        split.setRightComponent(createContentPane ());

        pane.add(split);
        pane.add(createToolBar (), BorderLayout.NORTH);
        pane.add(createStatusPane (), BorderLayout.SOUTH);

        setTitle ("NCGC Scaffold Hopper ["+webservice.getBase()+"]");
        getContentPane().add(pane);
        pack ();
        setDefaultCloseOperation (EXIT_ON_CLOSE);
    }

    void updateProperties () {
        Object[] selected = propCheckBoxList.getCheckBoxListSelectedValues();
        Set<String> remove = new HashSet<String>();
        remove.addAll(props.keySet());
        for (Object keep : selected) {
            remove.remove((String)keep);
        }
        for (RGroupTable rtab : getRGroup().getRGroupTables()) {
            rtab.removeProperties(remove);
        }
        ((ScaffoldTableModel)scaffoldTab.getModel())
        .removeProperties(remove);
        
        Map<String, Class> keep = new TreeMap<String, Class>();
        for (Object obj : selected) {
            String s = (String)obj;
            keep.put(s, props.get(s));
        }
        for (RGroupTable rtab : getRGroup().getRGroupTables()) {
            rtab.addProperties(keep);
        }
        ((ScaffoldTableModel)scaffoldTab.getModel()).addProperties(keep);
    }

    Component createStatusPane () {
        JPanel pane = new JPanel (new BorderLayout (5, 0));
        pane.setBorder(BorderFactory.createCompoundBorder
                       (BorderFactory.createEtchedBorder(),
                        BorderFactory.createEmptyBorder(1,1,1,1)));

        JSlider slider = new JSlider (50, 200, 100);
        //slider.setPreferredSize(new Dimension (100, 20));
        slider.addChangeListener(new ChangeListener () {
                public void stateChanged (ChangeEvent e) {
                    JSlider slider = (JSlider)e.getSource();
                    scaffoldTab.setRowHeight(slider.getValue());
                    rgroupTab.setRowHeight(slider.getValue());
                    instanceTab.setRowHeight(slider.getValue());
                    //pubchemPane.setRowHeight(slider.getValue());
                }
            });
        pane.add(slider, BorderLayout.WEST);

        statusField = new JTextField ();
        statusField.setEditable(false);
        pane.add(statusField);

        progress = new JProgressBar (0, 100);
        progress.setBorderPainted(true);
        progress.setStringPainted(false);
        pane.add(progress, BorderLayout.EAST);

        return pane;
    }

    Component createContentPane () {
        JideTabbedPane tab = new JideTabbedPane ();
        tab.setShowCloseButton(true);
        tab.setShowCloseButtonOnTab(true);
        tab.addChangeListener(new ChangeListener () {
                public void stateChanged (ChangeEvent e) {
                    //logger.info("contentTab changed");
                }
            });

        tab.addTab("References", createDocsPane ()); // REFERENCE_TAB
        tab.addTab("Compounds", createInstancesPane ()); // COMPOUND_TAB
        tab.addTab("R-group", createRGroupPane ()); // RGROUP_TAB
        //tab.addTab("PubChem", createPubChemPane ());

        tab.setTabLeadingComponent(new TabControl ());
        tab.setTabClosableAt(0, false); // first tab can't be closed!
        tab.setTabClosableAt(1, false); 
        tab.setTabClosableAt(2, false); 
        //tab.setTabClosableAt(3, false); 

        JPanel pane = new JPanel (new BorderLayout (0, 5));
        pane.setBorder(EMPTY);
        pane.add(tab);

        contentTab = tab;
        return pane;
    }

    Component createNavPane () {
        JPanel pane = new JPanel (new BorderLayout (0, 5));
        pane.setBorder(EMPTY);
        scaffoldTab = createTable ();
        //scaffoldTab.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
        scaffoldTab.getSelectionModel()
            .setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        scaffoldTab.getSelectionModel().addListSelectionListener(this);
        scaffoldTab.addMouseListener(new MouseAdapter () {
                public void mouseClicked (MouseEvent e) {
                    if (e.isPopupTrigger()) {
                        triggerPopup (e);
                    }
                }
                public void mousePressed (MouseEvent e) {
                    if (e.isPopupTrigger()) {
                        triggerPopup (e);
                    }
                }

                void triggerPopup (MouseEvent e) {
                    popup.show(scaffoldTab, e.getX(), e.getY());
                }
            });

        JSplitPane split  = new JSplitPane (JSplitPane.VERTICAL_SPLIT);
        split.setDividerSize(7);
        split.setResizeWeight(.75);
        split.setOneTouchExpandable(true);

        JPanel scafpane = new JPanel (new BorderLayout (0, 3));
        scaffoldLabel = new JLabel 
            ("<html><b>List of scaffolds", JLabel.CENTER);

        scaffoldTab.addPropertyChangeListener(this);
        scafpane.add(scaffoldLabel, BorderLayout.NORTH);
        scafpane.add(new JScrollPane (scaffoldTab));
        split.setTopComponent(scafpane);
        split.setBottomComponent(createMViewPane ());

        navtab = new JideTabbedPane ();

        JPanel sp = new JPanel (new BorderLayout (0, 1));
        sp.add(split);
        navtab.addTab("Scaffold", sp); // SCAFFOLD_TAB
        navtab.addTab("Network", network = new RGroupNetwork ());// NETWORK_TAB
        navtab.addTab("Singleton", createSingletonPane ());// SINGLETON_TAB
        network.addPropertyChangeListener(this);
        pane.add(navtab);
        
        JPanel bp = new JPanel (new GridLayout (1, 3));
        ((GridLayout)bp.getLayout()).setHgap(2);
        JButton btn = new JButton ("Add");
        btn.addActionListener(this);
        btn.setToolTipText("Add new scaffold");
        bp.add(btn);

        btn = new JButton ("Sort");
        btn.addActionListener(this);
        btn.setToolTipText
            ("Sort the scaffolds relative to the reference structure");
        bp.add(btn);

        btn = new JButton ("Extend");
        btn.addActionListener(this);
        btn.setToolTipText
            ("Extend the selected scaffold (if possible)");
        bp.add(btn);
        extBtn = btn;

        JPanel bp2 = new JPanel ();
        bp2.add(bp);
        sp.add(bp2, BorderLayout.SOUTH);

        return pane;
    }

    public void sortSimilar(JTable jtab, final Molecule ref, int column) {
        TableModel tm = jtab.getModel();
        if (tm.getColumnClass(column).isAssignableFrom(Molecule.class)) {
            TableRowSorter trs = new TableRowSorter(tm);                        
            final MolFpFactory fpFact = MolFpFactory.getInstance();
            trs.setComparator(column, new Comparator<Molecule>() {
                                  public int compare 
                                      (Molecule a, Molecule b) {
                                      double simA = fpFact.tanimotoSim
                                          (ref, a);
                                      double simB = fpFact.tanimotoSim
                                          (ref, b);
                                      int d = 0;
                                      if (simA > simB) d = -1;
                                      else if (simA < simB) d = 1;
                                      if (d == 0) {
                                          d = b.getAtomCount() - a.getAtomCount();
                                      }
                                      if (d == 0) {
                                          d = b.getBondCount() - a.getBondCount();
                                      }
                                      if (d == 0) {
                                          d = (int)Math.signum(b.getMass() - a.getMass());
                                      }
                                      return d;
                                  }
                              });
            jtab.setRowSorter(trs);
            trs.toggleSortOrder(column);
        }
    }

    public void actionPerformed (ActionEvent e) {
        if (e.getSource() instanceof JMenuItem) {
            JMenuItem item = (JMenuItem)e.getSource();
            RGroupTable rtab = getSelectedScaffold ();
            if (rtab != null) {
                new ExtendScaffold (rtab).execute();
            }
        }
        else {
            String cmd = e.getActionCommand();
            if (cmd.equalsIgnoreCase("sort")) {
                Molecule scaffold = mview.getM(0);
                if (scaffold != null) {
                    sortSimilar(scaffoldTab,scaffold,0);
                    //rgroup.sort(scaffold);
                    //((AbstractTableModel)scaffoldTab.getModel()).fireTableDataChanged();
                    
                    
                }
            }
            else if (cmd.equalsIgnoreCase("add")) {
                if (mview.getSelectedIndex() < 0) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, 
                         "There is nothing to add!", "Message", 
                         JOptionPane.INFORMATION_MESSAGE);
                    return;
                }

                mview.applyRotationMatrices();
                Molecule scaffold = mview.getM(0);
                logger.info("Adding user-defined scaffold... "
                            +scaffold.toFormat("cxsmarts"));
                statusField.setText("Adding user-defined scaffold...");
                new AddScaffold(scaffold.cloneMolecule()).execute();
            }
            else if (cmd.equalsIgnoreCase("extend")) {
                RGroupTable rtab = getSelectedScaffold ();                
                if (rtab == null) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, 
                         "Please select a scaffold to extend", "Message",
                         JOptionPane.INFORMATION_MESSAGE);
                }
                else if (!rtab.isExtensible()) {
                    JOptionPane.showMessageDialog
                        (RGroupTool.this, 
                         "Sorry, this scaffold can't be "
                         +"extended any further!", "Message", 
                         JOptionPane.INFORMATION_MESSAGE);
                }
                else {
                    new ExtendScaffold (rtab).execute();
                }
            }
            else if (cmd.equalsIgnoreCase("search")) {
                JTextField field = (JTextField)e.getSource();
                String text = field.getText();
                if (text != null && text.length() > 0) {
                    doSearch (text.trim());
                }
                else {
                    JOptionPane.showMessageDialog
                        (this, "No search terms specified!",
                         "Error", JOptionPane.ERROR_MESSAGE);
                }
            }
            else if (cmd.equalsIgnoreCase("open")) {
                openFile ();
            }
            else if (cmd.equalsIgnoreCase("go-back")) {
                goBack ();
            }
            else if (cmd.equalsIgnoreCase("go-forward")) {
                goForward ();
            }
            else if (cmd.equalsIgnoreCase("structure-search")) {
                getStructureSearchDialog().setVisible(true);
            }
            else {
                logger.log(Level.WARNING, "unknown command: " + cmd);
            }
        }
    }

    void updateRGroup () {
        updateRGroup (getRGroup ());
    }

    void updateRGroup (RGroupGenerator rgroup) {
        instanceTab.setModel
            (new InstanceTableModel (rgroup, props)); 
        scaffoldTab.setModel(new ScaffoldTableModel (props));
        //addMolecularSorter(scaffoldTab);
        //addMolecularSorter(instanceTab);
        for (RGroupTable rtab : rgroup.getRGroupTables()) {
            rtab.addProperties(props);
        }
        network.setRGroups(rgroup.getRGroupTables());
        scaffoldLabel.setText("<html><b>"+rgroup.getName());
                
        DefaultTableModel model = new DefaultTableModel 
            (new Object[]{"ID", "Structure"}, 0) {
                public Class getColumnClass (int col) {
                    if (col == 0) return String.class;
                    return Molecule.class;
                }
            };
                
        Enumeration<Molecule> singleton = rgroup.getSingletons();
        while (singleton.hasMoreElements()) {
            Molecule m = singleton.nextElement();
            model.addRow(new Object[]{m.getName(), m});
        }
        singletonTab.setModel(model);
        navtab.setTitleAt(SINGLETON_TAB, "Singletons ("
                          +rgroup.getSingletonCount()+")");
        navtab.setTitleAt(SCAFFOLD_TAB, "Scaffolds ("
                          +rgroup.getScaffoldCount()+")");

        // now regenerate all datasets
        Vector files = new Vector ();
        for (int i = contentTab.getTabCount(); --i > 0;) {
            final JComponent c = (JComponent)contentTab.getComponentAt(i);
            SearchRGroup rg = (SearchRGroup)c.getClientProperty
                ("rgroup.search");
            if (rg != null) {
                files.add(rg.getFile());
                contentTab.remove(i);
            }

            if (i == COMPOUND_TAB) {
                contentTab.setTitleAt
                    (i, "Compounds ("+rgroup.size()+")");
            }
        }

        for (Object file : files) {
            new SearchRGroup (rgroup.getRGroupTables(), 
                              file).execute();
        }
    }

    public void valueChanged (ListSelectionEvent e) {
        if (e.getValueIsAdjusting()) {
            return;
        }

        // scaffold selection
        RGroupTable rtab = getSelectedScaffold ();
        if (rtab != null) 
            showRGroupTable (rtab);
    }

    void showRGroupTable (final RGroupTable rtab) {
        contentTab.setTitleAt(RGROUP_TAB, "R-group ("+rtab.getRowCount()+")");
        rgroupTab.setModel(rtab);
        //addMolecularSorter(rgroupTab);
        Enumeration<TableColumn> tcI=rgroupTab.getColumnModel().getColumns();
        for(TableColumn tc=null;tcI.hasMoreElements();){
            tc=tcI.nextElement();
            int mindex=tc.getModelIndex();
            if(rtab.isVisible(mindex)){
                if(tc.getWidth()>0){
                    columnPreviousWidth.put(mindex,tc.getWidth());
                    tc.setMinWidth(0);
                    tc.setPreferredWidth(0);
                    tc.setResizable(false);
                }
            }
            else {
                if(tc.getWidth()<=0){
                    Integer restoreWidth=columnPreviousWidth.get(mindex);
                    tc.setResizable(true);
                    tc.setMinWidth(5);
                    tc.setPreferredWidth(restoreWidth);
                }
            }
        }
        rgroupTab.getTableHeader().repaint();
        EventQueue.invokeLater(new Runnable () {
                public void run () {
                    // clear all previous highlighting
                    for (int i = 0; i < getRGroup().size(); ++i) {
                        Molecule m = getRGroup().get(i);
                        for (MolAtom a : m.getAtomArray()) {
                            a.setSetSeq(0);
                        }
                    }
                    // do alignment and highlighting
                    rtab.align();

                    // repaint the two tables
                    rgroupTab.repaint();
                    instanceTab.repaint();
                }
            });

        // select the rest of the table
        for (int i = 1; i < contentTab.getTabCount(); ++i) {
            JComponent c = (JComponent)contentTab.getComponentAt(i);
            JPanel pane = (JPanel)c.getClientProperty("rgroup.pane");
            final JTable tab = (JTable)c.getClientProperty("rgroup.table");
            //addMolecularSorter(tab);
            SearchRGroup rg = (SearchRGroup)
                c.getClientProperty("rgroup.search");
            if (rg != null) {
                final RGroupTable rt = rg.getRGroupTable(rtab);
                if (rt != null) {
                    contentTab.setTitleAt
                        (i, rg.getName()+" ("+rt.getRowCount()+")");
                    tab.setModel(rt);
                    ((CardLayout)pane.getLayout()).show(pane, "rgroup.data");
                    EventQueue.invokeLater(new Runnable () {
                            public void run () {
                                rt.align();
                                tab.repaint();
                            }
                        });
                }
                else {
                    contentTab.setTitleAt(i, rg.getName());
                    ((CardLayout)pane.getLayout()).show(pane, "rgroup.empty");
                }
            }
        }
        
        extBtn.setEnabled(rtab.isExtensible());

        for (int i = 0; i < popup.getComponentCount(); ++i) {
            Component c = popup.getComponent(i);
            if (c instanceof JMenuItem) {
                JMenuItem item = (JMenuItem)c;
                if (item.getText().equalsIgnoreCase("extend")) {
                    item.setEnabled(rtab.isExtensible());
                }
            }
        }

        Molecule scaffold = rtab.getScaffold();
        if (scaffold != null) {
            mview.setM(0, scaffold);
        }

        contentTab.setTitleAt(REFERENCE_TAB, "References...");
        docsPane.search(rtab.getCore());
        //contentTab.setTitleAt(PUBCHEM_TAB, "PubChem...");
        //pubchemPane.search(rtab.getCore());
    }

    void dumpXml (PrintStream ps, RGroupTable rtab) {
        ps.println("<?xml version=\"1.0\"?>");
        ps.println("<rgroup-radial>");
        
        ps.println("</rgroup-radial>");
    }

    protected StructureSearchDialog getStructureSearchDialog () {
        if (_strucSearchDialog == null) {
            _strucSearchDialog = new StructureSearchDialog (this) {
                    protected void search (Molecule query, 
                                           SearchParams params) {
                        contentTab.setTitleAt(REFERENCE_TAB, "References...");
                        docsPane.search(query, params);
                    }
                };
            Dimension s0 = getSize ();
            Dimension s1 = _strucSearchDialog.getSize();
            _strucSearchDialog.setLocation
                (getX() + (s0.width - s1.width)/2, 
                 getY() + (s0.height - s1.height)/2);
        }
        return _strucSearchDialog;
    }

    protected Component createMViewPane () {
        JPanel pane = new JPanel (new BorderLayout (0, 2));
        pane.add(new JLabel ("Double-click to edit scaffold", JLabel.CENTER), 
                 BorderLayout.NORTH);
        
        mview = new MViewPane ();
        /*
          mview.getVisibleCellComponent(0)
          .addMouseMotionListener(new MouseMotionAdapter () {
          public void mouseDragged (MouseEvent e) {
          System.out.println
          ("mouse dragged: " + e.getX() + " " + e.getY());
          }
          });
        */
        mview.getVisibleCellComponent(0)
            .addMouseListener(new MouseAdapter () {
                    public void mouseReleased (MouseEvent e) {
                        updateSelectedScaffold ();
                    }
                });
        mview.addPropertyChangeListener(this);
        mview.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    System.out.println(e.getActionCommand());
                }
            });
        mview.setEditable(2);
        pane.add(mview);

        return pane;
    }

    void updateSelectedScaffold () {
        Molecule m = mview.getM(0);
        if (m == null) return;

        mview.applyRotationMatrices();
        RGroupTable rtab = getSelectedScaffold ();

        if (rtab != null) {
            Molecule core = rtab.getCore();
            Molecule scaffold = rtab.getScaffold();
        }

        //System.out.println("scaffold: " + scaffold + " edited: " + m);
    }

        RGroupTable getSelectedScaffold() {

                int row = scaffoldTab.getSelectedRow();
                return row < 0 ? null
                                : ((ScaffoldTableModel) scaffoldTab.getModel()).getScaffold(scaffoldTab.convertRowIndexToModel(row));
        }

        protected JTable createTable() {
                JTable tab = new RTable();
                tab.setDefaultRenderer(RGroupTable.CRC.class,
                                new CRCCellRenderer());
                return tab;
        }

    Component createRGroupPane () {
        JPanel pane = new JPanel (new BorderLayout ());
        /*
          pane.setBorder(BorderFactory.createCompoundBorder
          (BorderFactory.createTitledBorder("R-group"), EMPTY));
        */
        pane.setBorder(EMPTY);
        rgroupTab = createTable ();
        pane.add(new JScrollPane (rgroupTab));
        return pane;
    }

    Component createInstancesPane () {
        JPanel pane = new JPanel (new BorderLayout ());
        /*
          pane.setBorder(BorderFactory.createCompoundBorder
          (BorderFactory.createTitledBorder("R-group"), EMPTY));
        */
        pane.setBorder(EMPTY);
        instanceTab = createTable ();
        instanceTab.getSelectionModel()
            .setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        instanceTab.getSelectionModel().addListSelectionListener
            (new ListSelectionListener () {
                    public void valueChanged (ListSelectionEvent e) {
                        // this is in view coordinate; if sorting is available
                        //  the row should be converted into model's 
                        //  coordinate!
                        ScaffoldTableModel model = 
                            (ScaffoldTableModel)scaffoldTab.getModel();
                        int row = instanceTab.getSelectedRow();
                        row= (row < 0)?-1:instanceTab.convertRowIndexToModel(row);                        
                        if (row >= 0) {
                            RGroupTable[] rgs = 
                                getRGroup().getRGroupTablesFor(row);
                            model.setScaffolds(rgs);
                        }
                        else {
                            model.setScaffolds((RGroupTable[])null); // clear
                        }
                    }
                });

        pane.add(new JScrollPane (instanceTab));
        //addMolecularSorter(instanceTab);
        return pane;
    }

    Component createDocsPane () {
        docsPane = new RGroupDocsPane ();
        docsPane.setWS(webservice);
        docsPane.addPropertyChangeListener(this);
        return docsPane;
    }

    Component createSingletonPane () {
        JPanel pane = new JPanel (new BorderLayout ());
        singletonTab = createTable ();
        pane.add(new JScrollPane (singletonTab));
        return pane;
    }

    Component createPubChemPane () {
        pubchemPane = new RGroupPubChemPane ();
        pubchemPane.addPropertyChangeListener(this);

        return pubchemPane;
    }

    void makeRGroupStats() throws IOException {
        // read in data
        BufferedReader reader = new BufferedReader(new FileReader(activityDataFileName));
        HashMap<String, Double> activity = new HashMap<String, Double>();
        String line;
        while ((line = reader.readLine()) != null) {
            String[] toks = line.trim().split(",");
            activity.put(toks[0], Double.parseDouble(toks[1]));
        }
        System.out.println("Loaded activity data from " + activityDataFileName);
        System.out.println("Will process " + getRGroup().getScaffoldCount() + " scaffolds");

        int nskip = 0;
        for (int i = 0; i < getRGroup().getScaffoldCount(); i++) {
            RGroupTable rtab = getRGroup().getRGroupTable(i);
            double ratio = (double) rtab.getRGroupCount() / rtab.getRowCount();
            if (ratio >= 1) {
                //                System.out.println("Won't fit model for scaffold " + (i + 1) + " (ratio = " + ratio + ")");
                nskip++;
                continue;
            }

            Molecule[] members = rtab.getMembers();
            int ndesc = 2;
            double[][] x = new double[rtab.getRowCount()][rtab.getRGroupCount() * ndesc];
            double[] y = new double[rtab.getRowCount()];
            String[] molids = new String[rtab.getRowCount()];

            for (int row = 0; row < rtab.getRowCount(); row++) {
                Molecule[] groups = rtab.getRGroups(row);

                // get molecule id for this row
                String moleculeId = members[row].getName();
                molids[row] = moleculeId;
                y[row] = activity.get(moleculeId);

                // generate descriptors
                for (int col = 0; col < rtab.getRGroupCount(); col += ndesc) {
                    Molecule m = groups[col];
                    if (i == 11) System.out.println("m = " + m);
                    if (m != null) {
                        TopologicalIndices ti = new TopologicalIndices(m);
                        x[row][col + 0] = ti.TIkierflex();
                        x[row][col + 1] = m.getAtomCount();
                    } else {
                        for (int j = 0; j < ndesc; j++)
                            x[row][col + j] = 0.0;
                    }
                }
                System.out.println();
            }

            // dump x, y for testing purposes
            dumpRGroupDesc("scaf" + i + ".csv", molids, y, x);
            // fit model

        }
        System.out.println("Processed = " 
                           + (getRGroup().getScaffoldCount() - nskip) 
                           + " scaffolds, skipped " + nskip);
    }


    void addSearchRGroup (SearchRGroup rgroup) {
        JPanel rgpane = new JPanel (new BorderLayout ());
        rgpane.setBorder(EMPTY);
        final JTable rgtab = createTable ();
        rgpane.add(new JScrollPane (rgtab));
        
        JPanel spane = new JPanel (new BorderLayout ());
        JTable stab = createTable ();
        stab.setModel(rgroup.getSingletons());
        spane.add(new JScrollPane (stab));

        spane.setBorder(BorderFactory.createCompoundBorder
                        (BorderFactory.createTitledBorder
                         ("Singletons" + (stab.getRowCount() == 0 ? ""
                                          : (" ("+stab.getRowCount()+")"))),
                         EMPTY));

        JPanel pane = new JPanel (new CardLayout ());
        pane.add(rgpane, "rgroup.data");

        JPanel empty = new JPanel (new BorderLayout ());
        empty.add(new JLabel ("<html><b>This scaffold is not in data set "
                              + rgroup.getName(), JLabel.CENTER));
        pane.add(empty, "rgroup.empty");

        JSplitPane split = new JSplitPane (JSplitPane.VERTICAL_SPLIT);
        split.setDividerSize(7);
        split.setOneTouchExpandable(true);
        split.setTopComponent(pane);
        split.setBottomComponent(spane);
        split.setResizeWeight(.75);
        split.putClientProperty("rgroup.table", rgtab);
        split.putClientProperty("rgroup.search", rgroup);
        split.putClientProperty("rgroup.pane", pane);

        contentTab.addTab(rgroup.getName(), split);

        RGroupTable rtab = getSelectedScaffold ();
        if (rtab != null) {
            final RGroupTable tab = rgroup.getRGroupTable(rtab);
            if (tab != null) {
                rgtab.setModel(tab);
                contentTab.setTitleAt
                    (contentTab.getTabCount()-1, rgroup.getName()
                     +" ("+tab.getRowCount()+")");

                EventQueue.invokeLater(new Runnable () {
                        public void run () {
                            tab.align();
                            rgtab.repaint();
                        }
                    });
            }
            else {
                ((CardLayout)pane.getLayout()).show(pane, "rgroup.empty");
            }
        }
    }

    protected JToolBar createToolBar () {
        JToolBar toolbar = new JToolBar ();
        JButton btn = new JButton (BACK_ICON);
        btn.setBorderPainted(false);
        btn.setContentAreaFilled(false);
        btn.setRolloverEnabled(true);
        btn.setToolTipText("Go to previous dataset");
        btn.setActionCommand("go-back");
        btn.addActionListener(this);
        btn.setEnabled(false);
        navBackBtn = btn;
        toolbar.add(btn);

        btn = new JButton (FORWARD_ICON);
        btn.setBorderPainted(false);
        btn.setContentAreaFilled(false);
        btn.setRolloverEnabled(true);
        btn.setToolTipText("Go to next dataset");
        btn.setActionCommand("go-forward");
        btn.addActionListener(this);
        btn.setEnabled(false);
        navForwardBtn = btn;
        toolbar.add(btn);
        
        btn = new JButton (OPEN_ICON);
        btn.setBorderPainted(false);
        btn.setContentAreaFilled(false);
        btn.setRolloverEnabled(true);
        btn.setToolTipText("Open dataset");
        btn.setActionCommand("open");
        btn.addActionListener(this);
        toolbar.add(btn);

        btn = new JButton (COMPOUND_ICON);
        btn.setBorderPainted(false);
        btn.setContentAreaFilled(false);
        btn.setRolloverEnabled(true);
        btn.setToolTipText("Structure search references");
        btn.setActionCommand("structure-search");
        btn.addActionListener(this);
        toolbar.add(btn);

        toolbar.addSeparator();
        toolbar.add(new JLabel ("Search: "));
        searchField = new JTextField (30);
        searchField.setToolTipText
            ("Enter search terms (e.g., clk4, pde4 inhibitor)");
        searchField.setActionCommand("search");
        searchField.addActionListener(this);
        toolbar.add(searchField);

        return toolbar;
    }

    protected JMenuBar createMenuBar () {
        JMenuBar menubar = new JMenuBar ();
        JMenu menu;
        JMenuItem item;
        menubar.add(menu = new JMenu ("File"));

        JMenu open = new JMenu ("Open");
        menu.add(open);
        
        open.add(item = new JMenuItem ("File"));
        item.setToolTipText("Open input file to generate scaffolds");
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    openFile ();
                }
            });
        open.add(item = new JMenuItem ("URL"));
        item.setToolTipText("Open input stream from a URL");
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    openURL ();
                }
            });
        open.add(item = new JMenuItem ("DOI/PubMed"));
        item.setToolTipText("Open input from DOI/PubMed references");
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    openWS ();
                }
            });

                menu.addSeparator();
                JMenu exportDataMenu = new JMenu("Export dataset");
                JMenuItem exportJsonMenuItem = new JMenuItem("JSON");
                JMenuItem exportTabMenuItem = new JMenuItem("TSV");
                exportDataMenu.add(exportJsonMenuItem);
                exportDataMenu.add(exportTabMenuItem);
                menu.add(exportDataMenu);
                exportJsonMenuItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                exportDatasetAsJSON();
                        }
                });
                exportTabMenuItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                exportDatasetAsTSV();
                        }
                });

        menu.add(item = new JMenuItem ("Export image"));
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    exportImage ();
                }
            });
        menu.addSeparator();
        
        menu.add(item = new JMenuItem ("Quit"));
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    quit ();
                }
            });

        menubar.add(menu = new JMenu ("Options"));
        menu.add(item = new JMenuItem ("Add/remove data columns"));
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    propDialog.setVisible(true);
                }
            });

        menubar.add(menu = new JMenu ("Help"));
        menu.add(item = new JMenuItem ("About"));
        item.addActionListener(new ActionListener () {
                public void actionPerformed (ActionEvent e) {
                    about ();
                }
            });
        
        return menubar;
    }

    void about () {
        JOptionPane.showMessageDialog
            (this, "NCGC Scaffold Hopper, 2017\n"
             +"Please send questions and/or comments to\n"
             +"nguyenda@mail.nih.gov and tyler.peryea@nih.gov",
             "About", JOptionPane.INFORMATION_MESSAGE);
    }

    void newRGroupGenerator (RGroupWorker worker) {
        newRGroupGenerator (false, worker);
    }

    void newRGroupGenerator (boolean clearHistory, RGroupWorker worker) {
        if (clearHistory) {
            for (RGroupGenerator rg : history) {
                rg.removePropertyChangeListener(this);
            }
            history.clear();

            for (RGroupGenerator rg : forward) {
                rg.removePropertyChangeListener(this);
            }
            forward.clear();
            navForwardBtn.setEnabled(false);
        }

        RGroupGenerator rgroup = new RGroupGenerator ();
        rgroup.addPropertyChangeListener(this);
        history.push(rgroup);
        navBackBtn.setEnabled(history.size() > 1);
        progress.setStringPainted(true);
        worker.execute();
    }

    RGroupGenerator getRGroup () { 
        return history.peek();
    }

    void loadFile (String... argv) {
        if (argv.length > 0) {
            /*
              rgroup = new RGroupGenerator ();
              rgroup.addPropertyChangeListener(this);
              progress.setStringPainted(true);
              new RGroupWorker (argv).execute();
            */
            newRGroupGenerator (true, new RGroupWorker (argv));
        }
    }

    void loadFile (File... files) {
        if (files == null || files.length == 0) {
            JOptionPane.showMessageDialog
                (this, "No file(s) selected!", "Error", 
                 JOptionPane.ERROR_MESSAGE);
        }
        else {
            /*
              rgroup = new RGroupGenerator ();
              rgroup.addPropertyChangeListener(this);
              progress.setStringPainted(true);
              new RGroupWorker (files).execute();
            */
            newRGroupGenerator (true, new RGroupWorker (files));
        }
    }

    void loadFile (FileContents... files) {
        if (files == null || files.length == 0) {
            JOptionPane.showMessageDialog
                (this, "No file(s) selected!", "Error", 
                 JOptionPane.ERROR_MESSAGE);
        }
        else {
            /*
              rgroup = new RGroupGenerator ();
              rgroup.addPropertyChangeListener(this);
              progress.setStringPainted(true);
              new RGroupWorker (files).execute();
            */
            newRGroupGenerator (true, new RGroupWorker (files));
        }       
    }

    private void dumpRGroupDesc(String fname, String[] molids, double[] y, double[][] x) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(fname));
        for (int j = 0; j < y.length; j++) {
            StringBuffer sb = new StringBuffer();
            if (molids != null) sb.append(molids[j] + ",");
            sb.append(y[j] + ",");
            for (int k = 0; k < x[0].length; k++) {
                if (k != x[0].length - 1) sb.append(x[j][k] + ",");
                else sb.append(x[j][k]);
            }
            sb.append("\n");
            writer.write(sb.toString());
        }
        writer.close();
    }

    public void propertyChange (PropertyChangeEvent e) {
        String name = e.getPropertyName();
        //logger.info(e.getSource() + ": name="+name+" old="+e.getOldValue() + " new="+e.getNewValue());

        if (name.equals("progress")) {
            progress.setValue((Integer)e.getNewValue());
        }
        else if (name.equals("docs")) {
            JsonNode[] docs = (JsonNode[])e.getNewValue();
            for (int i = 0; i < contentTab.getTabCount(); ++i) {
                if (i == REFERENCE_TAB) {
                    contentTab.setTitleAt
                        (i, "References"+(docs != null 
                                          ? (" ("+docs.length+")") : ""));
                }
            }
        }
        else if (name.equals("pubchem")) {
            contentTab.setTitleAt
                (PUBCHEM_TAB, "PubChem ("+pubchemPane.getCount()+")");
        }
        else if (name.equals("load")) {
            // we should keep a history of this
            URL url = getWS().getCompoundURL(e.getNewValue()+"/pubmed");
            newRGroupGenerator (false, new RGroupWorker (url));
        }
        else if (name.equalsIgnoreCase("majr") 
                 || name.equalsIgnoreCase("mh")) {
            String mesh = "\""+e.getNewValue()+"\"["+name+"]";
            searchField.setText(mesh);
            doSearch (mesh);
        }
        else if (name.equalsIgnoreCase("tid")) {
            URL url = getWS().getDocumentURL(e.getNewValue()+"/tid?max=100");
            doLoad (url);
        }
        else if (name.equals("model")) {
            ((AbstractTableModel)scaffoldTab.getModel()).addTableModelListener
                (new TableModelListener () {
                        public void tableChanged (TableModelEvent e) {
                            navtab.setTitleAt
                                (SCAFFOLD_TAB, "Scaffolds ("
                                 +((TableModel)e.getSource()).getRowCount()
                                 +")");
                        }
                    });
        }
        else if (name.equals("rgroup")) {
            if (e.getSource() == network) {
                RGroupTable rtab = (RGroupTable)e.getNewValue();
                if (rtab != null) {
                    selectScaffold (rtab);
                }
            }
        }
        else if(name.equals("status")){
            statusField.setText((String)e.getNewValue());
        }

        if (mview != e.getSource()) {
            return;
        }

        /*
          if (name.equals("selectedIndex")) {
          mview.applyRotationMatrices();
          Molecule m = mview.getM(0).cloneMolecule();
          // remove the R-group...
          Set<MolAtom> remove = new HashSet<MolAtom>();
          for (MolAtom a : m.getAtomArray()) {
          if (a.getAtno() == MolAtom.RGROUP) {
          remove.add(a);
          }
          }
          for (MolAtom a : remove) {
          m.removeNode(a);
          }
          mview.setM(0, m);
          }
        */
    }

    void selectScaffold (RGroupTable rtab) {
        /*
         * this code doesn't work when the scaffoldTab is reduced
         * from the Instance selection. So to make life easy, we 
         * just clear the instanceTab selection.
         */
        instanceTab.clearSelection();
        RGroupGenerator rgroup = getRGroup ();
        for (int r = 0; r < rgroup.getScaffoldCount(); ++r) {
            if (rgroup.getRGroupTable(r) == rtab) {
                if(r<scaffoldTab.getRowCount()&&r>=0){
                    int vr = scaffoldTab.convertRowIndexToView(r);
                    scaffoldTab.setRowSelectionInterval(vr, vr);
                    int scroll = vr*scaffoldTab.getRowHeight();
                    JScrollPane sp = (JScrollPane)
                        SwingUtilities.getAncestorOfClass
                        (JScrollPane.class, scaffoldTab);
                    if (sp != null) {
                        sp.getVerticalScrollBar().setValue(scroll);
                    }
                    logger.info(vr+"th scaffold selected!");
                    break;
                }
            }
        }
    }

    void openURL () {
        String value = JOptionPane.showInputDialog
            (this, "Please enter URL", null);
        if (value != null) {
            try {
                URL url = new URL (value);
                /*
                  rgroup = new RGroupGenerator ();
                  rgroup.addPropertyChangeListener(this);
                  progress.setStringPainted(true);
                  new RGroupWorker (url).execute();
                */
                newRGroupGenerator (false, new RGroupWorker (url));
            }
            catch (Exception ex) {
                JOptionPane.showMessageDialog
                    (this, "Bad URL "+value, "Error", 
                     JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    void openWS () {
        String value = JOptionPane.showInputDialog
            (this, "Please enter DOIs and/or PubMed IDs", null);

        if (value != null) {
            String[] toks = value.split("[\\s]+");
            ArrayList<URL> urls = new ArrayList<URL>();
            for (String t : toks) {
                if (t.startsWith("10.")) {
                    urls.add(getWS().getCompoundURL(t+"/doi"));
                }
                else {
                    urls.add(getWS().getCompoundURL(t+"/pubmed"));
                }                
            }
        
            if (urls.isEmpty()) {
                JOptionPane.showMessageDialog
                    (this, "No valid DOI/PubMed IDs specified!", "Error", 
                     JOptionPane.ERROR_MESSAGE);
            }
            else {
                newRGroupGenerator (true, new RGroupWorker
                                    (urls.toArray(new URL[0])));
            }
        }
    }

    void doSearch (String query) {
        logger.info("## Searching for \""+query+"\"...");
        contentTab.setTitleAt(REFERENCE_TAB, "References...");
        docsPane.search(query);
    }

    void doLoad (URL url) {
        logger.info("## Loading "+url+"...");
        contentTab.setTitleAt(REFERENCE_TAB, "References...");
        docsPane.load(url);        
    }

    void goBack () {
        forward.push(history.pop());
        docsPane.push();
        navBackBtn.setEnabled(history.size() > 1);
        navForwardBtn.setEnabled(true);
        updateRGroup ();
    }

    void goForward () {
        history.push(forward.pop());
        docsPane.pop();
        navBackBtn.setEnabled(true);
        navForwardBtn.setEnabled(!forward.isEmpty());
        updateRGroup ();
    }

    void openFile () {
        if (UnsignedWebstart) {
            try {
                FileOpenService fos = (FileOpenService)ServiceManager.lookup
                    ("javax.jnlp.FileOpenService");
                FileContents[] files = fos.openMultiFileDialog
                    (".", new String[]{"sdf","mol","smi","smiles",
                                       "cml","sd"});
                loadFile (files);
            }
            catch (Exception ex) {
                JOptionPane.showMessageDialog
                    (RGroupTool.this, "Your Java Webstart doesn't have "
                     +"proper permission to open files.", 
                     "Fatal Error", JOptionPane.ERROR_MESSAGE);
            }
        }
        else {
            chooser.setDialogTitle("Please select input file(s)");
            if (JFileChooser.APPROVE_OPTION == chooser.showOpenDialog(this)) {
                loadFile (chooser.getSelectedFiles());
            }
        }
    }

    RGroupWebResources getWS () { return webservice; }

    void quit () { 
        System.exit(0);
    }


    private JsonNode _addMembersToScaffold(RGroupTable rgt, ObjectNode scaffoldNode, ObjectMapper mapper) {
        ArrayNode mols = mapper.createArrayNode();

        // add the rgroup labels to the scaffold node
        ArrayNode labels = mapper.createArrayNode();
        for (int i = 0; i < rgt.getRGroupCount(); i++) labels.add(rgt.getRGroupLabel(i));
        scaffoldNode.put("rgroupLabels", labels);

        // molecules for this r-group table
        for (int idx = 0; idx < rgt.rows.length; idx++) {
            int serial = rgt.rows[idx];
            ObjectNode molnode = mapper.createObjectNode();

            Molecule m = rgt.molseq.get(serial);
            molnode.put("id", m.getName());
            molnode.put("smiles", m.toFormat("cxsmiles"));

            // r-groups for this molecule
            ArrayNode rgroupNodes = mapper.createArrayNode();
            Molecule[] rgroups = rgt.rgroups[idx];
            for (int rpos = 0; rpos < rgroups.length; rpos++) {
                String label = rgt.getRGroupLabel(rpos);
                Molecule rgm = rgroups[rpos];
                if (rgm != null) {
                    ObjectNode tmp = mapper.createObjectNode();
                    tmp.put(label, rgm.toFormat("cxsmiles"));
                    rgroupNodes.add(tmp);
                }
            }
            molnode.put("rgroups", rgroupNodes);
            mols.add(molnode);
        }
        return mols;
    }

    void exportDatasetAsTSV() {
        RGroupGenerator rgen = getRGroup();
        if (scaffoldTab.getModel() == null || scaffoldTab.getModel() instanceof DefaultTableModel) return;
        ScaffoldTableModel model = (ScaffoldTableModel) scaffoldTab.getModel();

        StringBuilder sb = new StringBuilder();
        StringBuilder rsb = new StringBuilder();
        int RGROUP_MAX = 21;

        RGroupTable selectedScaffold = getSelectedScaffold();
        if (selectedScaffold != null) {

        } else {

            String delim = "";

            // set up header
            sb.append(delim).append("ScaffoldId");
            delim = "\t";
            sb.append(delim).append("Structure");
            // add in R-group labels
            sb.append(delim).append("RgroupLabels");

            for (int col = 1; col < model.getColumnCount(); col++) {
                Object data = model.getValueAt(0, col);
                String colName = model.getColumnName(col);
                if (data instanceof ScaffoldData) {
                    sb.append(delim).append(colName).append(delim).append(colName+"_SD");
                } else sb.append(delim).append(colName);
            }
            sb.append("\n");

            // rgroup file header
            rsb.append("ScaffoldID").
                append(delim).
                append("MolID").
                append(delim).
                append("Structure");
            for (int i = 1; i <= RGROUP_MAX; i++) rsb.append(delim).append("R"+i);
            rsb.append("\n");

            delim = "";
            for (int row = 0; row < model.getRowCount(); row++) {
                sb.append(delim);
                delim = "\t";

                Molecule scaffold = (Molecule) model.getValueAt(row, 0);
                sb.append(delim).append(row);
                sb.append(delim).append(scaffold.toFormat("cxsmiles"));

                // get set of R? labels for this scaffold
                RGroupTable rgt = rgen.getRGroupTable(row);
                /*
                String rlabels = IntStream.
                    range(0, rgt.getRGroupCount()).
                    mapToObj(i -> rgt.getRGroupLabel(i)).
                    collect(Collectors.joining(","));
                sb.append(delim).append(rlabels);
                */
                sb.append(delim);
                for (int i = 0; i < rgt.getRGroupCount(); ++i) {
                    if (i > 0) sb.append(",");
                    sb.append(rgt.getRGroupLabel(i));
                }

                // get property values
                for (int col = 1; col < model.getColumnCount(); col++) {
                    Object data = model.getValueAt(row, col);
                    if (data instanceof ScaffoldData) {
                        ScaffoldData sd = (ScaffoldData) data;
                        sb.append(delim).append(sd.getMean()).append(delim).append(sd.getStd());
                    } else {
                        sb.append(delim).append(data);
                    }
                }
                sb.append("\n");

                // get r-group table for this scaffold. It will be written to a separate file
                String rdelim = "";
                for (int idx = 0; idx < rgt.rows.length; idx++) {
                    int serial = rgt.rows[idx];
                    Molecule m = rgt.molseq.get(serial);

                    rsb.append(rdelim).append(row);
                    rdelim = "\t";
                    rsb.append(rdelim).append(m.getName()).append(delim).append(m.toFormat("cxsmiles"));

                    // r-groups for this molecule
                    int nrg = rgt.getRGroupCount();
                    Molecule[] rgroups = rgt.rgroups[idx];
                    for (int rpos = 0; rpos < RGROUP_MAX; rpos++) {
                        if (rpos >= nrg)
                            rsb.append(delim).append("");
                        else {
                            Molecule rgm = rgroups[rpos];
                            if (rgm != null) rsb.append(delim).append(rgm.toFormat("cxsmiles"));
                            else rsb.append(delim).append("");
                        }
                    }
                    rsb.append("\n");
                    rdelim = "";
                }
            }
        }

        chooser.setDialogTitle("Please select output TSV file");
        if (JFileChooser.APPROVE_OPTION == chooser.showSaveDialog(this)) {
            try {
                File file = chooser.getSelectedFile();
                FileWriter w = new FileWriter(file);
                w.write(sb.toString());
                w.close();

                String fname = file.getParent()+ File.separator+"rgoup-"+file.getName();
                w = new FileWriter(new File(fname));
                w.write(rsb.toString());
                w.close();
            } catch (Exception ex) {
                ex.printStackTrace();
                JOptionPane.showMessageDialog
                    (this, "Unable to export dataset as TSV", "Error",
                     JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    void exportDatasetAsJSON() {
        ObjectMapper mapper = new ObjectMapper();
        ArrayNode scaffoldList = mapper.createArrayNode();

        RGroupGenerator rgen = getRGroup();
        if (scaffoldTab.getModel() == null || scaffoldTab.getModel() instanceof DefaultTableModel) return;
        ScaffoldTableModel model = (ScaffoldTableModel) scaffoldTab.getModel();

        RGroupTable selectedScaffold = getSelectedScaffold();

        if (selectedScaffold != null) {
            ObjectNode scaffoldNode = mapper.createObjectNode();
            Molecule scaffold = selectedScaffold.getScaffold();
            scaffoldNode.put("scaffold",scaffold.toFormat("cxsmiles"));
            scaffoldNode.put("Scaffold Score", selectedScaffold.getScaffoldScore());
            scaffoldNode.put("Complexity", selectedScaffold.getScaffoldComplexity());
            scaffoldNode.put("Count", selectedScaffold.getMemberCount());
            scaffoldNode.put("members", _addMembersToScaffold(selectedScaffold, scaffoldNode, mapper));
            scaffoldList.add(scaffoldNode);
        } else {
            for (int row = 0; row < model.getRowCount(); row++) {

                ObjectNode scaffoldNode = mapper.createObjectNode();
                Molecule scaffold = (Molecule) model.getValueAt(row, 0);
                scaffoldNode.put("scaffold", scaffold.toFormat("cxsmiles"));

                for (int col = 1; col < model.getColumnCount(); col++) {
                    Object data = model.getValueAt(row, col);
                    if (data instanceof ScaffoldData) {
                        ScaffoldData sd = (ScaffoldData) data;
                        ObjectNode sdNode = mapper.createObjectNode();
                        sdNode.put("mean", sd.getMean());
                        sdNode.put("sd", sd.getStd());
                        scaffoldNode.put(model.getColumnName(col), sdNode);
                    } else {
                        Class colClass = model.getColumnClass(col);
                        if (colClass.getName().contains("Integer"))
                            scaffoldNode.put(model.getColumnName(col), (Integer) data);
                        else if (colClass.getName().contains("Double"))
                            scaffoldNode.put(model.getColumnName(col), (Double) data);
                        else if (colClass.getName().contains("Float"))
                            scaffoldNode.put(model.getColumnName(col), (Float) data);
                        else if (colClass.getName().contains("Long"))
                            scaffoldNode.put(model.getColumnName(col), (Long) data);
                        else if (colClass.getName().contains("String"))
                            scaffoldNode.put(model.getColumnName(col), (String) data);
                    }
                }

                // Pull in r-group table for this scaffold
                RGroupTable rgt = rgen.getRGroupTable(row);
                JsonNode mols = _addMembersToScaffold(rgt, scaffoldNode, mapper);
                scaffoldNode.put("members", mols);
                scaffoldList.add(scaffoldNode);
            }
        }

        chooser.setDialogTitle("Please select output JSON file");
        if (JFileChooser.APPROVE_OPTION == chooser.showSaveDialog(this)) {
            try {
                _dumpJson(chooser.getSelectedFile(), scaffoldList, mapper);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog
                    (this, "Unable to export dataset as JSON!", "Error",
                     JOptionPane.ERROR_MESSAGE);
            }
        }

        //              try {
        //                      String json = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(scaffoldList);
        //                      System.out.println(json);
        //              } catch (IOException e) {
        //                      e.printStackTrace();
        //              }
    }

        private void _dumpJson(File out, JsonNode doc, ObjectMapper mapper) throws IOException {
                mapper.writeValue(out, doc);
        }

    void exportImage () {
        if (UnsignedWebstart) {
            try {
                FileSaveService fss = (FileSaveService)ServiceManager.lookup
                    ("javax.jnlp.FileSaveService");

                PipedInputStream pis = new PipedInputStream ();
                final PipedOutputStream pos = new PipedOutputStream (pis);
                Executors.newSingleThreadExecutor().execute(new Runnable () {
                        public void run () {
                            try {
                                exportViewImage (pos);
                                pos.close();
                            }
                            catch (Exception ex) {
                                ex.printStackTrace();
                            }
                        }
                    });

                FileContents save = fss.saveFileDialog
                    (null, null, pis, "scaffold.png");
                if (save != null) {
                    System.err.println("saved image: "+save.getName());
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
                JOptionPane.showMessageDialog
                    (RGroupTool.this, "Your Java Webstart doesn't have "
                     +"proper permission to write files.", 
                     "Fatal Error", JOptionPane.ERROR_MESSAGE);         
            }
        }
        else {
            chooser.setDialogTitle("Please select output image file");
            if (JFileChooser.APPROVE_OPTION == chooser.showSaveDialog(this)) {
                try {
                    exportViewImage (chooser.getSelectedFile());
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog
                        (this, "Unable to export image!", "Error", 
                         JOptionPane.ERROR_MESSAGE);
                }
            }
        }
    }

    void exportViewImage (File out) throws Exception {
        exportViewImage (new FileOutputStream (out));
    }

    void exportViewImage (OutputStream os) throws Exception {
        int row = scaffoldTab.getSelectedRow();
        row= (row < 0)?-1:scaffoldTab.convertRowIndexToModel(row);   
        if (row < 0) {
            JOptionPane.showMessageDialog
                (this, "No scaffold selected; no image generated!", 
                 "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        int htop = mview.getHeight() + 5;

        int index = contentTab.getSelectedIndex();
        JTable tab = rgroupTab; // the main r-group table
        if (index > 0) {
            JComponent c = (JComponent)contentTab.getTabComponentAt(index);
            SearchRGroup rg = (SearchRGroup)
                c.getClientProperty("rgroup.search");
            if (rg != null) {
                RGroupTable rtab = rg.getRGroupTable(getSelectedScaffold ());
                if (rtab == null) {
                    JOptionPane.showMessageDialog
                        (this, "There is no such scaffold exists in "
                         + rg.getName()+"!", "Error", 
                         JOptionPane.ERROR_MESSAGE);
                    return;
                }
                tab = (JTable)c.getClientProperty("rgroup.table");
            }
            else {
                JOptionPane.showMessageDialog
                    (this, "No scaffold selected; no image generated!", 
                     "Error", JOptionPane.ERROR_MESSAGE);
                return;
            }
        }

        int width = Math.max(mview.getWidth(), tab.getWidth());
        int height = htop + tab.getHeight() 
            + tab.getTableHeader().getHeight();
        BufferedImage img = new BufferedImage 
            (width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = img.createGraphics();
        g.setPaint(Color.white);
        g.fillRect(0, 0, width, height);
        int x = (width - mview.getWidth())/2;
        g.translate(x, 0);
        mview.paint(g);
        g.translate(-x, htop);
        tab.getTableHeader().paint(g);
        g.translate(0, tab.getTableHeader().getHeight());
        tab.paint(g);
        g.dispose();
        ImageIO.write(img, "png", os);
    }

    public static void main (final String[] argv) throws Exception {
        try {
            UIManager.setLookAndFeel(new Plastic3DLookAndFeel());
        } catch (Exception e) {}
        
        EventQueue.invokeLater(new Runnable () {
                public void run () {
                    RGroupWebResources web = null;
                    if (argv.length == 2) {
                        logger.info("Web Resource: host="
                                    +argv[0]+" context="+argv[1]);
                        web = new RGroupWebResources (argv[0], argv[1]);
                    }

                    RGroupTool tool = new RGroupTool (web);
                    tool.setSize(800, 600);
                    tool.setVisible(true);
                    if (web == null && argv.length > 0) {
                        tool.loadFile(argv);
                    }
                }
            });
    }

    public static class JnlpLaunch {
        public static void main (final String[] argv) throws Exception {
            UnsignedWebstart = true;
            try {
                UIManager.setLookAndFeel(new Plastic3DLookAndFeel());
            } catch (Exception e) {}
            
            EventQueue.invokeLater(new Runnable () {
                    public void run () {
                        RGroupTool tool = new RGroupTool ();
                        tool.setSize(800, 600);
                        tool.setVisible(true);
                    }
                });
        }
    }
}
