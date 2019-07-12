// $Id: RGroupTable.java 3864 2009-12-18 22:09:45Z nguyenda $
package gov.nih.ncgc.rgroup;


import gov.nih.ncgc.model.DataSeq;
import gov.nih.ncgc.util.MolAligner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.table.AbstractTableModel;

import chemaxon.marvin.util.MolExportException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;


public class RGroupTable extends AbstractTableModel{
    private static final String APP_PREFIX = "!$RGROUP$";
    private static final String CRC_SUFFIX = "-CRC";
    public static final String PREFIX_RGROUP= "!$RGROUP_R$";
    public static final String PREFIX_SCAFFOLD= "!$RGROUP_SCAFF$";
    public static final String SINGLETON_VALUE="$ST$";
	
    private static Logger logger = Logger.getLogger
	(RGroupTable.class.getName());

    protected DataSeq<Molecule> molseq; // input 
    protected Molecule scaffold; // decorated scaffold
    protected Molecule core; // original core (undecorated)
    protected int complexity; // scaffold complexity
    protected int[] fp; // scaffold fingerprint
    protected int[] rows = {}; // rows.length == rgroups.length
    protected int[][] hits = {}; // atom mapping against the scaffold
    
    protected int rgroupCount;
    protected Molecule[][] rgroups;
    protected String[] rgroupLabels = {}; // label for each group
    
    protected int xgroupCount=0;
    protected Molecule[][] xgroups;
    protected String[] xgroupLabels = {}; // label for each group
    
    
    protected Map[] ligands; // unique ligands for each r-group
   
    // boolean flag indicating whether it might be possible to extend
    //   further at the correpsonding r-group
    protected BitSet extension = new BitSet (); 
    protected double score; // R-Group score

    protected ArrayList<String> columns = new ArrayList<String>();
    protected Map<String, Class> props = new TreeMap<String, Class>();
    protected Set<String> generatedProps = new LinkedHashSet<String>();
    
    
    /*
     * This is the number of R-group columns to pad up to. This is only done to
     * preserve changes to a JTable's column ordering/resizing when there are a
     * different number of rgroups
     */
    private static final int PADDING_COLUMNS_RGROUP = 25;
    private static final int PADDING_COLUMNS_XGROUP = 3;
	
    
    Map<Integer, CRC[]> crcCache = new HashMap<Integer, CRC[]>();
    Set<Integer> invisibleColumns = new HashSet<Integer>();

    public static class CRC {
	static final double ln10 = 2.30258509299404568401;

	int size;
	double[][] xy;
	// four-parameter hill equation
	double zero, inf;
	double logac50, slope;

	CRC (int size) {
	    this.size = size;
	    xy = new double[2][size];
	}

	public int size () { return size; }
	public double[] x () { return xy[0]; }
	public double[] y () { return xy[1]; }
	public double[][] xy () { return xy; }

	public double[][] hillfit () {
	    if (zero != 0. || inf != 0. || logac50 != 0. || slope != 0.) {
		double[] x = upsample (upsample (upsample (upsample (xy[0]))));
		double[][] curve = new double[2][x.length];
		for (int i = 0; i < x.length; ++i) {
		    curve[0][i] = x[i];
		    curve[1][i] = hillf (x[i]);
		}

		return curve;
	    }
	    return null;
	}

	public double hillf (double x) {
	    return zero + ((inf - zero)
			   / (1. + Math.exp(ln10 * slope * 
					    (logac50 - x))));
	}

	private static double[] upsample (double[] x) {
	    if (x.length == 0) {
		return new double[]{};
	    }
	    
	    double[] x2 = new double[2*x.length-1];
	    for (int i=0; i<x2.length; i++) {
		if (i%2 == 0)
		    x2[i] = x[i/2];
		else
		    x2[i] = 0.5*(x[(i-1)/2]+x[(i+1)/2]);
	    }
	    return x2;
	}
    }

    protected RGroupTable (DataSeq<Molecule> molseq) {
	this.molseq = molseq;
	
    }

    public int getRGroupCount () { return rgroupCount; }
    public String getRGroupLabel (int g) { return rgroupLabels[g]; }
    public Map getRGroupLigands (int g) { return ligands[g]; }
    public Molecule getScaffold () { return scaffold; }
    public Molecule getCore () { return core; }
    public Molecule[] getRGroups (int row) { return rgroups[row]; }
    public Molecule getRGroup (int row, int rindex) {
	return rgroups[row][rindex];
    }
    public int getXGroupCount() {return xgroupCount;}
    public Molecule[] getXGroups (int row) { return xgroups[row]; }
    public Molecule getXGroup (int row, int rindex) {
	return xgroups[row][rindex];
    }
    public String getXGroupLabel(int j) {return xgroupLabels[j];}

    public boolean isExtensible () { return !extension.isEmpty(); }

    public int getScaffoldComplexity () { return complexity; }
    public double getScaffoldScore () { return score; }
    public int[] getScaffoldFp () { return fp; }

    public void align () {
    	//if(true)return;
	//MolHandler mh = new MolHandler (core);
	for (int i = 0; i < rows.length; ++i) {
		//System.out.println(i+":" +rows[i]);
	    Molecule m = molseq.get(rows[i]);
	    if (m.getDim() < 2) {
		m.clean(2, null);
	    }
	    //mh.align(m, hits[i]);
	    for (MolAtom a : m.getAtomArray()) {
                a.setSetSeq(0); // clear any previous coloring
                a.setAtomMap(0);
	    }
	    int[] h = hits[i];
	    
	    for (int j = 0; j < h.length; ++j) {
		MolAtom a = m.getAtom(h[j]);
		a.setSetSeq(3); // 1 - read, 2 - green, 3 - blue
		//a.setAtomMap(h[j]+1);
	    }
	    MolAligner.align(core, m, h);
	}
    }

    synchronized public void setupColumns(){
    	setupColumns(null);
    }

    /*
     * Setting up the default columns now requires a list of names to avoid, to
     * prevent collisions
     */
    synchronized protected void setupColumns (Set<String> avoidNames) {
    	columns.clear();
	if (columns.isEmpty()) {
            setColumnClass(APP_PREFIX+ "ID",String.class,avoidNames);
            setColumnClass(APP_PREFIX+"Structure",Molecule.class,avoidNames);
            for (int i = 0; i < xgroupCount; ++i) {
	    	setColumnClass(APP_PREFIX+ xgroupLabels[i],Molecule.class,avoidNames);
	    }
            for (int i=xgroupCount;i<PADDING_COLUMNS_XGROUP;i++){
                setColumnClass(APP_PREFIX+"_X" + i,Void.class,avoidNames);
                invisibleColumns.add(columns.size()-1);
            }
            for (int i = 0; i < rgroupCount; ++i) {
	    	setColumnClass(APP_PREFIX+ rgroupLabels[i],Molecule.class,avoidNames);
	    }		
            for (int i=rgroupCount;i<PADDING_COLUMNS_RGROUP;i++){
                setColumnClass(APP_PREFIX+"_R" + i,Void.class,avoidNames);
                invisibleColumns.add(columns.size()-1);
            }
            
	    fireTableStructureChanged ();
	}
    }
    private boolean setColumnClass(String lbl,Class clazz,Set<String> avoidNames){
    	String label=lbl+"";
    	boolean changed=false;
    	if(avoidNames!=null){
            while(avoidNames.contains(label)){
                label+="*";
                changed=true;
            }
    	}
        props.put(label, clazz);
        columns.add(label);
        return changed;
    }

    public void addProperties (Map<String, Class> props) {
	setupColumns (props.keySet());
	for (Map.Entry<String, Class> e : props.entrySet()) {
	    String name = e.getKey();
	    Class clazz = e.getValue();
	    if (name.endsWith(CRC_SUFFIX)) {
		clazz = CRC.class;
	    }
	    Class val = this.props.put(name, clazz);
	    if (val == null) {
		// new value
		columns.add(name);
	    }
	}
	fireTableStructureChanged ();
    }

    public void removeProperties (Collection<String> props) {
	for (String s : props) {
	    Class c = this.props.remove(s);
	    if (c != null) {
		columns.remove(s);
	    }
	}
	fireTableStructureChanged ();
    }

    public int getMemberCount () { return getRowCount (); }
    public int[] getMemberHits (int row) { return hits[row]; }
    public Molecule[] getMembers () {
	Vector<Molecule> mb = new Vector<Molecule>();
	for (int i = 0; i < rows.length; ++i) {
	    mb.add(molseq.get(rows[i]));
	}
	return mb.toArray(new Molecule[0]);
    }
    public BitSet getMemberBitset(){
    	BitSet bs = new BitSet();
    	for (int i = 0; i < rows.length; ++i) {
    	    bs.set(rows[i]);
    	}
    	return bs;
    }

    /**
     * TableModel interface
     */
    public int getRowCount () { 
	return rows.length; 
    }
    
    public int getColumnCount () { 
	return columns.size(); 
    }

    public Class getColumnClass (int col) {
    	if(props.get(columns.get(col)) == Molecule.class){
    		return Molecule.class;
    	}
    	return props.get(columns.get(col));
    }
    public String getColumnName (int col) {
        String cname = columns.get(col);
        if (cname.startsWith(this.APP_PREFIX)) {
            return cname.substring(this.APP_PREFIX.length());
        }
        return cname;
    }
	
    public Object getValueAt (int row, int col) {
	if (col == 0) {
	    return molseq.get(rows[row]).getName();
	}
	if (col == 1) {
	    return molseq.get(rows[row]);
	}
	if (col < xgroupCount+2){
		return xgroups[row][col-2];
	}
	if (col <PADDING_COLUMNS_XGROUP+2){
		return null;
	}
	if (col < rgroupCount+2+PADDING_COLUMNS_XGROUP) {
	    return rgroups[row][col-2-PADDING_COLUMNS_XGROUP];
	}

	Molecule m = molseq.get(rows[row]);
	if (m.getPropertyCount() == 0) {
	    return null;
	}

	String prop = columns.get(col);
	Object val = m.getPropertyObject(prop);
	if (val != null) {
	    Class clazz = props.get(prop);
	    if (clazz == String.class) {
		val = val.toString();
	    }
	    else if (clazz == Double.class) {
		try {
		    val = Double.parseDouble(val.toString());
		}
		catch (NumberFormatException ex) {
		}
	    }
	    else if (clazz == Long.class) {
		try {
		    val = Long.parseLong(val.toString());
		}
		catch (NumberFormatException ex) {
		}
	    }
	    else if (clazz == CRC.class) {
		CRC[] crc = crcCache.get(val.hashCode());
		if (crc == null) {
		    try {
			crc = parseCRC (val.toString());
			crcCache.put(val.hashCode(), crc);
		    }
		    catch (Exception ex) {
			logger.log(Level.WARNING, 
				   "Invalid CRC format; data ignored!");
		    }
		}
		val = crc;
	    }
	}
	return val;
    }

    public void setValueAt (Object value, int row, int col) {
	if (value == null || !(value instanceof Molecule)) {
	    return;
	}

	if (col == 1) {
	    Molecule mol = molseq.get(rows[row]);
	    if (mol != null) {
		((Molecule)value).clonecopy(mol);
	    }
	}
	else if (col < xgroupCount+2){
		Molecule xg = xgroups[row][col-2];
		if (xg != null) {
			((Molecule)value).clonecopy(xg);
		 }
	}
	else if (col < rgroupCount+2+PADDING_COLUMNS_XGROUP) {
	    Molecule rg = rgroups[row][col-2-PADDING_COLUMNS_XGROUP];
	    if (rg != null) {
		((Molecule)value).clonecopy(rg);
	    }
	}
    }

    CRC[] parseCRC (String str) {
	//System.out.println("parsing CRC => " + str);
	String[] lines = str.split(";");
	CRC[] empty = new CRC[0];
	if (lines.length == 0) {
		return empty;
	}
	int tok = 0;
	System.out.println("line "+tok+": "+lines[tok]);
	String[] tokens = lines[tok].trim().split("\\s");
	int nc = Integer.parseInt(tokens[1]);
	System.out.println(nc);
	CRC[] crc = new CRC[nc];
	for (int i = 0; i < nc; ++i) {
	    tokens = lines[++tok].trim().split("="); // length
	    if (!tokens[0].equals("length")) return empty;
	    CRC c = new CRC (Integer.parseInt(tokens[1]));

	    tokens = lines[++tok].trim().split("="); // zero
	    if (!tokens[0].equals("zero")) return empty;
	    if (!tokens[1].equals("null")) {
		c.zero = Double.parseDouble(tokens[1]);
	    }

	    tokens = lines[++tok].trim().split("="); // inf
	    if (!tokens[0].equals("inf")) return empty;
	    if (!tokens[1].equals("null")) {
		c.inf = Double.parseDouble(tokens[1]);
	    }

	    tokens  = lines[++tok].trim().split("="); // half
	    if (!tokens[0].equals("half")) return empty;
	    if (!tokens[1].equals("null")) {
		c.logac50 = Double.parseDouble(tokens[1]);
	    }

	    tokens = lines[++tok].trim().split("="); // slope
	    if (!tokens[0].equals("slope")) return empty;
	    if (!tokens[1].equals("null")) {
		c.slope = Double.parseDouble(tokens[1]);
	    }
	    
	    for (int j = 0; j < c.size; ++j) {
		tokens = lines[++tok].trim().split("\\s");
		c.xy[0][j] = Double.parseDouble(tokens[0]);
		c.xy[1][j] = Double.parseDouble(tokens[1]);
	    }
	    crc[i] = c;
	}
	return crc;
    }

    public boolean isCellEditable (int row, int col) { 
	return getColumnClass (col) == Molecule.class;
    }
    public double calculateScore(){
    	double score=0;
    	double sig = 0., size = core.getAtomCount();
        for (int i = 0; i < rows.length; ++i) {
            Molecule m = molseq.get(rows[i]);
            double x = m.getAtomCount() - size;
            sig += x*x;
        }
        if (sig > 0.) {
            
            //score = Math.log10(rgd.rows.length)* Math.sqrt(snr);
            score = -Math.log10
                (Math.sqrt(Math.pow(size,1)*//bigger is usually better
                           Math.pow(((double)rows.length)/molseq.size(),1)*//coverage is better
                           Math.pow(sig,-2)*//the members should be fairly close to the scaffold
                           Math.pow(getRGroupCount(),-1)//and there shouldn't be many r-groups
                ));
        }
        this.score=score;
        return score;
    }

    public boolean isVisible(int col){
    	return invisibleColumns.contains(col);
    }
}
