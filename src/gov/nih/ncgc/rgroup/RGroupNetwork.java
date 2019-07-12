package gov.nih.ncgc.rgroup;

import java.beans.*;
import java.util.*;
import java.util.List;
import java.io.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.net.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.border.*;

import org.jdesktop.swingworker.SwingWorker;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.DPoint3;
import chemaxon.sss.search.MolSearch;
import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.MolPrinter;
import chemaxon.marvin.paint.DispOptConsts;
import chemaxon.marvin.util.MolExportException;
import chemaxon.marvin.beans.MSketchPane;

import gov.nih.ncgc.viz.MViewRenderer;

import prefuse.*;
import prefuse.action.*;
import prefuse.action.filter.VisibilityFilter;
import prefuse.action.animate.*;
import prefuse.action.assignment.*;
import prefuse.action.layout.*;
import prefuse.action.layout.graph.*;
import prefuse.activity.*;
import prefuse.controls.*;
import prefuse.data.*;
import prefuse.data.util.Sort;
import prefuse.data.event.TupleSetListener;
import prefuse.data.io.GraphMLReader;
import prefuse.data.query.SearchQueryBinding;
import prefuse.data.search.PrefixSearchTupleSet;
import prefuse.data.search.SearchTupleSet;
import prefuse.data.tuple.DefaultTupleSet;
import prefuse.data.tuple.TupleSet;
import prefuse.data.expression.*;
import prefuse.render.*;
import prefuse.render.Renderer;
import prefuse.util.*;
import prefuse.util.display.DisplayLib;
import prefuse.util.ui.*;
import prefuse.visual.*;
import prefuse.visual.expression.InGroupPredicate;
import prefuse.visual.sort.TreeDepthItemSorter;

public class RGroupNetwork extends Display implements ActionListener {
    private static final Logger logger = 
        Logger.getLogger(RGroupNetwork.class.getName());

    /**
     * Data groups
     */
    static final String Graph = "graph";
    static final String GraphNodes = "graph.nodes";
    static final String GraphEdges = "graph.edges";
    static final String GraphAttachments = "graph.attachments";
    static final String GraphScaffolds = "graph.scaffolds";
    static final String Linear = "linear";


    static final String RGroup = "rgroup"; // r-group table
    static final String Score = "score";
    static final String Complexity = "complexity";
    static final String Scaffold = "scaffold";
    static final String Attachment = "attachment";
    static final String Structure = "structure";
    static final String Label = "label";
    static final String RootDistance = "rootdistance";

    /**
     * Set node fill colors
     */
    static class NodeColorAction extends ColorAction {
        public NodeColorAction(String group) {
            super(group, VisualItem.FILLCOLOR, ColorLib.gray(220,230));
            add("ingroup('_hover')", ColorLib.rgba(200,200,200,200));
            add("ingroup('_search_')", ColorLib.rgb(0, 0, 255));
            add("ingroup('_focus_')", ColorLib.rgba(198,229,229, 220));
        }
    } // end of inner class NodeColorAction
    
    static class NodeStrokeColorAction extends ColorAction {
	public NodeStrokeColorAction (String group) {
	    super (group, VisualItem.STROKECOLOR, ColorLib.rgb(200,200,200));
	    add (VisualItem.HIGHLIGHT, ColorLib.rgb(255,200,125));
            add ("ingroup('_hover')", ColorLib.rgb(200,200,200));
	}
    }

    static class EdgeStrokeColorAction extends ColorAction {
	static int highlightColor = ColorLib.rgb(255,200,125);

	public EdgeStrokeColorAction (String group) {
	    super (group, VisualItem.STROKECOLOR, ColorLib.rgb(20,20,20));
	    add("ingroup('_search_')", ColorLib.rgb(0, 0, 255));
	}

	public int getColor (VisualItem item) {
	    if (item.isHighlighted() || item.isHover()) {
		return highlightColor;
	    }

            return super.getColor(item);
        }
    }

    static class EdgeStrokeAction extends StrokeAction {
	public EdgeStrokeAction (String group) {
	    super (group);
	    add (VisualItem.HIGHLIGHT, new BasicStroke (2.f));
	}
    }

    /**
     * Set node text colors
     */
    static class TextColorAction extends ColorAction {
        public TextColorAction(String group) {
            super(group, VisualItem.TEXTCOLOR, ColorLib.gray(0));
            add("_hover", ColorLib.rgb(255,0,0));
        }
    } // end of inner class TextColorAction

    /**
     * Switch the root of the tree by requesting a new spanning tree
     * at the desired root
     */
    static class TreeRootAction extends GroupAction {
        public TreeRootAction (String graphGroup) {
            super(graphGroup);
        }
        public void run(double frac) {
            TupleSet focus = m_vis.getGroup(Visualization.FOCUS_ITEMS);
            if ( focus==null || focus.getTupleCount() == 0 ) return;
            
            Graph g = (Graph)m_vis.getGroup(m_group);
            Node f = null;
            Iterator tuples = focus.tuples();
            while (tuples.hasNext()) {
		Tuple tuple = (Tuple)tuples.next();
		if (!(tuple instanceof Node)) {
		    return;
		}
		if (!g.containsTuple(f = (Node)tuple)) {
		    f = null;
		}
	    }
            if ( f == null ) return;
            g.getSpanningTree(f);
        }
    }


    protected EdgeRenderer edgeRenderer;
    protected JPopupMenu popup;
    protected VisualItem currentSelectedItem = null;
    protected RadialTreeLayout graphLayout;
    protected PropertyChangeSupport pcs = new PropertyChangeSupport (this);
    protected Graph network;
    int largestRgroupCount=0;
    List<Integer> largestCounts;
    public RGroupNetwork () {
        super (new Visualization ());
        initUI ();
    }

    protected void initUI () {
        popup = createPopupMenu ();
        //m_vis.setInteractive(treeEdges, null, false);

        // -- set up renderers --
        MViewRenderer nodeRenderer2 = new MViewRenderer ();
        //nodeRenderer.setMolField(Structure);
        nodeRenderer2.setRenderType
            (AbstractShapeRenderer.RENDER_TYPE_DRAW_AND_FILL);
        nodeRenderer2.setRoundedCorner(15, 15);
        edgeRenderer = new EdgeRenderer ();
        
      
        final DefaultRendererFactory rf = new DefaultRendererFactory (nodeRenderer2);
        rf.add(new InGroupPredicate (GraphEdges), edgeRenderer);

        final MViewRenderer nodeRenderer = new MViewRenderer ();
        nodeRenderer.setMolField(Scaffold);
        nodeRenderer.setShape(new Arc2D.Double
                              (0., 0., 60., 60., 0., 360., Arc2D.OPEN));
        nodeRenderer.setRenderType
            (AbstractShapeRenderer.RENDER_TYPE_DRAW_AND_FILL);
        
        RendererFactory rf2 = new RendererFactory(){
        	double minSize = 30.;
			@Override
			public Renderer getRenderer(VisualItem arg0) {
				RGroupTable rgt=(RGroupTable)arg0.get(RGroup);
				if(rgt != null){
					double frac=Math.sqrt(((double)(rgt.getRowCount()))/largestRgroupCount);
					nodeRenderer.setShape(new Arc2D.Double
                            (0., 0., minSize+frac*120, minSize+frac*120, 0., 360., Arc2D.OPEN));
					return nodeRenderer;
				}
				return rf.getRenderer(arg0);
			}
        };
        rf.add(new AbstractPredicate () {
                public boolean getBoolean (Tuple t) {
                    return t.get(Scaffold) != null;
                }
            }, nodeRenderer);

        m_vis.setRendererFactory(rf2);
               
        // -- set up processing actions --
        
        // colors
        ItemAction nodeColor = new NodeColorAction(GraphNodes);
        ItemAction textColor = new TextColorAction(GraphNodes);
        m_vis.putAction("textColor", textColor);
        
        ItemAction edgeHighlightColor = 
            new EdgeStrokeColorAction (GraphEdges);
        ItemAction nodeHighlightColor = 
            new NodeStrokeColorAction (GraphNodes);
        ItemAction edgeStroke = new EdgeStrokeAction (GraphEdges);
        
        FontAction fonts = new FontAction
            (GraphNodes, FontLib.getFont("Tahoma", 10));
        fonts.add("ingroup('_focus_')", FontLib.getFont("Tahoma", 11));

        // recolor
        ActionList recolor = new ActionList();
        recolor.add(edgeStroke);
        recolor.add(nodeColor);
        recolor.add(edgeHighlightColor);
        recolor.add(nodeHighlightColor);
        recolor.add(textColor);
        m_vis.putAction("recolor", recolor);
        
        // repaint
        ActionList repaint = new ActionList();
        repaint.add(recolor);
        repaint.add(new RepaintAction());
        m_vis.putAction("repaint", repaint);
        
        // animate paint change
        ActionList animatePaint = new ActionList(400);
        animatePaint.add(new ColorAnimator(GraphNodes));
        animatePaint.add(new RepaintAction());
        m_vis.putAction("animatePaint", animatePaint);
        
        // create the tree layout action
        graphLayout = new RadialTreeLayout (Graph, 150);
        //graphLayout.setAngularBounds(-Math.PI/2, Math.PI);
        m_vis.putAction("graphLayout", graphLayout);
        
        Layout subLayout = 
            new CollapsedSubtreeLayout (Graph);
            //new BalloonTreeLayout (Graph, 100);
            //new NodeLinkTreeLayout (Graph);
            //new RadialTreeLayout(Graph);
        //new FruchtermanReingoldLayout(Graph);
        m_vis.putAction("subLayout", subLayout);

        // recentre and rezoom on reload
        prefuse.action.Action resizeAction = new prefuse.action.Action() {
                public void run(double frac) {
                    // animate reset zoom to fit the data (must run only AFTER layout)
                    Rectangle2D bounds = m_vis.getBounds(Graph);
			
                    if (bounds.getWidth() == 0) return;
                    GraphicsLib.expand
                        (bounds, 10+(int)(1/m_vis.getDisplay(0).getScale()));
                    DisplayLib.fitViewToBounds
                        (m_vis.getDisplay(0), bounds, (long) 1250);
                }
            };
        m_vis.putAction("resize", resizeAction);

        
        // create the filtering and layout
        ActionList filter = new ActionList();
        filter.add(new TreeRootAction (Graph));
        filter.add(fonts);
        filter.add(graphLayout);
        filter.add(subLayout);
        filter.add(textColor);
        filter.add(nodeColor);
        filter.add(edgeHighlightColor);
        filter.add(nodeHighlightColor);
        filter.add(edgeStroke);
        m_vis.putAction("filter", filter);


        // animated transition
        ActionList animate = new ActionList (1250);
        animate.setPacingFunction(new SlowInSlowOutPacer ());
        animate.add(new QualityControlAnimator ());
        animate.add(new VisibilityAnimator (Graph));
        animate.add(new PolarLocationAnimator (GraphNodes, Linear));
        animate.add(new ColorAnimator (GraphNodes));
        animate.add(new RepaintAction());
        m_vis.putAction("animate", animate);
        m_vis.alwaysRunAfter("filter", "animate");
        m_vis.alwaysRunAfter("animate", "resize");
        
        // ------------------------------------------------
        
        // initialize the display
        //setSize(600,600);
        setItemSorter(new TreeDepthItemSorter(true));

        addControlListener(new DragControl());
        addControlListener(new ZoomToFitControl());
        addControlListener(new ZoomControl());
        addControlListener(new PanControl());
        addControlListener(new WheelZoomControl ());
        addControlListener(new FocusControl(1, "filter"));
        addControlListener(new HoverActionControl("repaint"));
        addControlListener(new NeighborHighlightControl () {
                @Override
                public void itemEntered(VisualItem item, MouseEvent e) {
                    doHighlight (item, true);
                }

                @Override
                public void itemExited(VisualItem item, MouseEvent e) {
                    doHighlight (item, false);
                }

                @Override
                public void itemPressed (VisualItem item, MouseEvent e) {
                    if (currentSelectedItem != null) {
                        doHighlight (currentSelectedItem, false);
                    }
                    doHighlight (item, true);
                }

                void doHighlight (VisualItem item, boolean highlight) {
                    if (item instanceof NodeItem) {
                        setNeighborHighlight((NodeItem)item, highlight);
                    }
                    else if (item instanceof EdgeItem) {
                        item.setHighlighted(highlight);
                        ((EdgeItem)item).getSourceItem()
                            .setHighlighted(highlight);
                        ((EdgeItem)item).getTargetItem()
                            .setHighlighted(highlight);
                    }
                }
            });
        addControlListener(new ToolTipControl ("name") {
                public void itemEntered (VisualItem item, MouseEvent e) {
                    Display dis = (Display)e.getSource();
                }
            });
        addControlListener (new ControlAdapter () {
                public void itemPressed (VisualItem item, MouseEvent e) {
                    VisualItem old = currentSelectedItem;
                    currentSelectedItem = item;

                    if (e.isPopupTrigger()) {
                        showPopup (e);
                    }
                    else {
                        RGroupTable rg = (RGroupTable)item.get(RGroup);
                        if (rg != null) {
                            RGroupTable tab = old != null 
                                ? (RGroupTable)old.get(RGroup) : null;
                            pcs.firePropertyChange("rgroup", tab, rg);
                            logger.info("## selected "+rg);
                        }
                    }
                }

                public void itemEntered (VisualItem item, MouseEvent e) {
                    //logger.info(text);
                    //pcs.firePropertyChange("hover", null, text);
                    
                }

                public void itemExited (VisualItem item, MouseEvent e) {
                }
            });
        
        // ------------------------------------------------
        // maintain a set of items that should be interpolated linearly
        // this isn't absolutely necessary, but makes the animations nicer
        // the PolarLocationAnimator should read this set and act accordingly
        m_vis.addFocusGroup(Linear, new DefaultTupleSet());
        m_vis.getGroup(Visualization.FOCUS_ITEMS).addTupleSetListener
            (new TupleSetListener() {
                    public void tupleSetChanged
                        (TupleSet t, Tuple[] add, Tuple[] rem) {

                        TupleSet linearInterp = m_vis.getGroup(Linear);
                        if ( add.length < 1 ) return; 
                        linearInterp.clear();
                        if (add[0] instanceof Node) {
                            for ( Node n = (Node)add[0]; 
                                  n!=null; n=n.getParent() )
                                linearInterp.addTuple(n);
                        }
                    }
                });

        addComponentListener (new ComponentAdapter () {
                public void componentResized (ComponentEvent e) {
                    m_vis.run("resize");
                }
            });
    }

    
    JPopupMenu createPopupMenu () {
        JPopupMenu menu = new JPopupMenu ();
        JMenuItem item;
        menu.add(item = new JMenuItem ("Hide"));
        item.setToolTipText("Hide node");
        item.addActionListener(this);

        menu.addSeparator();
        menu.add(item = new JMenuItem ("Copy to editor"));
        item.setToolTipText("Copy this structure to the editor");
        item.addActionListener(this);

        return menu;
    }

    void showPopup (MouseEvent e) {
        popup.show(this, e.getX(), e.getY());
    }

    public void actionPerformed (ActionEvent e) {
        String cmd = e.getActionCommand().toLowerCase();
        logger.info("Action "+cmd);

        if (currentSelectedItem == null) {
            logger.info("No item selected!");
            return;
        }
	    
        if (cmd.equals("hide")) {
            // hide all edge too
            SwingUtilities.invokeLater(new Runnable () {
                    public void run () {
                        NodeItem node = (NodeItem)currentSelectedItem;
                        for (Iterator iter = node.edges();  
                             iter.hasNext(); ) {
                            EdgeItem edge = (EdgeItem)iter.next();
                            edge.setVisible(false);
                        }
                        currentSelectedItem.setVisible(false);
                        RGroupTable rgt;
                        if((rgt=(RGroupTable) node.get(RGroup))!=null){
                        	if(rgt.getRowCount()>=largestRgroupCount){
                        		largestCounts.remove(largestCounts.size()-1);
                        		largestRgroupCount=largestCounts.get(largestCounts.size()-1);
                        	}
                        }
                        currentSelectedItem = null;
                    }
                });
        }
        else if (cmd.startsWith("copy")) {
            Molecule mol = (Molecule)currentSelectedItem.get("_structure");
            if (mol == null) {
                mol = (Molecule)currentSelectedItem.get("_frag");
            }

            if (mol != null) {
                mol = mol.cloneMolecule();
                mol.aromatize();
            }
            pcs.firePropertyChange("copyMol", null, mol);
        }
    }

    void selectItem (VisualItem item, MouseEvent e) {
        currentSelectedItem = item;
        logger.info("Item selected: "+item);
    }

    protected void setGraph (Graph network) {
        this.network = network;
        m_vis.reset();
        m_vis.addGraph(Graph, network);
        m_vis.run("filter");
    }

    public void reload () {
        Graph g = (Graph)m_vis.getSourceData(Graph);
        setGraph (g);
    }

    protected void run (String action) {
        m_vis.run(action);
    }

    public void setRadialRadius (int radius) {
        graphLayout.setRadiusIncrement(radius);
        TupleSet ts = m_vis.getSourceData(Graph);
        if (ts != null && ts.getTupleCount() > 0) {
            m_vis.run("filter");
        }
    }

    public void unhideAll () {
        SwingUtilities.invokeLater(new Runnable () {
                public void run () {
                    VisualGraph vg = 
                        (VisualGraph)m_vis.getVisualGroup(Graph);
                    for (Iterator iter = vg.nodes(); iter.hasNext(); ) {
                        NodeItem n = (NodeItem)iter.next();
                        n.setVisible(true);
                    }
                    for (Iterator iter = vg.edges(); iter.hasNext(); ) {
                        EdgeItem e = (EdgeItem)iter.next();
                        e.setVisible(true);
                    }
                }
            });
    }

    protected Visualization getVis () { return m_vis; }

    public void addPropertyChangeListener (PropertyChangeListener l) {
        pcs.addPropertyChangeListener(l);
    }

    public void removePropertyChangeListener (PropertyChangeListener l) {
        pcs.removePropertyChangeListener(l);
    }

    public void setRGroups (final RGroupTable... rgroups) {
        SwingUtilities.invokeLater(new Runnable () {
                public void run () {
                    setGraph (createNetwork3 (rgroups));
                }
            });
    }

    // a simple hash for the attachment
    static int getAttachmentHash (Molecule mol) {
        BitSet bs = new BitSet ();
        for (MolAtom a : mol.getAtomArray()) {
            int h = a.getAtno();
            bs.set(h);
            h *= 31;
            for (int i = 0; i < a.getBondCount(); ++i) {
                MolAtom xa = a.getBond(i).getOtherAtom(a);
                bs.set(h+a.getBond(i).getType()*xa.getAtno());
            }
        }
        return bs.hashCode();
    }

    /**
     * This network consists of scaffolds connecting to their r-group
     * unique substituents
     */
    protected Graph createNetwork1 (RGroupTable... rgroups) {
        Graph network = new Graph ();
        network.addColumn(RGroup, RGroupTable.class);
        network.addColumn(Score, double.class);
        network.addColumn(Complexity, int.class);
        network.addColumn(Scaffold, Molecule.class);
        network.addColumn(Attachment, Molecule.class);

        Map<Integer, Set<Node>> adj = 
            new HashMap<Integer, Set<Node>>();
        Map<Integer, Molecule> attachments = new HashMap<Integer, Molecule>();

        for (RGroupTable rg : rgroups) {
            Node node = network.addNode();
            node.set(RGroup, rg);
            node.set(Score, rg.getScaffoldScore());
            node.set(Complexity, rg.getScaffoldComplexity());
            node.set(Scaffold, rg.getScaffold());
            for (int i = 0; i < rg.getRowCount(); ++i) {
                for (int j = 0; j < rg.getRGroupCount(); ++j) {
                    Molecule rij = rg.getRGroup(i, j);
                    if (rij != null) {
                        int hash = getAttachmentHash (rij);
                        Set<Node> parents = adj.get(hash);
                        if (parents == null) {
                            adj.put(hash, parents = new HashSet<Node>());
                            attachments.put(hash, rij);
                        }
                        parents.add(node);
                    }
                }
            }
        }

        // now create the nodes for the attachments
        for (Map.Entry<Integer, Set<Node>> me : adj.entrySet()) {
            try {
                MolHandler mh = new MolHandler (attachments.get(me.getKey()));
                Molecule att = mh.getMolecule();
                Node node = network.addNode();
                node.set(Attachment, att);
                for (Node parent : me.getValue()) {
                    Edge e = network.addEdge(node, parent);
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        return network;
    }

    protected Graph createNetwork2 (RGroupTable... rgroups) {
        Graph network = new Graph ();
        network.addColumn(RGroup, RGroupTable.class);
        network.addColumn(Score, double.class);
        network.addColumn(Complexity, int.class);
        network.addColumn(Scaffold, Molecule.class);
        network.addColumn(Attachment, Molecule.class);
        network.addColumn(Structure, Molecule.class);

        Map<Molecule, Node> members = new HashMap<Molecule, Node>();

        for (RGroupTable rg : rgroups) {
            Node node = network.addNode();
            node.set(RGroup, rg);
            node.set(Score, rg.getScaffoldScore());
            node.set(Complexity, rg.getScaffoldComplexity());
            node.set(Scaffold, rg.getScaffold());

            for (Molecule m : rg.getMembers()) {
                Node n = members.get(m);
                if (n == null) {
                    members.put(m, n = network.addNode());
                    n.set(Structure, m);
                }
                Edge e = network.addEdge(node, n);
            }

        }

        return network;
    }
    /**
     * This network consists of scaffolds connecting to other scaffolds
     * based on sub/super structure. It's directed.
     */
    protected Graph createNetwork3 (RGroupTable... rgroups) {
        Graph network = new Graph (true);
        network.addColumn(RGroup, RGroupTable.class);
        network.addColumn(Score, double.class);
        network.addColumn(Complexity, int.class);
        network.addColumn(Scaffold, Molecule.class);
        network.addColumn(Attachment, Molecule.class);
        network.addColumn(RootDistance, int.class);
        largestCounts=new ArrayList<Integer>();
        Map<Integer,Node> nodes = new HashMap<Integer,Node>();
        largestRgroupCount=0;
        for (int k=0;k<rgroups.length;k++) {
        	RGroupTable rg = rgroups[k];
            Node node = network.addNode();
            node.set(RGroup, rg);
            node.set(Score, rg.getScaffoldScore());
            node.set(Complexity, rg.getScaffoldComplexity());
            node.set(Scaffold, rg.getScaffold());
            largestRgroupCount=Math.max(largestRgroupCount, rg.getRowCount());
            largestCounts.add(rg.getRowCount());
            nodes.put(k,node);
        }
        Collections.sort(largestCounts);
        for (int i=0;i<rgroups.length;i++) {
        	RGroupTable rg = rgroups[i];
        	int[] fp = rg.getScaffoldFp();

        	int scount = 0;
        	for(int k=0;k<fp.length;k++){
				scount+=Integer.bitCount(fp[k]);
			}
        	
        	double minDist = Double.MAX_VALUE;
        	int minIndex=-1;
        	for(int j=0;j<rgroups.length;j++){
        		if(i!=j){
        			RGroupTable rg2 = rgroups[j];
        			int[] fp2 = rg2.getScaffoldFp();
        			int dcount = 0;
        			int ccount = 0;
        			int ocount = 0;
        			for(int k=0;k<fp2.length;k++){
        				ccount+=Integer.bitCount(fp[k]&fp2[k]);
        				ocount+=Integer.bitCount(fp2[k]);
        				dcount+=Integer.bitCount(fp[k]|fp2[k]);
        			}
        			double dist = 1 - ((double)ccount/(double)dcount);
        			if(dist<minDist){
        				minDist= dist;
        				minIndex=j;
        			}
        			if(ccount==scount && ccount!=ocount){
        				Edge e = network.addEdge(nodes.get(i), nodes.get(j));
        			}
        		}
        	}
        	//Edge e = network.addEdge(nodes.get(i), nodes.get(minIndex));
        }
        pruneGraph(network);
        return network;
    }
    public static void pruneGraph(Graph g){
    	Stack<Node> grandparents = new Stack<Node>();
    	Set<Edge> toRemove = new HashSet<Edge>();
    	Set<Node> roots = new HashSet<Node>();
    	
    	for(int i=0;i<g.getNodeCount();i++){
    		Node n=g.getNode(i);
    		if(n.getInDegree()==0){
    			roots.add(n);
    		}
    	}
    	Node realRoot = g.addNode();
    	for(Node n:roots){
    		g.addEdge(realRoot,n);
    	}
    	pruneGraph(g,grandparents,realRoot,toRemove);
    	/*
    	 for(Node n:roots){
    			
    			pruneGraph(g,grandparents,n,toRemove);
    			try {
					System.out.println("PARENT NODE" + ((Molecule)n.get(Scaffold)).exportToFormat("cxsmarts:q"));
				} catch (MolExportException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}*/
    	//}
    	for(Edge e:toRemove){
    		g.removeEdge(e);
    	}
    	
    	
    		
    }
    public static void pruneGraph(Graph g, Stack<Node> grandparents,Node parent,Set<Edge> toremove){
    	Integer mdist =(Integer) parent.get(RootDistance);
    	if(grandparents.contains(parent)){
    		System.out.println("NO SENSE!!!!");
    		for(Node n : grandparents){
    			try {
					System.out.println("CHILD:" + ((Molecule)n.get(Scaffold)).exportToFormat("cxsmarts:q"));
				} catch (MolExportException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
    		}
    	}
    	if(mdist==null)
    		mdist=0;
    	parent.set(RootDistance,Math.max(grandparents.size(),mdist));
    	
    	Node child;
    	Iterator<Node> neighbors=parent.outNeighbors();
    	for(;neighbors.hasNext();){
    		child=neighbors.next();
    		for(Node j:grandparents){
        		Edge e= g.getEdge(child, j);
        		if(e!=null){
        			toremove.add(e);
        		}
        	}
    		
    		grandparents.push(parent);
    		pruneGraph(g, grandparents,child,toremove);
    		grandparents.pop();
    	}
    }
}
