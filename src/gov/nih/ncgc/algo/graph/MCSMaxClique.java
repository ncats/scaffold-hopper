// $Id: MCSMaxClique.java 3567 2009-11-15 06:22:20Z nguyenda $

package gov.nih.ncgc.algo.graph;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;
import java.util.logging.Level;

import java.io.PrintStream;
import java.io.OutputStream;
import java.io.IOException;
import java.io.FileOutputStream;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolFormatException;
import chemaxon.util.MolHandler;

import gov.nih.ncgc.util.ChemUtil;
import gov.nih.ncgc.util.AtomComparator;
import gov.nih.ncgc.util.DefaultAtomComparator;

public class MCSMaxClique extends MCSAbstract {
    static private final Logger logger = Logger.getLogger
        (MCSMaxClique.class.getName());

    // minimum clique size
    static final int MIN_CLIQUE = 2;

    // debugging level
    static private int DEBUG = 0;
    static {
        try {
            DEBUG = Integer.getInteger("maxclique.debug", 0);
        }
        catch (Exception ex) {}
    }

    static enum MatchDirection {
        FORWARD,
            REVERSE,
            BIDIR,
            NONE
            };

    // association node
    static class ANode {
        // n1 & n2 are bonds belonging to G1 (query) and G2 (target),
        //   respecitvely
        public int n1, n2;
        // (u1,v1) is an edge in G1 that corresponds to bond n1
        // similarly for (u2,v2)
        public int u1, v1, u2, v2;
        public MatchDirection dir;

        public ANode (int n1, int n2, int u1, int v1, int u2, int v2,
                      MatchDirection dir) {
            this.n1 = n1;
            this.n2 = n2;
            this.u1 = u1;
            this.v1 = v1;
            this.u2 = u2;
            this.v2 = v2;
            this.dir = dir;
        }

        public String toString () {
            String d = "?";
            switch (dir) {
            case FORWARD: d = "+"; break;
            case REVERSE: d = "-"; break;
            case BIDIR: d = "="; break;
            }
            return "[" + (u1+1) +"," + (v1+1)
                + "]" + d + "[" + (u2+1) + "," + (v2+1) + "]";
        }
    }

    /**
     * Maximal clique
     */
    public class MClique implements Comparable<MClique> {
        AGraph graph; // reference the corresponding graph
        BitSet component; // input subset of graph vertices
        int size; // clique size
        Molecule core; // corresponding core
        long elapsed; // time spent 

        // all maximal cliques
        Map<BitSet, int[]> cliques = new HashMap<BitSet, int[]>(); 
        
        MClique (AGraph graph, BitSet component) {
            this.graph = graph;
            this.component = component;
        }

        public int size () { return size; } 
        public Molecule core () { return core; }
        public int[][] getHits () {
            return cliques.values().toArray(new int[0][]);
        }

        public int compareTo (MClique mc) {
            return mc.core.getBondCount() - core.getBondCount();
        }

        boolean isClique (int vertex) {
            for (BitSet bv : cliques.keySet()) 
                if (bv.get(vertex)) 
                    return true;
                
            return false;
        }

        public String toString () {
            StringBuilder sb = new StringBuilder ("MClique {\n");
            //sb.append(" component: "+component.cardinality()+component+"\n");
            sb.append(" size: "+size+"\n");
            sb.append(" core: "+core.getAtomCount()+"/"
                      +core.getBondCount()+" "+core.toFormat("smiles:q")+"\n");
            sb.append(" elapsed: "+elapsed+"ms\n");
            sb.append(" cliques: "+cliques.size()+" = [\n");
            for (Map.Entry<BitSet, int[]> me : cliques.entrySet()) {
                sb.append("  ");
                int[] v = me.getValue();
                for (int i = 0; i < v.length; ++i) {
                    if (v[i] >= 0) {
                        sb.append(" "+(i+1)+":"+(v[i]+1));
                    }
                }
                sb.append("\n");
            }
            sb.append(" ]\n}\n");
            return sb.toString();
        }

        public void debugG (String name, OutputStream os) throws IOException {
            PrintStream ps = new PrintStream (os);  
            ps.println(" <graph id=\""+name+"\" edgedefault=\"undirected\">");
            for (int i = 0; i < graph.size(); ++i) {
                ANode n = graph.vertex(i);
                ps.println("  <node id=\"n"+(i+1)+"\">");
                ps.println("   <data key=\"u1\">"+(n.u1+1)+"</data>");
                ps.println("   <data key=\"v1\">"+(n.v1+1)+"</data>");
                ps.println("   <data key=\"u2\">"+(n.u2+1)+"</data>");
                ps.println("   <data key=\"v2\">"+(n.v2+1)+"</data>");
                ps.println("   <data key=\"b1\">"
                           +(query.getBond(n.n1).getType())+"</data>");
                ps.println("   <data key=\"b2\">"
                           +(target.getBond(n.n2).getType())+"</data>");
                ps.println("   <data key=\"ns\">"+(graph.mapped(i)?1:0)
                           +"</data>");
                ps.println("   <data key=\"nc\">"+(isClique (i)?1:0)
                           +"</data>");
                ps.println("  </node>");
            }

            // do d-edge
            debugD (ps);

            // do c-edge
            debugC (ps);

            ps.println(" </graph>");
            ps.flush();
        }

        public void debugC (String name, OutputStream os) throws IOException {
            PrintStream ps = new PrintStream (os);  
            ps.println(" <graph id=\""+name+"\" edgedefault=\"undirected\">");
            for (int i = 0; i < graph.size(); ++i) {
                ANode n = graph.vertex(i);
                if (!graph.cEdges(i).isEmpty()) {
                    ps.println("  <node id=\"n"+(i+1)+"\">");
                    ps.println("   <data key=\"u1\">"+(n.u1+1)+"</data>");
                    ps.println("   <data key=\"v1\">"+(n.v1+1)+"</data>");
                    ps.println("   <data key=\"u2\">"+(n.u2+1)+"</data>");
                    ps.println("   <data key=\"v2\">"+(n.v2+1)+"</data>");
                    ps.println("   <data key=\"b1\">"
                               +(query.getBond(n.n1).getType())+"</data>");
                    ps.println("   <data key=\"b2\">"
                               +(target.getBond(n.n2).getType())+"</data>");
                    ps.println("   <data key=\"nc\">"+(isClique (i)?1:0)
                               +"</data>");
                    ps.println("   <data key=\"ns\">"+(graph.mapped(i)?1:0)
                               +"</data>");
                    ps.println("  </node>");
                }
            }

            // do c-edge
            debugC (ps);

            ps.println(" </graph>");
        }

        public void debugC (PrintStream ps) {
            int edges = 0;
            for (int i = 0; i < graph.size(); ++i) {
                for (int j = graph.cEdges(i).nextSetBit(0); 
                     j >= 0; j = graph.cEdges(i).nextSetBit(j+1)) {
                    if (j > i) { // undirected
                        ps.println("  <edge source=\"n"
                                   +(i+1)+"\" target=\"n"+(j+1)+"\">");
                        ps.println("   <data key=\"type\">1</data>");
                        ps.println("   <data key=\"es\">"
                                   +(graph.mapped(i) 
                                     && graph.mapped(j) ? 1:0)+"</data>");
                        ps.println("   <data key=\"ec\">"
                                   +(isClique (i) && isClique (j)
                                     ? 1 : 0)+"</data>");
                        ps.println("  </edge>");
                    }
                }
                edges += graph.cEdges(i).cardinality();
            }
            logger.info("## c-edges="+edges);
        }

        public void debugD (String name, OutputStream os) throws IOException {
            PrintStream ps = new PrintStream (os);  
            ps.println(" <graph id=\""+name+"\" edgedefault=\"undirected\">");
            for (int i = 0; i < graph.size(); ++i) {
                if (!graph.dEdges(i).isEmpty()) {
                    ANode n = graph.vertex(i);
                    ps.println("  <node id=\"n"+(i+1)+"\">");
                    ps.println("   <data key=\"u1\">"+(n.u1+1)+"</data>");
                    ps.println("   <data key=\"v1\">"+(n.v1+1)+"</data>");
                    ps.println("   <data key=\"u2\">"+(n.u2+1)+"</data>");
                    ps.println("   <data key=\"v2\">"+(n.v2+1)+"</data>");
                    ps.println("   <data key=\"b1\">"
                               +(query.getBond(n.n1).getType())+"</data>");
                    ps.println("   <data key=\"b2\">"
                               +(target.getBond(n.n2).getType())+"</data>");
                    ps.println("  </node>");
                }
            }

            // do d-edge
            debugD (ps);
            ps.println(" </graph>");
        }

        public void debugD (PrintStream ps) {
            int edges = 0;
            for (int i = 0; i < graph.size(); ++i) {
                for (int j = graph.dEdges(i).nextSetBit(0); 
                     j >= 0; j = graph.dEdges(i).nextSetBit(j+1)) {
                    if (j > i) { // undirected
                        ps.println("  <edge source=\"n"
                                   +(i+1)+"\" target=\"n"+(j+1)+"\">");
                        ps.println("   <data key=\"type\">2</data>");
                        ps.println("   <data key=\"es\">"
                                   +(graph.mapped(i) && graph.mapped(j) 
                                     ? 1:0)+"</data>");
                        ps.println("   <data key=\"ec\">"
                                   +(isClique (i) && isClique (j)
                                     ? 1 : 0)+"</data>");
                        ps.println("  </edge>");
                    }
                }
                edges += graph.dEdges(i).cardinality();
            }
            logger.info("## d-edges="+edges);
        }
    } // MClique

    /**
     * Association graph
     */
    class AGraph {
        ANode[] V; // vertices of graph
        BitSet[] G; // adjancy graph; G = C U D
        BitSet[] C; // c-edge graph; C.length = V.length
        BitSet[] D; // d-edge graph; D.length = V.length
        BitSet[] components;

        // atom mapping
        Map<Integer, Integer> mapping;
        BitSet mapped; // mapped ANode's

        int cEdges, dEdges;

        AGraph (ANode[] V) {
            this (V, new HashMap<Integer, Integer>(), new BitSet (V.length));
        }

        AGraph (ANode[] V, Map<Integer, Integer> mapping, BitSet mapped) {
            G = new BitSet[V.length]; // association graph
            D = new BitSet[V.length]; // d-graph of G
            C = new BitSet[V.length]; // c-graph of G

            this.mapping = mapping;
            this.mapped = mapped;

            if (DEBUG > 1) {
                System.err.println("** Association graph");
            }

            for (int i = 0; i < V.length; ++i) {
                C[i] = new BitSet (V.length);
                D[i] = new BitSet (V.length);
                BitSet bs = new BitSet (V.length);
                for (int j = 0; j < V.length; ++j) {
                    if (i != j && isConnected (V[i], V[j])) {
                        if (isCEdge (V[i], V[j])) {
                            C[i].set(j);
                        }
                        else { // d-edge
                            D[i].set(j);
                        }
                        bs.set(j);
                    }
                }
                G[i] = bs;

                cEdges += C[i].cardinality();
                dEdges += D[i].cardinality();

                //System.out.println(i + ": " + bs.cardinality() + " " + bs);
                if (DEBUG > 1) {
                    System.err.println(String.format("[%1$3d]", i) 
                                       + ": C=" + C[i].cardinality() + C[i] 
                                       + " D=" + D[i].cardinality() + D[i]);
                }
            }

            BitSet all = new BitSet (V.length);
            all.set(0, V.length);

            List<BitSet> concomp = new ArrayList<BitSet>();
            for (int u = all.nextSetBit(0); u >= 0; u = all.nextSetBit(u+1)) {
                BitSet seen = new BitSet (V.length);
                dfs (u, C, seen);
                concomp.add((BitSet)seen.clone());
                all.andNot(seen);
            }

            Collections.sort(concomp, new Comparator<BitSet>() {
                                 public int compare (BitSet a, BitSet b) {
                                     return b.cardinality() - a.cardinality();
                                     //return a.cardinality() - b.cardinality();
                                 }
                             });
            components = concomp.toArray(new BitSet[0]);
            if (DEBUG > 0) {
                Molecule q = query.cloneMolecule();
                Molecule t = target.cloneMolecule();
                for (int i = 0; i < q.getAtomCount(); ++i)
                    q.getAtom(i).setAtomMap(i+1);
                for (int i = 0; i < t.getAtomCount(); ++i)
                    t.getAtom(i).setAtomMap(i+1);
                System.err.println("***  query: "+q.toFormat("smiles:q"));
                System.err.println("*** target: "+t.toFormat("smiles:q"));
                System.err.println
                    ("## "+components.length+" connected component(s):");
                for (int i = 0; i < components.length; ++i) {
                    System.err.print(">> " + components[i].cardinality());
                    for (int j = components[i].nextSetBit(0); j >= 0;
                         j = components[i].nextSetBit(j+1)) {
                        System.err.print(" "+j+":"+V[j]);
                    }
                    System.err.println();
                }
            }

            this.V = V;            
        }

        public boolean isNeighbor (int v, int u) {
            return G[v].get(u);
        }
        public ANode vertex (int n) { return V[n]; }
        public int edgeCount () { return cEdges+dEdges; }
        public int size () { return V.length; }
        public boolean isMapped (int n) { return mapped.get(n); }
        public int cEdges () { return cEdges; }
        public int dEdges () { return dEdges; }
        public BitSet[] components () { return components; }
        public BitSet mapped () { return mapped; }
        public boolean mapped (int n) { return mapped.get(n); }
        public BitSet dEdges (int vertex) { return D[vertex]; }
        public BitSet cEdges (int vertex) { return C[vertex]; }
        public BitSet edges (int vertex) { return G[vertex]; }
    } // class AGraph

    static class Edge {
        int u, v;
        Edge (int u, int v) {
            this.u = u;
            this.v = v;
        }
    }

    private AGraph G; // association/product graph
    private MolBond[] qBonds, tBonds;

    private boolean allcliques = false; // generate all cliques?
    private long maxtime = -1; // maxtime allowed (in milli seconds) 
    // should the mcs be evaluated by its chemical relevance?
    private boolean evalByScore = false; 

    /*
     * all maximal cliques and their resulting MCSes
     */ 
    private List<MClique> maxcliques = new ArrayList<MClique>();

    // generate atom maps in mcs output:
    //  Query - mapping from query molecule
    //  Target - mapping from target molecule
    //  None - no mapping in output
    public enum OutputAtomMap {
        Query (-1),
            Target (1),
            None (0);
        final int ord;
        OutputAtomMap (int ord) {
            this.ord = ord;
        }
        public int ord () { return ord; }
    }
    private OutputAtomMap outmap = OutputAtomMap.None;

    // should the output atom map be restricted to the seeded atoms
    // only (i.e., atoms that are specified via setMapping)
    private boolean atomMapFromMappingOnly = false;

    /*
     * threading
     */
    private final ReentrantLock lock = new ReentrantLock ();
    private int nthreads = 0; // default number of threads (0 - no threading)

    /*
     * temporary variables
     */
    private long start, stop; // start & stop time
    private int maxDepth = 0; // control the maximum total depth to explore

    private volatile long ctime = 0; // time spent searching for maxclique
    private volatile int cdepth = 0; // current max clique depth
    private volatile int tdepth = 0; // total depth

    private volatile BitSet maxclique = new BitSet ();
    private volatile Edge[] maxcpath; // max clique path (if any)
    private volatile int[] forward = {}, reverse = {}; // maximal mapping
    private volatile int maxcsize = 0; // maximum clique size found thus far
    // contain all maximal hits if allcliques is true
    private volatile Map<BitSet, int[]> hits = new HashMap<BitSet, int[]>(); 
    private volatile Molecule maxcore = null; // max mces core
    private volatile int maxscore = 0; // score of maxcore

    // user-specified mapping
    private Map<Integer, Integer> mapping = new HashMap<Integer, Integer>();
    // exclude query and/or target atoms from matching
    private BitSet xQuery, xTarget; 

    /*
     * Seeding (if any) for the mcs
     */
    private Molecule seed;

    public MCSMaxClique () {
    }

    public MCSMaxClique (MCSAbstract mcs) {
        super (mcs);
    }

    private MatchDirection getMatchDirection (MolBond a, MolBond b) {
        MatchDirection dir = MatchDirection.NONE;

        if (atomComparator.match(a.getAtom1(), b.getAtom1())
            && atomComparator.match(a.getAtom2(), b.getAtom2())) {
            dir = MatchDirection.FORWARD;
        }
        else if (atomComparator.match(a.getAtom1(), b.getAtom2())
                 && atomComparator.match(a.getAtom2(), b.getAtom1())) {
            dir = MatchDirection.REVERSE;
        }

        if (dir != MatchDirection.NONE
            && atomComparator.match(a.getAtom1(), a.getAtom2())) {
            dir = MatchDirection.BIDIR;
        }

        return dir;
    }

    static private boolean isConnected (MolBond a, MolBond b) {
        return getCommonAtom (a, b) != null;
    }

    static private MolAtom getCommonAtom (MolBond a, MolBond b) {
        MolAtom a1 = a.getAtom1(), a2 = a.getAtom2();
        MolAtom b1 = b.getAtom1(), b2 = b.getAtom2();
        if (a1.equals(b1) || a1.equals(b2)) return a1;
        if (a2.equals(b1) || a2.equals(b2)) return a2;
        return null;
    }

    private boolean isConnected (ANode a, ANode b) {
        MolAtom u = getCommonAtom (qBonds[a.n1], qBonds[b.n1]);
        MolAtom v = getCommonAtom (tBonds[a.n2], tBonds[b.n2]);
        return ((u != null && v != null && atomComparator.match(u,v))
                || (u == null && v == null));
    }

    private boolean isCEdge (ANode a, ANode b) {
        return (isConnected (qBonds[a.n1], qBonds[b.n1])
                || isConnected (tBonds[a.n2], tBonds[b.n2]));
    }


    private int pickVertex (BitSet P) {
        int maxc = -1, max = -1, vmax = -1;
        for (int v = P.nextSetBit(0); v >= 0; v = P.nextSetBit(v+1)) {
            // have preference for vertex with higher c-edge
            int c = G.cEdges(v).cardinality();
            if (c > maxc) {
                maxc = c;
                max = G.edges(v).cardinality(); // update total max
                vmax = v;
            }
            else if (c == maxc && G.edges(v).cardinality() > max) {
                max = G.edges(v).cardinality();
                vmax = v;
            }
        }
        return vmax;
    }

    private int pickMinVertex (BitSet P) {
        int min = 1000000, vmin = -1;
        for (int v = P.nextSetBit(0); v >= 0; v = P.nextSetBit(v+1)) {
            int c = G.edges(v).cardinality();
            if (c < min) {
                min = c;
                vmin = v;
            }
        }
        return vmin;
    }

    private ANode[] createANodesNoMapping () {
        List<ANode> nodes = new ArrayList<ANode>();

        for (int i = 0; i < qBonds.length; ++i) {
            MolBond qb = qBonds[i];
            int u1 = query.indexOf(qb.getAtom1());
            int v1 = query.indexOf(qb.getAtom2());
            if (xQuery != null && (xQuery.get(u1) || xQuery.get(v1))) {
                // exclude this from query
            }
            else {
                for (int j = 0; j < tBonds.length; ++j) {
                    MolBond tb = tBonds[j];
                    int u2 = target.indexOf(tb.getAtom1());
                    int v2 = target.indexOf(tb.getAtom2());
                    if (xTarget != null
                        && (xTarget.get(u2) || xTarget.get(v2))) {
                        // exclude this
                    }
                    else if (bondComparator.match(qb, tb)) {
                        MatchDirection dir = getMatchDirection (qb, tb);
                        if (dir != MatchDirection.NONE) {
                            ANode n = new ANode
                                (i, j, u1, v1, u2, v2, dir);
                            nodes.add(n);
                        }
                    }
                }
            }
        } // endfor each query bond

        return nodes.toArray(new ANode[0]);
    }


    private ANode[] createANodesMapping 
        (Map<Integer, Integer> map, BitSet mapped) {

        List<ANode> nodes = new ArrayList<ANode>();
        BitSet set1 = new BitSet();
        BitSet set2 = new BitSet();

        for (int i = 0; i < qBonds.length; ++i) {
            MolBond qb = qBonds[i];

            int u1 = query.indexOf(qb.getAtom1());
            int v1 = query.indexOf(qb.getAtom2());
            boolean bu1 = map.containsKey(u1);
            boolean bv1 = map.containsKey(v1);

            if (bu1 && bv1) {
                set1.set(i);
            }

            for (int j = 0; j < tBonds.length; ++j) {
                MolBond tb = tBonds[j];

                int u2 = target.indexOf(tb.getAtom1());
                int v2 = target.indexOf(tb.getAtom2());
                boolean bu2 = map.containsValue(u2);
                boolean bv2 = map.containsValue(v2);
                if (bu2 && bv2) {
                    set2.set(j);
                }

                if (bondComparator.match(qb, tb)
                    && (((bu1 || bv1) && (bu2 || bv2))
                        || (!bu1 && !bv1 && !bu2 && !bv2))) {
                    MatchDirection dir = getMatchDirection (qb, tb);
                    if (dir != MatchDirection.NONE) {
                        ANode n = new ANode (i, j, u1, v1, u2, v2, dir);
                        Integer u = map.get(u1);
                        Integer v = map.get(v1);

                        if (u != null && v != null) {
                            if (u == u2 && v == v2) {
                                n.dir = MatchDirection.FORWARD;
                            }
                            else if (u == v2 && v == u2) {
                                n.dir = MatchDirection.REVERSE;
                            }
                            else {
                                n = null;
                            }

                            if (n != null) {
                                mapped.set(nodes.size());
                                nodes.add(n);
                            }
                        }
                        else if (u != null && v == null) {
                            if (u == u2 && !bv2) {
                                n.dir = MatchDirection.FORWARD;
                            }
                            else if (u == v2 && !bu2) {
                                n.dir = MatchDirection.REVERSE;
                            }
                            else {
                                n = null;
                            }

                            if (n != null) {
                                /* Note that this is not considered as a
                                 * mapped node, since only half of it is
                                 * mapped
                                 */
                                nodes.add(n);
                            }
                        }
                        else if (u == null && v != null) {
                            if (v == u2 && !bv2) {
                                n.dir = MatchDirection.REVERSE;
                            }
                            else if (v == v2 && !bu2) {
                                n.dir = MatchDirection.FORWARD;
                            }
                            else {
                                n = null;
                            }

                            if (n != null) {
                                // similarly, this is not a mapped node
                                nodes.add(n);
                            }
                        }
                        else {
                            nodes.add(n);
                        }
                    }
                }
            }
        } // endfor each query bond

        // now create association nodes for the mapping
        for (int i = set1.nextSetBit(0); i >= 0;
             i = set1.nextSetBit(i+1)) {
            MolBond qb = qBonds[i];
            int u1 = query.indexOf(qb.getAtom1());
            int v1 = query.indexOf(qb.getAtom2());

            int u2 = map.get(u1), v2 = map.get(v1);
            int j = tBtab[u2][v2];
            if (j < 0 || !set2.get(j)) {
                logger.warning("Bogus mapping found "+(u1+1)+":"+(u2+1)+" "
                               +(v1+1)+":"+(v2+1)
                               +"; MCS might not be correct!");
            }
        }

        if (DEBUG > 0) {
            logger.info("Mapping vertices: " + mapped.cardinality());
            System.err.print("$$");
            for (int i = mapped.nextSetBit(0);
                 i >= 0; i = mapped.nextSetBit(i+1)) {
                System.err.print(" "+i+nodes.get(i));
            }
            System.err.println();
            System.err.println("** " + nodes.size() + " Vertices");
            int dir = 0;
            for (int i = 0; i < nodes.size(); ++i) {
                ANode n = nodes.get(i);
                System.err.print(String.format("[%1$3d]: ",i) + n);
                if (mapped.get(i)) {
                    System.err.print(" **");
                }
                System.err.println();
                if (n.dir != MatchDirection.BIDIR) {
                    ++dir;
                }
            }
            System.err.println("## "+dir+" direct mapping(s)!");
        }

        return nodes.toArray(new ANode[0]);
    }

    protected AGraph createAGraphFromSeed (Molecule seed) {
        AGraph G = null;

        /*
         * generate all possible isomorphisms from query
         */
        VFLib2 vf1 = VFLib2.subgraphIsomorphism
            (seed, getQuery (), getAtomComparator (), getBondComparator ());
        int[][] qhits = vf1.findAll(false); // all isomorphisms

        /*
         * only unique ones from target
         */
        VFLib2 vf2 = VFLib2.subgraphIsomorphism
            (seed, getTarget (), getAtomComparator (), getBondComparator ());
        int[][] thits = vf2.findAll(true); // only unique isomorphisms

        int k = 0, size = getQuery().getAtomCount();
        for (int[] qh : qhits) {
            for (int[] th : thits) {
                Map<Integer, Integer> map = new HashMap<Integer, Integer>();
                if (DEBUG > 0) {
                    System.err.print("## Seed mapping "+k+":");
                }

                for (int i = 0; i < qh.length; ++i) {
                    if (DEBUG > 0) {
                        System.err.print(" "+(qh[i]+1)+":"+(th[i]+1));
                    }
                    map.put(qh[i], th[i]);
                }

                if (DEBUG > 0) {
                    System.err.println();
                }

                BitSet bs = new BitSet ();
                ANode[] nodes = createANodesMapping (map, bs);

                AGraph ag = new AGraph (nodes, map, bs);
                /*
                 * TODO: this is not complete since we're simply selecting
                 * the best seed here. This means that if allcliques is set,
                 * only the non-seed mappings will be enumerated.
                 */
                if (G == null || G.cEdges() < ag.cEdges()) {
                    G = ag;
                }
                ++k;
            }
        }

        return G;
    }

    private boolean initialize () {
        initVars ();

        if (seed != null) {
            G = createAGraphFromSeed (seed);
            // have the mapping reflective of the actual mappings generated
            //  by the seed
            mapping.clear(); 
            mapping.putAll(G.mapping);
        }
        else {
            ANode[] nodes;
            BitSet mapped = new BitSet ();
            if (!mapping.isEmpty()) {
                if (!validateMapping ()) {
                    return false;
                }

                nodes = createANodesMapping (mapping, mapped);
            }
            else {
                nodes = createANodesNoMapping ();
            }
            G = new AGraph (nodes, mapping, mapped);
        }

        return true;
    }

    /*
     * initialize variables
     */
    private void initVars () {
        start = System.currentTimeMillis();
        stop = 0;
        ctime = 0;
        cdepth = 0;
        tdepth = 0;

        qBonds = query.getBondArray();
        tBonds = target.getBondArray();

        forward = new int[qAtoms.length];
        reverse = new int[tAtoms.length];

        G = null;
        maxcliques.clear();

        resetVars ();
    }

    void resetVars () {
        hits.clear();
        maxcore = null;
        maxscore = Integer.MIN_VALUE;
        maxcsize = 0;
        maxclique.clear();
    }

    protected boolean validateMapping () {
        if (mapping.isEmpty()) {
            return true;
        }

        List<Integer> remove = new ArrayList<Integer>();
        // check to make sure the mapping are within bounds
        for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
            Integer q = e.getKey();
            Integer t = e.getValue();
            if (q >= 0 && q < qAtoms.length
                && t >= 0 && t < tAtoms.length) {
                if (!atomComparator.match(qAtoms[q], tAtoms[t])) {
                    logger.warning("Mapping "+(q+1)+" => "+(t+1)
                                   +" doesn't match atom comparator; "
                                   +"mapping ignored!");
                    System.err.print(" Query");
                    _atom(System.err, q+1, qAtoms[q]);
                    System.err.print("Target");
                    _atom(System.err, t+1, tAtoms[t]);

                    remove.add(q);
                }
            }
            else {
                logger.warning
                    ("Remove bogus mapping: "+(q+1)+" => "+(t+1));
                remove.add(q);
            }
        }

        if (!remove.isEmpty()) {
            System.err.println("** Total of "+remove.size()
                               + " bogus mapping(s) removed!");
            System.err.println
                ("**  QUERY: "+getQuery().toFormat("smiles:q")
                 +"\t"+getQuery().getName());
            System.err.println
                ("** TARGET: "+getTarget().toFormat("smiles:q")
                 +"\t"+getTarget().getName());
            System.err.print("** MAPPING:");
            for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
                System.err.print(" "+e.getKey()+":"+e.getValue());
            }
            System.err.println();
        }

        for (Integer key : remove) {
            mapping.remove(key);
        }

        return remove.isEmpty();
    }

    /*
     * approx coloring routine based on Tomita & Kameda's paper:
     * An efficient branch-and-bound algorithm for finding a maximum clique
     *  with computational experiments
     * doi: 10.1007/s10898-006-9039-7
     */
    protected Map<Integer, Integer> coloring (BitSet S) {
        int maxno = 0;

        BitSet[] K = new BitSet[S.cardinality()];
        for (int v = 0; v < K.length; ++v) {
            K[v] = new BitSet (G.size());
        }
        final Map<Integer, Integer> color = new TreeMap<Integer, Integer>();

        for (int v = S.nextSetBit(0); v >= 0; v = S.nextSetBit(v+1)) {
            int k = 0;
            while (K[k].intersects(G.cEdges(v))) {
                ++k;
            }

            if (k > maxno) {
                maxno = k;
                K[k].clear();
            }
            color.put(v, k+1);
            K[k].set(v);
        }

        return color;
    }

    protected static void _atom (java.io.PrintStream ps, int index, MolAtom a) {
        ps.println("  " + String.format("%1$3d", index)
                   +"[a="+a.getAtno()
                   +",m="+a.getAtomMap()
                   +",c="+a.getCharge()
                   +",h="+a.getImplicitHcount()
                   +",H="+a.getExplicitHcount()
                   +",r="+a.getRadical()
                   +",q="+a.getQuerystr()
                   +",l="+a.getQueryLabel()
                   +",x="+a.getExtraLabel()
                   +",v="+a.getValence()
                   +",s="+a.getSymbol()
                   +",b="+a.hasAromaticBond()
                   +"]");
    }

    // Bron-Kerbosch algorithm
    private void bronKerbosch (BitSet C, BitSet P, BitSet S) {
        if (P.isEmpty() && S.isEmpty()) {
            outputClique (C, null);
        }
        else {
            for (int u = P.nextSetBit(0); u >=0 ; u = P.nextSetBit(u+1)) {
                P.clear(u);
                BitSet PP = (BitSet)P.clone();
                BitSet SS = (BitSet)S.clone();
                PP.and(G.edges(u));
                SS.and(G.edges(u));
                C.set(u);
                bronKerbosch (C, PP, SS);
                C.clear(u);
                S.set(u);
            }
        }
    }

    // modified Bron-Kerbosh algorithm
    private void bronKerbosh2 (BitSet C, BitSet P, BitSet S) {
        if (P.isEmpty()) {
            if (S.isEmpty()) {
                outputClique (C, null);
            }
        }
        else {
            int v = pickVertex (P);
            for (int u = P.nextSetBit(0); u >=0 ; u = P.nextSetBit(u+1)) {
                if (!G.isNeighbor(u, v)) {
                    P.clear(u);
                    BitSet PP = (BitSet)P.clone();
                    BitSet SS = (BitSet)S.clone();
                    PP.and(G.edges(u));
                    SS.and(G.edges(u));
                    C.set(u);
                    bronKerbosh2 (C, PP, SS);
                    C.clear(u);
                    S.set(u);
                }
            }
        }
    }

    private boolean genMaxCliques (BitSet P, BitSet R, BitSet clique) {
        boolean done = false;
        if (P.isEmpty()) {
            done = outputClique (R, null);
        }
        else {
            while (!P.isEmpty() && !done && extendable (P, R, clique)) {
                int v = pickMinVertex (P);

                R.set(v);
                P.clear(v);
                BitSet PP = (BitSet)P.clone();
                PP.and(G.edges(v));

                done = genMaxCliques (PP, R, clique);
                R.clear(v);
            }
        }
        return done;
    }

    /*
     * FIXME: this is a very weak bound.... use graph coloring here?
     */
    private boolean extendable (BitSet P, BitSet C, BitSet clique) {
        return ((C.cardinality() + P.cardinality())
                >= clique.cardinality()) && checkTime ();
    }

    private boolean checkTime () {
        if (maxtime > 0) {
            long time = System.currentTimeMillis();
            if ((time - start) > maxtime) {
                if (stop == 0) {
                    logger.info("** MaxClique search truncated");
                    stop = time;
                }
                return false;
            }
        }
        return true;
    }

    // assume initialize is called first!
    private void genMaxCliques () {
        BitSet P = new BitSet (G.size());
        for (int i = 0; i < G.size(); ++i) {
            P.set(i);
        }
        BitSet C = new BitSet (G.size());

        //System.out.println("\n** MaxClique");
        genMaxCliques (P, C, maxclique);
    }

    private void genCliques () {
        BitSet P = new BitSet (G.size());
        for (int i = 0; i < G.size(); ++i) {
            P.set(i);
        }
        BitSet C = new BitSet (G.size());
        BitSet S = new BitSet (G.size());

        //System.out.println("\n** AllCliques");
        bronKerbosh2 (C, P, S);
    }

    static void dfs (int node, BitSet[] g, BitSet seen) {
        seen.set(node);
        BitSet nbr = g[node];

        for (int u = nbr.nextSetBit(0); u >= 0; u = nbr.nextSetBit(u+1)) {
            if (!seen.get(u)) {
                dfs (u, g, seen);
            }
        }
    }

    /*
     * Main entry point for generating c-cliques
     */
    private void genC_Cliques () {
        BitSet[] components = G.components();

        for (int i = 0; i < components.length; ++i) {
            BitSet A = components[i];

            if (DEBUG > 0) {
                System.err.println("[**** component " + (i+1) + " *****]");
                System.err.println("%% " + A.cardinality() + A);
            }

            if (A.cardinality() < maxclique.cardinality()) {
                if (DEBUG > 0) {
                    System.err.println
                        ("## component "+(i+1)
                         +" is less than max clique "+maxclique.cardinality()
                         +"; skipping this component!");
                }
                break;
            }

            if (G.mapped().isEmpty()) {
                if (nthreads <= 1) {
                    genC_Cliques (A);
                }
                else {
                    genC_CliquesThreaded (A, nthreads);
                }
            }
            else if (A.intersects(G.mapped())) {
                genC_CliquesMapped (A);
            }
            else if (DEBUG > 0) {
                System.err.println
                    ("## skipping component "+(i+1)
                     +" since it doesn't contain mapping vertices!");
            }

            if (maxcore != null) {
                MClique clique = new MClique (G, A);
                clique.core = maxcore;
                clique.size = maxcsize;
                clique.cliques.putAll(hits);
                clique.elapsed = searchTimeMillis ();
                //System.out.println(clique);
                maxcliques.add(clique);
            }

            // now reset 
            resetVars ();
        }
    }

    private void genC_CliquesThreaded (BitSet A, int n) {
        final BitSet T = new BitSet (G.size());

        ExecutorService threadPool = Executors.newFixedThreadPool(n);
        List<Future> tasks = new ArrayList<Future>();

        final AtomicBoolean done = new AtomicBoolean (false);
        for (int u = A.nextSetBit(0);
             u >= 0 && !done.get(); u = A.nextSetBit(u+1)) {

            final BitSet P = new BitSet (G.size());
            final BitSet Q = new BitSet (G.size());
            final BitSet X = new BitSet (G.size());
            final BitSet Y = new BitSet (G.size());
            final BitSet R = new BitSet (G.size());

            Q.or(A);
            Q.andNot(T);
            Q.and(G.dEdges(u)); // Q = (V - T) + D[u]
            Y.or(G.dEdges(u));
            Y.and(T);

            P.or(A);
            P.andNot(T);
            P.and(G.cEdges(u)); // P = (V - T) + C[u]
            X.or(G.cEdges(u));
            X.and(T);

            if (DEBUG > 2) {
                System.err.println("^^ Node " + u + " "
                                   + G.edges(u) + " P = "
                                   + P.cardinality() + P
                                   + " Q = " + Q.cardinality() + Q);
            }

            R.set(u);
            T.set(u);

            final int v = u;
            Future task = threadPool.submit
                (new Runnable () {
                        public void run () {
                            LinkedList<Edge> path = new LinkedList<Edge>();
                            done.set(genC_Cliques (R, P, Q, X, Y, v, path));
                            if (DEBUG > 0) {
                                logger.info(Thread.currentThread()
                                            +": Search Vertex = "+v+"; C="
                                            +G.cEdges(v).cardinality()
                                            +" D="+G.dEdges(v).cardinality()
                                            +" "+(System.currentTimeMillis()-start)
                                            +"ms");
                            }
                        }
                    });
            tasks.add(task);
        }

        for (Future f : tasks) {
            if (!f.isDone()) {
                try {
                    f.get(); // wait for thread to finish
                }
                catch (Exception ex) {
                    logger.log(Level.WARNING, "Thread "+f+" interrupted!", ex);
                }
            }
        }
        threadPool.shutdown();
    }

    private void genC_CliquesMapped (BitSet A) {
        BitSet P = new BitSet (G.size());
        BitSet Q = new BitSet (G.size());
        BitSet X = new BitSet (G.size());
        BitSet Y = new BitSet (G.size());
        BitSet R = new BitSet (G.size());

        BitSet Cu = new BitSet (G.size());
        BitSet Du = new BitSet (G.size());
        BitSet mapped = G.mapped();
        for (int u = mapped.nextSetBit(0);
             u >= 0; u = mapped.nextSetBit(u+1)) {
            Cu.or(G.cEdges(u));
            Du.or(G.dEdges(u));
        }

        Cu.andNot(mapped);
        Du.andNot(mapped);

        R.or(mapped);

        Q.or(A);
        Q.and(Du);

        P.or(A);
        P.and(Cu);

        if (DEBUG > 0) {
            System.err.print("## R = ");
            for (int i = R.nextSetBit(0); i >= 0; i = R.nextSetBit(i+1)) {
                System.err.print(" " + i + G.vertex(i));
            }
            System.err.println();

            System.err.print("## P = ");
            for (int i = P.nextSetBit(0); i >= 0; i = P.nextSetBit(i+1)) {
                System.err.print(" " + i + G.vertex(i));
            }
            System.err.println();

            System.err.print("## Q = ");
            for (int i = Q.nextSetBit(0); i >= 0; i = Q.nextSetBit(i+1)) {
                System.err.print(" " + i + G.vertex(i));
            }
            System.err.println();
        }

        LinkedList<Edge> path = new LinkedList<Edge>();
        genC_Cliques (R, P, Q, X, Y, -1, path);
    }


    private void genC_Cliques (BitSet A) {
        final BitSet T = new BitSet (G.size());
        for (int u = A.nextSetBit(0); u >= 0; u = A.nextSetBit(u+1)) {

            final BitSet P = new BitSet (G.size());
            final BitSet Q = new BitSet (G.size());
            final BitSet X = new BitSet (G.size());
            final BitSet Y = new BitSet (G.size());
            final BitSet R = new BitSet (G.size());

            Q.or(A);
            Q.andNot(T);
            Q.and(G.dEdges(u)); // Q = (V - T) + D[u]
            Y.or(G.dEdges(u));
            Y.and(T);

            P.or(A);
            P.andNot(T);
            P.and(G.cEdges(u)); // P = (V - T) + C[u]
            X.or(G.cEdges(u));
            X.and(T);

            if (DEBUG > 2) {
                System.err.println("^^ Node " + u + " "
                                   + G.vertex(u) + " P = "
                                   + P.cardinality() + P
                                   + " Q = " + Q.cardinality() + Q);
            }

            R.set(u);
            T.set(u);

            LinkedList<Edge> path = new LinkedList<Edge>();
            if (genC_Cliques (R, P, Q, X, Y, u, path)) {
                break;
            }

            if (DEBUG > 0) {
                logger.info(Thread.currentThread()
                            +": Search Vertex = "+u+"; C="
                            +G.cEdges(u).cardinality()
                            +" D="+G.dEdges(u).cardinality()
                            +" "+(System.currentTimeMillis()-start)
                            +"ms");
            }
        }
    }

    Collection<Integer> orderSet (BitSet S) {
        List<Integer> order = new ArrayList<Integer>();
        for (int i = S.nextSetBit(0); i >= 0; i = S.nextSetBit(i+1)) {
            order.add(i);
        }
        Collections.sort
            (order, new Comparator<Integer>() {
                public int compare (Integer i, Integer j) {
                    return G.cEdges(j).cardinality()
                        - G.cEdges(i).cardinality();
                }
            });
        return order;
    }

    // Bron-Kerbosch algorithm for connected cliques
    private boolean genC_Cliques (BitSet R, BitSet P, 
                                  BitSet Q, BitSet X, 
                                  BitSet Y, int v, 
                                  LinkedList<Edge> path) {
        boolean done = false;
        if (P.isEmpty() && X.isEmpty()) {
            done = outputClique (R, path);
        }
        else if ((allcliques || evalByScore
                  // very loose upper bound for max clique
                  || ((R.cardinality() + P.cardinality()
                       + Q.cardinality())
                      >= Math.max(mapping.size(), maxclique.cardinality())))
                 && checkTime() // time constraint
                 // depth constraint
                 && (maxDepth == 0 || tdepth < maxDepth)) {

            if (DEBUG > 3) {
                System.err.println(">> clique = " + R.cardinality() + R
                                   + " P = "+P.cardinality() + P
                                   + " Q = "+Q.cardinality() + Q);
            }

            for (int u = P.nextSetBit(0);
                 u >= 0 && !done; u = P.nextSetBit(u+1)) {

                P.clear(u);
                BitSet PP = (BitSet)P.clone();
                BitSet QQ = (BitSet)Q.clone();
                BitSet XX = (BitSet)X.clone();
                BitSet YY = (BitSet)Y.clone();

                QQ.and(G.dEdges(u));
                YY.and(G.dEdges(u));
                PP.and(G.edges(u));
                XX.and(G.edges(u));
                {
                    BitSet q = (BitSet)Q.clone();
                    q.and(G.cEdges(u));
                    PP.or(q);
                    BitSet x = (BitSet)Y.clone();
                    x.and(G.cEdges(u));
                    XX.or(x);
                }

                R.set(u);
                path.push(new Edge (v, u));

                done = genC_Cliques (R, PP, QQ, XX, YY, u, path);

                Edge e = path.pop();
                R.clear(u);
                X.set(u);
            }
            if (DEBUG > 3) {
                System.err.println("<< clique = " + R.cardinality() + R
                                   + " P = "+P.cardinality() + P
                                   + " Q = "+Q.cardinality() + Q);
            }
        }

        return done;
    }

    void setupMapping (int[] fwd, int[] rev) {
        /*
          if (!mapping.isEmpty()) {
          for (int i = 0; i < rev.length; ++i) {
          rev[i] = -1;
          }
          for (int i = 0; i < fwd.length; ++i) {
          Integer m = mapping.get(i);
          if (m != null) {
          fwd[i] = m;
          rev[m] = i;
          }
          else {
          fwd[i] = -1;
          }
          }
          }
          else*/ {
              for (int i = 0; i < rev.length; ++i) {
                  rev[i] = -1;
              }
              for (int i = 0; i < fwd.length; ++i) {
                  fwd[i] = -1;
              }
          }
    }

    /*
     * this assumes that the clique is connected!!!!
     */
    private boolean generateMapping (BitSet C, int[] fwd, int[] rev) {
        setupMapping (fwd, rev);

        /*
         * first consider unambiguous matches first
         */
        LinkedList<ANode> queue = new LinkedList<ANode>();

        Map<ANode, Integer> amap = new HashMap<ANode, Integer>();
        for (int i = C.nextSetBit(0); i >= 0; i = C.nextSetBit(i+1)) {
            ANode n = G.vertex(i);
            amap.put(n, i);

            switch (n.dir) {
            case FORWARD: // consider forward & reverse first
                if (fwd[n.u1] < 0 && rev[n.u2] < 0) {
                    fwd[n.u1] = n.u2;
                    rev[n.u2] = n.u1;
                }
                if (fwd[n.v1] < 0 && rev[n.v2] < 0) {
                    fwd[n.v1] = n.v2;
                    rev[n.v2] = n.v1;
                }
                break;

            case REVERSE:
                if (fwd[n.u1] < 0 && rev[n.v2] < 0) {
                    fwd[n.u1] = n.v2;
                    rev[n.v2] = n.u1;
                }
                if (fwd[n.v1] < 0 && rev[n.u2] < 0) {
                    fwd[n.v1] = n.u2;
                    rev[n.u2] = n.v1;
                }
                break;

            case BIDIR:
                queue.add(n);
                break;

            default:
                throw new IllegalStateException ("Bogus node " + n);
            }
        }

        ANode anchor = null;
        if (queue.size() == C.cardinality()) {
            // if there isn't any anchor matches, we simply pick
            //  the first one as the seed
            ANode n = queue.poll();
            /*if (rev[n.u2] < 0)*/ {
                fwd[n.u1] = n.u2;
                rev[n.u2] = n.u1;
            }
            /*if (rev[n.v2] < 0)*/ {
                fwd[n.v1] = n.v2;
                rev[n.v2] = n.v1;
            }
            anchor = n;
        }

        // now do a second pass to update the rest of the matches
        LinkedList<ANode> aux = new LinkedList<ANode>();
        int[] maps = new int[qBonds.length];
        int nflips = 0, niters = 0;

        // some cutoff in case we have bogus mapping that results
        //  in infinite loop
        int maxiters = G.cEdges() * 4;

        while (!queue.isEmpty()) {
            ANode n = queue.poll();

            //System.out.println(">> polling " + n);
            boolean flip = false;

            if (fwd[n.u1] < 0 && fwd[n.v1] >= 0) {
                if (fwd[n.v1] == n.v2) { // forward
                    if (rev[n.u2] < 0) {
                        fwd[n.u1] = n.u2;
                        rev[n.u2] = n.u1;
                    }
                }
                else if (fwd[n.v1] == n.u2) { // reverse
                    if (rev[n.v2] < 0) {
                        fwd[n.u1] = n.v2;
                        rev[n.v2] = n.u1;
                    }
                }
                else if (anchor != null) { // reverse the anchor
                    flip = true;
                }
            }
            else if (fwd[n.u1] >= 0 && fwd[n.v1] < 0) {
                if (fwd[n.u1] == n.u2) { // forward
                    if (rev[n.v2] < 0) {
                        fwd[n.v1] = n.v2;
                        rev[n.v2] = n.v1;
                    }
                }
                else if (fwd[n.u1] == n.v2) { // reverse
                    if (rev[n.u2] < 0) {
                        fwd[n.v1] = n.u2;
                        rev[n.u2] = n.v1;
                    }
                }
                else if (anchor != null) {
                    flip = true;
                }
            }
            else if (fwd[n.u1] < 0 && fwd[n.v1] < 0) {
                //System.out.println("<< putting " + n + " back");
                /*
                  ++maps[n.n1];
                  if (maps[n.n1] > 1) {
                  // look like we have disconnected components
                  anchor = n;
                  fwd[n.u1] = n.u2;
                  rev[n.u2] = n.u1;
                  fwd[n.v1] = n.v2;
                  rev[n.v2] = n.v1;
                  maps[n.n1] = 0;
                  }
                  else*/ {
                      queue.add(n);
                      n = null;
                  }
            }

            if (flip) {
                if (DEBUG > 0) {
                    System.err.println("## flipping anchor "
                                       + anchor + " due to node " + n);
                    System.err.print("##");
                    for (int i = 0; i < fwd.length; ++i) {
                        if (fwd[i]>= 0) {
                            System.err.print
                                (" " + (i+1) + "-" + (fwd[i]+1));
                        }
                    }
                    System.err.println();
                }

                if (++nflips == C.cardinality()) {
                    logger.warning("Unable to generate proper "
                                   +"mapping after "+nflips
                                   +" flips; offending node "+n);
                    for (ANode a : queue) {
                        int u = amap.get(a);
                        C.clear(u); // remove this bogus node
                    }
                    return false;
                }

                setupMapping (fwd, rev);
                niters = 0;

                // reverse the anchor
                fwd[anchor.u1] = anchor.v2;
                rev[anchor.v2] = anchor.u1;
                fwd[anchor.v1] = anchor.u2;
                rev[anchor.u2] = anchor.v1;

                queue.add(n);
                queue.addAll(aux);
                aux.clear();
                n = null;
            }

            if (n != null) {
                aux.add(n);
            }

            if (false && ++niters > maxiters) {
                logger.warning
                    ("Number of iterations exceeds max ("+maxiters+"); "
                     +"aborting with "+queue.size() + " dangling node(s)!");
                System.err.println("** clique: "+C);
                for (ANode a : queue) {
                    int u = amap.get(a);
                    C.clear(u); // remove this bogus node
                    System.err.println
                        ("** "+u+""+a + " C = "
                         +G.cEdges(u) + " D = "+G.dEdges(u));
                }
                if (!mapping.isEmpty()) {
                    System.err.print("** mapping:");
                    for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
                        System.err.print(" "+e.getKey()+":"+e.getValue());
                    }
                    System.err.println();
                }
                System.err.println("**  query: "+getQuery().toFormat("smiles:q"));
                System.err.println("** target: "+getTarget().toFormat("smiles:q"));
                System.err.println();

                return false;
            }
        } // foreach queue

        return true;
    }

    private int getCliqueSize (BitSet C, BitSet qatoms, BitSet qbonds) {
        BitSet t = new BitSet ();

        if (qbonds == null) {
            qbonds = new BitSet ();
        }

        // The size of the clique is not quite the cardinality of C, because
        //   there might be matches where the same bond is matched multiple
        //   times.  So the true clique size is the number of unique bonds
        //   in either direction.
        for (int i = C.nextSetBit(0); i >= 0; i = C.nextSetBit(i+1)) {
            ANode n = G.vertex(i);
            qbonds.set(n.n1); // query bonds
            t.set(n.n2); // target bonds

            if (qatoms != null) {
                qatoms.set(n.u1);
                qatoms.set(n.v1);
            }
        }

        return Math.min(qbonds.cardinality(), t.cardinality());
    }

    protected boolean evaluateCore (Molecule core) {
        boolean ok = false;

        if (maxcore == null
            || (maxcore.getBondCount() < core.getBondCount())) {
            ok = true;
        }

        return ok;
    }

    boolean checkBonds (Molecule m, int[][]  btab, int[] map) {
        for (MolBond b : m.getBondArray()) {
            int a1 = m.indexOf(b.getAtom1()),
                a2 = m.indexOf(b.getAtom2());
            int i1 = map[a1], i2 = map[a2];
            if (i1 >= 0 && i2 >= 0) {
                if (btab[i1][i2] < 0) {
                    // make sure the correspond mapping
                    if (DEBUG > 0) {
                        System.err.println
                            ("@@@@ atoms "
                             +(a1+1)+" and "+(a2+1)+" are "
                             +"connected but their mapping aren't!");
                    }
                    return false;
                }
            }
        }
        return true;
    }

    private synchronized 
        boolean outputClique (BitSet C, LinkedList<Edge> path) {
            
            long time = System.currentTimeMillis();
            tdepth += 1;

            BitSet qatoms = new BitSet ();
            BitSet qbonds = new BitSet ();
            int cliqueSize = getCliqueSize (C, qatoms, qbonds);
            
            if (DEBUG > 0) {
                System.err.println(">> clique = " + C.cardinality() 
                                   + " ("+cliqueSize+"/"+maxcsize+") " + C);
            }

            if (evalByScore) {
            }
            else if (cliqueSize < maxcsize) {
                return false;
            }

            // now check to see if this clique contains the mapped clique
            if (!G.mapped().isEmpty()) {
                BitSet c = (BitSet)C.clone();
                c.and(G.mapped());
                if (c.equals(G.mapped())) {
                    if (DEBUG > 0) {
                        logger.info("Max clique contains mapped clique");
                    }
                }
                else {
                    if (DEBUG > 0) {
                        logger.info("   Max clique: "+C);
                        logger.info("Mapped clique: "+G.mapped());
                    }
                    return false;
                }
            }

            // generate atom mappings
            int[] fwd = new int[qAtoms.length];
            int[] rev = new int[tAtoms.length];

            if (!generateMapping (C, fwd, rev)) {
                qatoms.clear();
                qbonds.clear();

                // recalculate the clique size if there is something gone
                //  wrong in trying to generate the atom mapping.  this can
                //  happen when running in mapped mode.
                cliqueSize = getCliqueSize (C, qatoms, qbonds);
            }

            if (DEBUG > 0) {
                System.err.println
                    ("** max clique found in " + (time-start) + "ms "
                     + cliqueSize + " (" + C.cardinality() + ") " + C);
                System.err.print("##");
                for (int i = C.nextSetBit(0); i >= 0; i = C.nextSetBit(i+1)) {
                    System.err.print(" " + i + G.vertex(i));
                }
                System.err.println();
            }

            // extract mcs core
            //Molecule core = extractFragment (qatoms, qbonds, fwd);
            Molecule core = ChemUtil.createMolecule
                (query, target, fwd, bondComparator, outmap.ord());
            int score = ChemUtil.calcMolScore(core);

            // score mcs core
            if (DEBUG > 0) {
                System.err.print("++");
                for (int i = 0; i < fwd.length; ++i) {
                    if (fwd[i] >= 0) {
                        System.err.print(" " + (i+1) + "-" + (fwd[i]+1));
                    }
                }
                System.err.println();
                System.err.println("=> " + core.toFormat("cxsmarts:q") 
                                   + " " + core.getAtomCount() 
                                   + " " + core.getBondCount()
                                   + " " + score+"/"+maxscore);
            }

            /* TODO: turning this on will not allow matching of partial ring
            if (!checkBonds (query, tBtab, fwd)
                || !checkBonds (target, qBtab, rev))
                return false;
            */
            
            boolean done = false;
            if ((evalByScore 
                 && ((maxscore < score)
                     || (maxscore == score && evaluateCore (core))))
                || (!evalByScore && evaluateCore (core))) {
                maxcore = core;
                maxscore = score;
                maxclique.clear();
                maxclique.or(C);
                maxcsize = cliqueSize;
                if (path != null) {
                    maxcpath = path.toArray(new Edge[0]);
                }

                if (isOutputAtomMapFromMappingOnly ()) {
                    switch (outmap) {
                    case Query:
                        for (MolAtom a : core.getAtomArray()) {
                            int map = a.getAtomMap() - 1;
                            if (map >= 0 && !mapping.containsKey(map)) {
                                // removing this mapping since it's not part of
                                //   the original seed
                                a.setAtomMap(0); 
                            }
                        }
                        break;

                    case Target: 
                        {  Set<Integer> maps = 
                                new HashSet<Integer>(mapping.values());
                            for (MolAtom a : core.getAtomArray()) {
                                int map = a.getAtomMap() - 1;
                                if (map >= 0 && !maps.contains(map)) {
                                    // removing this mapping since it's not part of
                                    //   the original seed
                                    a.setAtomMap(0); 
                                }
                            }
                        }
                        break;
                    }
                }

                ctime = time; // time to the first maximal clique
                cdepth += 1;

                if (DEBUG > 0) {
                    logger.info("Max clique at depth "+cdepth+"/"+tdepth
                                +"; "+(time - start)+"ms");
                }

                System.arraycopy(fwd, 0, forward, 0, fwd.length);
                System.arraycopy(rev, 0, reverse, 0, rev.length);

                BitSet mBonds = new BitSet (qBonds.length);
                mBonds.flip(0, qBonds.length-1);
                mBonds.andNot(qbonds);

                // bail out early if one graph is a subgraph of another
                done = mBonds.isEmpty() && !(allcliques || evalByScore);

                hits.clear();
                hits.put((BitSet)C.clone(), (int[])fwd.clone());

                if (DEBUG > 0) {
                    System.err.print("## unmatched query bonds:");  
                    for (int i = mBonds.nextSetBit(0); 
                         i >= 0; i = mBonds.nextSetBit(i+1)) {
                        System.err.print
                            (" "+(query.indexOf(qBonds[i].getAtom1())+1)
                             +"-" + (query.indexOf(qBonds[i].getAtom2())+1));
                    }
                    System.err.println();
                    if (done) {
                        logger.info(">> search stopped; max mcs found! <<");
                    }
                    else {
                        logger.info(">> maxcore updated <<");
                    }
                }
            }
            else if (maxcsize == cliqueSize 
                     && maxcore.getBondCount() == core.getBondCount()) {
                // all else equal, try to retain the best score core
                if (score > maxscore) { 
                    maxscore = score;
                    maxcore = core;
                }
                hits.put((BitSet)C.clone(), (int[])fwd.clone());
            }

            return done;
        }

    public static void debugHead (OutputStream os) throws IOException {
        PrintStream ps = new PrintStream (os);
        ps.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
        ps.println("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"");
        ps.println(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
        ps.println(" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">");
        ps.println("<key id=\"u1\" for=\"node\" attr.name=\"u1\" attr.type=\"int\"/>");
        ps.println("<key id=\"v1\" for=\"node\" attr.name=\"v1\" attr.type=\"int\"/>");
        ps.println("<key id=\"u2\" for=\"node\" attr.name=\"u2\" attr.type=\"int\"/>");
        ps.println("<key id=\"v2\" for=\"node\" attr.name=\"v2\" attr.type=\"int\"/>");
        ps.println("<key id=\"dir\" for=\"node\" attr.name=\"dir\" attr.type=\"string\"/>");
        ps.println("<key id=\"b1\" for=\"node\" attr.name=\"qbond\" attr.type=\"int\"/>");
        ps.println("<key id=\"b2\" for=\"node\" attr.name=\"tbond\" attr.type=\"int\"/>");
        ps.println("<key id=\"nc\" for=\"node\" attr.name=\"maxclique\" attr.type=\"int\"/>");
        ps.println("<key id=\"ns\" for=\"node\" attr.name=\"seed\" attr.type=\"int\"/>");
        // edge type: 1 = c-edge; 2 = d-edge
        ps.println("<key id=\"type\" for=\"edge\" attr.name=\"type\" attr.type=\"int\"/>");
        ps.println("<key id=\"es\" for=\"edge\" attr.name=\"seed\" attr.type=\"int\"/>");
        ps.println("<key id=\"ec\" for=\"edge\" attr.name=\"maxclique\" attr.type=\"int\"/>");
    }


    public static void debugTail (OutputStream os) throws IOException {
        PrintStream ps = new PrintStream (os);
        ps.println("</graphml>");
        ps.flush();
    }

    public boolean search () {
        if (query == null) {
            throw new IllegalStateException ("No query molecule specified");
        }

        if (target == null) {
            throw new IllegalStateException ("No target molecule specified");
        }

        if (!initialize ()) {
            return false; // construct association graph
        }

        genC_Cliques ();
        Collections.sort(maxcliques);

        long time = System.currentTimeMillis();
        if (stop == 0) {
            stop = time;
        }

        if (DEBUG > 0) {
            System.err.println("** MCSMaxClique total search time: "
                               + String.format("%1$dms",
                                               (stop - start) /* 1e-3*/));
        }

        return !maxcliques.isEmpty();
    }

    public long totalTimeMillis () { return stop - start; }
    public long searchTimeMillis () { return ctime - start; }
    public int searchDepth () { return cdepth; }
    public int totalDepth () { return tdepth; }

    public void setOutputAtomMap (OutputAtomMap outmap) {
        this.outmap = outmap;
    }
    public OutputAtomMap getOutputAtomMap () { return outmap; }

    public int[] getResult () { 
        if (maxcliques.isEmpty()) {
            return new int[0];
        }

        MClique clique = maxcliques.get(0);
        int[][] hits = clique.getHits();
        return hits.length > 0 ? hits[0] : new int[0];
    }

    public Collection<MClique> maxCliques () {
        return Collections.unmodifiableCollection(maxcliques);
    }

    public MClique getMaxClique () {
        if (maxcliques.isEmpty()) return null;
        return maxcliques.get(0);
    }

    public void setMapping (int[] map) {
        if (seed != null) {
            throw new IllegalStateException
                ("Seed is not empty; please clear the seed first!");
        }
        mapping.clear();
        for (int i = 0; i < map.length; ++i) {
            if (map[i] < 0) {
            }
            else {
                mapping.put(i, map[i]);
            }
        }
        // force allcliques to be false since the mapping focuses
        //   the search in a specific direction
        //allcliques = false;
    }
    public void setOutputAtomMapFromMappingOnly (boolean only) {
        atomMapFromMappingOnly = only;
    }
    public boolean isOutputAtomMapFromMappingOnly () { 
        return atomMapFromMappingOnly;
    }

    public void addMapping (int query, int target) {
        if (target < 0) {
            mapping.remove(query);
        }
        else {
            if (seed != null) {
                throw new IllegalStateException
                    ("Seed is not empty; please clear the seed first!");
            }
            mapping.put(query, target);
        }
    }
    public void clearMapping () { mapping.clear(); }

    @Override
    public void setQuery (Molecule query) {
        super.setQuery(query);
        xQuery = null;
    }

    public void setQuery (Molecule query, int[] map) {
        if (map != null && map.length == query.getAtomCount()) {
            xQuery = new BitSet (query.getAtomCount());
            for (int i = 0; i < map.length; ++i)
                if (map[i] >= 0)
                    xQuery.set(i);
        }
        else {
            xQuery = null;
        }
        super.setQuery(query);
    }
    public void setQuery (Molecule query, Integer... xatoms) {
        if (xatoms != null && xatoms.length > 0) {
            xQuery = new BitSet (query.getAtomCount());
            for (int i = 0; i < xatoms.length; ++i)
                xQuery.set(xatoms[i]);
        }
        else {
            xQuery = null;
        }
        super.setQuery(query);
    }

    @Override
    public void setTarget (Molecule target) {
        super.setTarget(target);
        xTarget = null;
    }

    public void setTarget (Molecule target, int[] map) {
        if (map != null && map.length == target.getAtomCount()) {
            xTarget = new BitSet (target.getAtomCount());
            for (int i = 0; i < map.length; ++i)
                if (map[i] >= 0)
                    xTarget.set(i);
        }
        else {
            xTarget = null;
        }
        super.setTarget(target);
    }
    
    public void setTarget (Molecule target, Integer... xatoms) {
        if (xatoms != null && xatoms.length > 0) {
            xTarget = new BitSet (target.getAtomCount());
            for (int i = 0; i < xatoms.length; ++i)
                xTarget.set(xatoms[i]);
        }
        else {
            xTarget = null;
        }
        super.setTarget(target);
    }

    public void setThreads (int nthreads) { 
        if (nthreads < 0) {
            throw new IllegalArgumentException
                ("Bogus number of threads: "+nthreads);
        }
        this.nthreads = nthreads;
    }
    public int getThreads () { return nthreads; }
    public void setEvalByScore (boolean evalByScore) {
        this.evalByScore = evalByScore;
    }
    public boolean getEvalByScore () { return evalByScore; }

    /**
     * Seed the MCS with the specified core; this option is mutually 
     * exclusive with setMapping()!
     */
    public void setSeed (Molecule seed) {
        this.seed = seed;
        if (seed != null) {
            if (!mapping.isEmpty()) {
                throw new IllegalStateException
                    ("Atom mapping is not empty; please clear it first!");
            }
            //this.allcliques = true;
        }
    }

    public void setSeed (String core) {
        try {
            MolHandler mh = new MolHandler (core);
            mh.aromatize();
            setSeed (mh.getMolecule());
        }
        catch (MolFormatException ex) {
            throw new IllegalArgumentException
                ("Bogus molecule format for seed: "+core);
        }
    }

    public MClique maxClique () {
        return maxcliques.isEmpty() ? null : maxcliques.iterator().next();
    }

    public Molecule getSeed () { return seed; }

    public void setAllCliques (boolean allcliques) {
        this.allcliques = allcliques;
    }
    public boolean getAllCliques () { return allcliques; }
    public int[][] getAllHits () { 
        MClique clique = maxClique ();
        return clique != null ? clique.getHits() : new int[0][];
    }

    public int getResultScore () { return maxscore; }
    public int getResultSize () { 
        MClique clique = maxClique ();
        return clique != null ? clique.core.getAtomCount() : 0;
    }
    public Molecule getResultAsMolecule () {
        return getResultAsMolecule (false);
    }
    public Molecule getResultAsMolecule (boolean atomMapping) {
        MClique clique = maxClique ();
        return clique != null ? clique.core() : null;
    }
    public void setTimeoutLimit (long maxtime) { this.maxtime = maxtime; }
    public long getTimeoutLimit () { return maxtime; }


    /****************************************************************
     * TEST CASES
     ***************************************************************/
    public static void test (String[] argv) throws Exception {
        String[][] smi = {
            //{"CC1=C(C(C(=C(N1)C)C(=O)OC)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC", "NCGC00024983-01"},
            //{"CC1=C([C@@H](C(=C(N1)C)C(=O)OCCCN2CCC(CC2)(C3=CC=CC=C3)C4=CC=CC=C4)C=5C=CC=C(C5)[N+](=O)[O-])C(=O)OC", "NCGC00025014-01"},
            //{"CCOC(=O)C1=C(C)NC(C)=C(C1C2=C(Cl)C(Cl)=CC=C2)C(=O)OC","NCGC00093906-01"},
            //{"CCN1C(C)=C(C(C2=CC=CC=C2Cl)C(C([O-])=O)=C1C([O-])=O)C(=O)OC(C)C","NCGC00165748-01"},
            //{"COC(=O)C1=C(C)NC(C)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC(C)(C)CN(C)CCC(C3=CC=CC=C3)C4=CC=CC=C4", "NCGC00167492-01"},
            //{"COC(=O)C1=C(C)NC(C)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OCCN3CCN(CC3)C(C4=CC=CC=C4)C5=CC=CC=C5", "NCGC00167493-01"},
            {"OC(=O)C1=CC2=CC(F)=CC=C2N1", "NCGC00014995-02"},
            {"COC(=O)C1=CC2=C(NC3=C2C=CC=C3)C=N1", "NCGC00093917-01"},
            //{"OCC(=O)(N)", "foo"},
            //{"O=CC(=O)(N)", "bar"},
            //{"C1CC1", "delta"},
            //{"CC(C)C", "Y"},
            //{"[O-][N+](=O)C1=CC=CC=C1", "nitrobenzene"}
        };

        MCSMaxClique mcs = new MCSMaxClique ();
        // mcs.setTimeoutLimit(2000); // truncate search to 2s

        List<Molecule> mols = new ArrayList<Molecule>();
        for (int i = 0; i < smi.length; ++i) {
            Molecule target = new MolHandler (smi[i][0]).getMolecule();
            target.aromatize(Molecule.AROM_BASIC);
            target.hydrogenize(false);
            target.calcHybridization();
            target.setName(smi[i][1]);
            mcs.setTarget(target);
            mols.add(target);
            for (Molecule query : mols) {
                System.out.println
                    ("** " + query.getName() + " vs. " + target.getName());
                mcs.setQuery(query);
                mcs.search();
                System.out.println();
            }
        }
    }


    static Molecule setAtomMapping (Molecule mol) {
        for (int i = 0; i < mol.getAtomCount(); ++i) {
            MolAtom a = mol.getAtom(i);
            a.setAtomMap(i + 1);
        }
        return mol;
    }

    static MCSMaxClique testSeededMCS
        (String seed, String query, String target) throws Exception {

        MCSMaxClique mcs = new MCSMaxClique ();
        //mcs.setAllCliques(true);
        //mcs.setThreads(2);

        MolHandler mh = new MolHandler ();
        mh.setMolecule(seed);
        mh.aromatize();
        mcs.setSeed(mh.getMolecule());

        mh.setMolecule(query);
        mh.aromatize();
        mcs.setQuery(setAtomMapping (mh.getMolecule()));

        mh.setMolecule(target);
        mh.aromatize();
        mcs.setTarget(setAtomMapping (mh.getMolecule()));

        System.out.println("##   Seed: "
                           +mcs.getSeed().toFormat("smiles:q"));
        System.out.println("##  Query: "
                           +mcs.getQuery().getAtomCount()+"/"
                           +mcs.getQuery().getBondCount()+" "
                           +mcs.getQuery().toFormat("smiles:q"));
        System.out.println("## Target: "
                           +mcs.getTarget().getAtomCount()+"/"
                           +mcs.getTarget().getBondCount()+" "
                           +mcs.getTarget().toFormat("smiles:q"));

        if (!mcs.search()) {
            logger.warning("Dude, I got no matches!");
            return mcs;
        }

        System.out.println("## Seeded MCS search took "
                           +String.format("%1$.3fs", 
                                          1e-3*mcs.totalTimeMillis()));
        Molecule core = mcs.getResultAsMolecule();
        System.out.println("##   Core: "+core.getAtomCount()
                           +"/"+core.getBondCount()+" "
                           +core.toFormat("smiles:q"));
        int[][] hits = mcs.getAllHits();
        for (int j = 0; j < hits.length; ++j) {
            int[] hit = hits[j];
            System.out.print ("** Matched " + (j + 1) + ":");
            for (int i = 0; i < hit.length; ++i) {
                if (hit[i] >= 0) {
                    System.out.print (" " + (i + 1) + ":" + (hit[i] + 1));
                }
            }
            System.out.println ();
        }

        return mcs;
    }

    static void testSeededMCSes () throws Exception {
        String[][] tests = new String[][] {
            {"C1=CC=NC=C1", "CP(C)C1=CN=CC(O)=C1", "CCPC1=CC(OC)=CN=C1"},
            {"c1ccccc1", "O(CCCN1CCN(CC1)CCC)c1ccc(cc1)C(=O)CCCCC", "O(CCCN1C[C@H](N(C)C)CC1)c1ccc(cc1)-c1ccc(cc1)C#N"},
            {"c1ccccc1", "COC(=O)C1=C(C)NC(C)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC(C)(C)CN(C)CCC(C3=CC=CC=C3)C4=CC=CC=C4", "COC(=O)C1=C(C)NC(C)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OCCN3CCN(CC3)C(C4=CC=CC=C4)C5=CC=CC=C5"},
            {"C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4",
             "CC(=C)O[C@@]12CO[C@@H]1C[C@H](O)[C@]3(C)C2[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)C(NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](OC(=O)C)C3=O)C5(C)C)C",
             "CC(=O)O[C@H]1C[C@H]2OC[C@@]2(OC(=O)C)[C@H]3[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](O)C(=O)[C@]13C)C5(C)C)C"},
            {"C(CNCC1=CC=CC=C1)COCCC2CCCCCC(C3COC3)C2OCC4=CC=CC=C4",
             "CC(=C)O[C@@]12CO[C@@H]1C[C@H](O)[C@]3(C)C2[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)C(NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](OC(=O)C)C3=O)C5(C)C)C",
             "CC(=O)O[C@H]1C[C@H]2OC[C@@]2(OC(=O)C)[C@H]3[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](O)C(=O)[C@]13C)C5(C)C)C"},
            {"C(CNCCNCC(CC1=CC=CC=C1)NCNC2CCCCC2)NCCNCCC3=CC=CC=C3",
             "CC(C)C[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](Cc2ccccc2)NC(=O)NC34CC5CC(CC(C5)C3)C4)C(=O)N[C@@H](Cc6ccccc6)C(=O)O",
             "CC(C)C[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](Cc2ccccc2)NC(=O)NC34CC5CC(CC(C5)C3)C4)C(=O)N[C@@H](Cc6ccccc6)C(=O)O"}
        };

        for (int i = 0; i < tests.length; ++i) {
            System.out.println
                (">>>>>>>>>>>>>>> TEST "+(i+1)+" <<<<<<<<<<<<<<");
            MCSMaxClique mcs = testSeededMCS
                (tests[i][0], tests[i][1], tests[i][2]);
            Collection<MClique> cliques = mcs.maxCliques();
            int j = 1;
            for (MClique clique : cliques) {
                OutputStream os = new FileOutputStream
                    (new java.io.File 
                     ("maxclique_test_"+(i+1)+"-clique-"+j+".xml"));
                debugHead (os);
                clique.debugG("clique-"+j, os);
                debugTail (os);
                os.close();
                ++j;
            }
        }
    }

    static void testSubgraphIsomorphisms () throws Exception {
        String[][] tests = new String[][]{
            {"C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", "CC(=C)O[C@@]12CO[C@@H]1C[C@H](O)[C@]3(C)C2[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)C(NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](OC(=O)C)C3=O)C5(C)C)C"},
            {"C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", "CC(=O)O[C@H]1C[C@H]2OC[C@@]2(OC(=O)C)[C@H]3[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](O)C(=O)[C@]13C)C5(C)C)C"},
            {"C(C1=CC=CC=C1)C2=CC=CC=C2", "OC(=O)c1ccc(cc1)C2c3ccccc3c4ccccc24"},
        };

        for (int i = 0; i < tests.length; ++i) {
            System.out.println (">> Test " + (i + 1));
            testSubgraphIsomorphism (tests[i][0], tests[i][1]);
        }
    }

    static void testSubgraphIsomorphism (String Q, String T) 
        throws Exception {
        MolHandler mh = new MolHandler (Q);
        mh.aromatize();
        Molecule query = mh.getMolecule();

        mh.setMolecule(T);
        mh.aromatize();
        Molecule target = mh.getMolecule();

        setAtomMapping (target);
        setAtomMapping (query);

        MCSMaxClique mcs = new MCSMaxClique ();
        mcs.setThreads(2);
        mcs.setAllCliques(true);
        mcs.setMolecules(query, target);
        long start = System.currentTimeMillis ();
        mcs.search();
        long end = System.currentTimeMillis ();

        System.out.println ("---- subgraph isomorphism (MaxClique) ----");
        int[][] hits = mcs.getAllHits();
        for (int j = 0; j < hits.length; ++j) {
            int[] hit = hits[j];
            System.out.print ("** Matched " + (j + 1) + ":");
            for (int i = 0; i < hit.length; ++i) {
                System.out.print (" " + (i + 1) + ":" + (hit[i] + 1));
            }
            System.out.println ();
        }

        System.out.println ("## " + hits.length + " hit(s) found in "
                            + String.format ("%1$3dms", end - start));

        chemaxon.sss.search.MolSearch ms = 
            new chemaxon.sss.search.MolSearch ();
        ms.setQuery (query);
        ms.setTarget (target);
        System.out.println ("---- subgraph isomorphism (MolSearch) ----");
        start = System.currentTimeMillis ();
        hits = ms.findAll ();
        end = System.currentTimeMillis ();
        if (hits != null) {
            for (int j = 0; j < hits.length; ++j) {
                int[] hit = hits[j];
                System.out.print ("** Matched " + (j + 1) + ":");
                for (int i = 0; i < hit.length; ++i) {
                    System.out.print (" " + (i + 1) + ":" + (hit[i] + 1));
                }
                System.out.println ();
            }
        }
        System.out.println ("## " + hits.length + " hit(s) found in "
                            + String.format ("%1$3dms", end - start));

        System.out.println (">>  QUERY: " + query.getAtomCount()+"/"
                            +query.getBondCount()+" "
                            +query.toFormat ("smiles:q"));
        System.out.println (">> TARGET: " + target.getAtomCount()+"/"
                            +target.getBondCount()+" "
                            +target.toFormat ("smiles:q"));
    }

    public static void testMCS (String[] argv) throws Exception {
        if (argv.length == 0) {
            System.out.println("usage: MCSMaxClique FILES...");
            System.exit(1);
        }

        List<Molecule> mols = new ArrayList<Molecule>();
        for (int i = 0; i < argv.length; ++i) {
            MolImporter mi = new MolImporter (argv[i]);
            for (Molecule m; (m = mi.read()) != null
                     /*&& mols.size() < 2000*/; ) {
                //m.aromatize(Molecule.AROM_BASIC);
                m.aromatize();
                m.calcHybridization();
                m.hydrogenize(false);
                String name = m.getProperty("field_0");
                if (name != null) {
                    m.setName(name);
                }
                setAtomMapping (m);
                mols.add(m);
            }
            mi.close();
        }
        //System.out.println(mols.size() + " molecules read!");

        MCSMaxClique mcs = new MCSMaxClique ();
        for (int i = 0; i < mols.size(); ++i) {
            Molecule query = mols.get(i);
            mcs.setQuery(query);
            for (int j = i+1; j < mols.size(); ++j) {
                Molecule target = mols.get(j);
                mcs.setTarget(target);
                System.out.println
                    (">>>>>>>>>>>>>>> TEST "+query.getName()
                     +" vs "+target.getName()+" <<<<<<<<<<<<<<");
                if (mcs.search()) {
                    System.out.println("##  Query: "
                                       +query.getAtomCount()+"/"
                                       +query.getBondCount()+" "
                                       +query.toFormat("smiles:q"));
                    System.out.println("## Target: "
                                       +target.getAtomCount()+"/"
                                       +target.getBondCount()+" "
                                       +target.toFormat("smiles:q"));
                    Collection<MClique> cliques = mcs.maxCliques();
                    for (MClique clique : cliques) {
                        System.out.println(clique);
                    }
                }
            }
        }
    }

    public static void main (String[] argv) throws Exception {
        //testSubgraphIsomorphisms ();
        testSeededMCSes ();
        //testMCS (argv);
    }
}

