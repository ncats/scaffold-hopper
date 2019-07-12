// $Id: RGroupSolver.java 3864 2009-12-18 22:09:45Z nguyenda $

package gov.nih.ncgc.rgroup;

import gov.nih.ncgc.algo.graph.MCSMaxClique;
import gov.nih.ncgc.algo.graph.VFLib2;
import gov.nih.ncgc.model.DataSeq;
import gov.nih.ncgc.util.ChemUtil;
import gov.nih.ncgc.util.DefaultAtomComparator;
import gov.nih.ncgc.util.MolFpFactory;
import gov.nih.ncgc.util.SimpleNamer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;


/**
 * Perform R-group decomposition
 */
public class RGroupSolver {
    static private boolean debug = false;
    static private Molecule hydro;
    static private int[] hFP;
    private MolFpFactory fpFact;
    private ConcurrentHashMap<Integer,int[]> fingerprints = new ConcurrentHashMap<Integer,int[]>();
    private DataSeq<Molecule> molseq;
    private boolean DOSETUP=true;
    
    static {
        try {
            debug = Boolean.getBoolean("rgroup.debug");
        }
        catch (Exception ex) {
        }
        try {
            hydro=MolImporter.importMol("*[H]");
        } catch (Exception e) {
            e.printStackTrace();
        }
        MolFpFactory fpFact = MolFpFactory.getInstance();
        hFP=fpFact.generate(hydro);
    }

    static private final Logger logger = Logger.getLogger
        (RGroupSolver.class.getName());


    protected RGroupSolver (DataSeq<Molecule> molseq) {
        this(molseq,true);
    }
    //protected RGroupSolver () {}
    protected RGroupSolver (DataSeq<Molecule> molseq, boolean setupcolumns) {
        this.molseq=molseq;
        fpFact = MolFpFactory.getInstance();
        this.DOSETUP=setupcolumns;
    }
    public static RGroupTable solve (Molecule scaffold, 
                                     DataSeq<Molecule> molseq, 
                                     BitSet subset, boolean forceMatch)
        throws Exception {
        return (new RGroupSolver(molseq)).solve(scaffold,subset,forceMatch);
    }
    public static RGroupTable solve (Molecule scaffold, 
                                     DataSeq<Molecule> molseq, 
                                     BitSet subset) throws Exception {
        return solve(scaffold,molseq,subset,false);
    }
    public RGroupTable solve (Molecule scaffold, 
                              BitSet subset) throws Exception {
        return solve(scaffold,molseq,subset,false);
    }    
    /**
     * This method performs an R-group decomposition of the given
     * scaffold against the sequence of molecules.  Optionally, the
     * subset can be provided so that only a subset of the molecules
     * is used.  This is useful when a given scaffold is known to 
     * be part of only a specific subset of molecules, thereby saving
     * execution cycles.
     */
    public RGroupTable solve (Molecule scaffold,
                              BitSet subset, boolean forceMatch)
        throws Exception {
        boolean doAutomorph=true; //perform automorph to save decomposition time.
        boolean warnNonHit = true; //warn if no hits.
        int numXgroups=0; //number of variable atoms. Based on scaffold.
        
        //First, standardize the scaffold query a little.
        //1) remove all hydrogens
        //2) expand S-groups
        //3a) treat all psuedo atoms as "any", and an X-group
        //3b) remove all terminal RGROUPS, make note of R-group
        //3c) treat all non-terminal RGROUPS as "any", and an X-group
        //3d) treat all lists as lists, but add an X-group
        //4) aromatize
        //5a) if automorph option, determine automorphism of scaffold on self.
        //      (this will eliminate the need to enumerate all maps of scaffold
        //      to target, and re-run the same r-group decomposition on the same 
        //      sets of atoms, and can just permute on automorphs)
        //5b) Otherwise, set molsearch to enumerate all maps, including automorphs
        //6) if forcematch, make all bonds into "ANY" bonds, to deal with 
        //      certain problem tautomers. Will be added back after search. (needs work).
        
        BitSet usedLabels = new BitSet ();
        scaffold.hydrogenize(false);
        scaffold.expandSgroups();
        int[][] assignedGroups = null;
        {
            Set<MolAtom> remove = new HashSet<MolAtom>();
            Map<MolAtom, BitSet> labeledRgroup = 
                new HashMap<MolAtom, BitSet>();
            
            MolAtom[] atoms = scaffold.getAtomArray();
            for (int i = 0; i < atoms.length; ++i) {
                MolAtom a = atoms[i];
                switch (a.getAtno()) {
                case MolAtom.PSEUDO:
                    // force this atom to be any
                    a.setAtno(MolAtom.ANY);
                    numXgroups++;
                    break;
                case MolAtom.RGROUP:
                    if (a.isTerminalAtom()) {
                        // remove this atom... and make a note 
                        if (a.getBondCount() > 0) {
                            MolAtom xa = a.getBond(0).getOtherAtom(a);
                            BitSet bs = labeledRgroup.get(xa);
                            if (bs == null) {
                                labeledRgroup.put(xa, bs = new BitSet ());
                            }
                            bs.set(a.getRgroup());
                            usedLabels.set(a.getRgroup());
                            remove.add(a);
                        }
                    }
                    else { // treat this as a query atom
                        a.setAtno(MolAtom.ANY);
                        numXgroups++;
                    }
                    break;
                    
                case MolAtom.ANY:
                case MolAtom.LIST:
                    numXgroups++;
                    break;
                }
            }
        
            for (MolAtom a : remove) {
                scaffold.removeNode(a);
            }
        
            if (!labeledRgroup.isEmpty()) {
                assignedGroups = new int[scaffold.getAtomCount()][];
                for (Map.Entry<MolAtom, BitSet> e : labeledRgroup.entrySet()) {
                    int pos = scaffold.indexOf(e.getKey());
                    BitSet bs = e.getValue();
                    assignedGroups[pos] = new int[bs.cardinality()];
                    for (int i = bs.nextSetBit(0), j = 0; i >= 0;
                         i = bs.nextSetBit(i+1), ++j) {
                        assignedGroups[pos][j] = i;
                    }
                }
            }
        }
        scaffold.aromatize();
        MolSearch msearch = new MolSearch ();
        int[][] amorph;
        if(doAutomorph){
            amorph=new int[1][scaffold.getAtomCount()];
            for(int i=0;i<amorph[0].length;i++){
                amorph[0][i]=i;
            }
            msearch.setOrderSensitiveSearch(true);
        }else{
            amorph=VFLib2.automorphism(scaffold).findAll();
            msearch.setOrderSensitiveSearch(false);
        }

        Map<MolBond,Integer> oldBond = new HashMap<MolBond,Integer>();
        if (forceMatch) {
            msearch.setOption(MolSearch.OPTION_VAGUE_BOND,
                              MolSearch.VAGUE_BOND_LEVEL4);
            MolBond[] allbonds = scaffold.getBondArray();
            for (MolBond mb : allbonds) {
                oldBond.put(mb, mb.getType());
                mb.setType(MolBond.ANY);
            }
        }
        int[] queryFp = fpFact.generate(scaffold);
        msearch.setQuery(scaffold);
        
        
        // r-group count at each position (there should be at most 3 groups 
        // at a particular atom)
        int[] rgroup = new int[scaffold.getAtomCount()];
        
        if (subset == null) {
            // all bits
            subset = new BitSet (molseq.size());
            subset.flip(0, molseq.size());
            warnNonHit = false;
        }
        int max_attach = 3;
        
        
        HashMap<Integer,Molecule[][]> rtable3 = new HashMap<Integer,Molecule[][]>();
        Map<Integer,Integer> keyTranslate = new HashMap<Integer,Integer>();
        
        //This is a quick fix to bring down the previous size of the rtable.
        //Rather than making an array for every member of molseq, just do it for
        //possible hits. Make a map from molseq position to incremental number.
        for (int i = subset.nextSetBit(0),k=0; 
             i >=0; i = subset.nextSetBit(i+1)) {
            keyTranslate.put(i, k++);
        }
        // Mapping from query to each target 
        int[][] hits = new int[subset.cardinality()][];
        
        if (debug) {
            System.out.println("## scaffold decomposition");
            System.out.print
                (scaffold.toFormat("smarts:u") + "\t" + subset.cardinality());
        }
        
        // sorter for complexity
        Comparator<Molecule> sorter = new Comparator<Molecule>() {
            public int compare (Molecule m1, Molecule m2) {
                return ChemUtil.complexity(m1) - ChemUtil.complexity(m2);
            }
        };
        
        /*
          The idea here is that there are few cases where we don't have 
          exact ligands matching at a particular position, but we do want
          "similar" ligands to align.
          So, we chose the mapping that minimizes the sum of additional bit 
          additions to the set of R-group fingerprints so far.
          In the case of a tie, we chose the mapping that minimizes the number 
          of substituents with new bits.
        */
        int[][] alignments = new int[scaffold.getAtomCount()][];
        /*
          This stores the consensus chirality at each scaffoldatom for each mapping
          to target. If there's disagreement, it is considered unspecified.
        */
        int[] chirality=new int[scaffold.getAtomCount()];
        for(int i=0;i<chirality.length;i++){
            chirality[i]=-1; //unassigned
        }
        
        for (int i = 0; i < alignments.length; ++i) {
            int[] fp = new int[fpFact.getNumInts()];
            for (int j = 0; j < fp.length; ++j) {
                fp[j] = ~0;
            }
            alignments[i] = fp;
        }

        BitSet unmask = new BitSet ();
        for (int i = subset.nextSetBit(0); 
             i >=0; i = subset.nextSetBit(i+1)) {

            Molecule target=null;
            boolean substructure = true;
            // get from cache to avoid
            // regeneration
            // and re-fetching in case
            // of transient molseq.
            int[] targetFp = fingerprints.get(i);                       
            if(targetFp == null){
                target = molseq.get(i);
                target.aromatize();
                targetFp = fpFact.generate(target);
                fingerprints.put(i, targetFp);
            }
            for (int j = 0; j < queryFp.length && substructure; ++j) {
                substructure = (queryFp[j] & targetFp[j]) == queryFp[j];
            }   

            if (!substructure) {
                //subset.clear(i);
                unmask.set(i);
                continue;
            }
            if(target==null){
                target = molseq.get(i);
                target.aromatize();
            }

            msearch.setTarget(target);
            if (debug) {
                System.out.print(" "+ target.getName() + "{");
            }
            try {
                int[][] allHits = msearch.findAll();
                int[] hit = new int[scaffold.getAtomCount()];
                Molecule[] ligands = new Molecule[scaffold.getAtomCount()];
                //get the best mapping of query to hit. If there is one at all 
                if(getBestHitsAndLigands(target, amorph, allHits, alignments, ligands, hit, chirality)){
                    if(debug){System.out.println("Success...");}
                    // String[][] ligtable = new String[scaffold.getAtomCount()][max_attach];
                    Molecule[][] tmpRtable = new Molecule[scaffold.getAtomCount()][max_attach];
                    if(debug){System.out.println("putting rgroupTable ...");}
                    rtable3.put(keyTranslate.get(i),tmpRtable);
                    for (int j = 0; j < hit.length; ++j) {
                        if(debug){System.out.println("Looking at lig ..." + j + " of " + hit.length);}
                        Molecule lig = ligands[j];
                        //if it's not a null ligand, do some indexing
                        if (lig != null) {
                            if (debug) {
                                try{
                                    System.out.print
                                        (","+(j+1)+"-"+(hit[j]+1)+":"
                                         +lig.toFormat("cxsmarts:q"));
                                }catch(Exception e){
                                    System.out.println("prob");
                                }
                            }
                            
                            
                            
                            int nf = lig.getFragCount();
                            if (nf > 1) {
                                //break up multiple fragments into their own r-groups
                                Molecule[] frags = lig.convertToFrags();
                                Arrays.sort(frags, sorter);
                                if (frags.length > max_attach) {
                                    logger.log
                                        (Level.WARNING, target.getName() 
                                         +" generates "+frags.length
                                         +" posible R-groups at atom "
                                         +(hit[j]+1));
                                    nf = max_attach;
                                }       
                                        
                                //store each fragment as seperate in array
                                for (int k = 0; k < nf; ++k) {
                                    tmpRtable[j][k] = frags[k];
                                }
                            }
                            else {
                                //store 
                                tmpRtable[j][0] = lig;
                            }
                            
                            if (assignedGroups != null 
                                && assignedGroups[j] != null) {
                                // preserve use-defined groups
                                if (nf < assignedGroups[j].length) {
                                    nf = assignedGroups[j].length;
                                }
                            }
                            //if there are more attachments defined by user, 
                            //respect that.
                            if (nf > rgroup[j]) {
                                rgroup[j] = nf;
                            }
                        }
                    }
                    hits[keyTranslate.get(i)] = hit;
                    ligands = null;
                }
                else {
                    //subset.clear(i);
                    unmask.set(i);
                    // don't include this molecule in the final rgroup table
                    if (warnNonHit) {
                        logger.log(Level.WARNING,
                                   "scaffold " 
                                   + scaffold.toFormat("cxsmarts:q") 
                                   + " doesn't match " + target.getName() + " "
                                   + target.toFormat("smiles:q"));
                    }
                }
            }
            catch (chemaxon.sss.search.SearchException ex) {
                logger.log(Level.SEVERE, "substructure search error: query=\""
                           +scaffold.toFormat("cxsmarts") + "\" target=\""
                           +target.toFormat("smiles:q") + "\"", ex);
            }
            if (debug) {
                System.out.print("}");
            }
        }
        subset.andNot(unmask);
        
        if(forceMatch){
            for(MolBond mb:oldBond.keySet()){
                mb.setType(oldBond.get(mb));
            }
        }
        scaffold.dearomatize();
        alignments = null;
        int  mincardinality=3;
        if (subset.cardinality() < mincardinality) {
                
            if (false && debug) {
                        
                System.out.println
                    ("... scaffold doesn't have sufficient hit rate!");
                try{
                    System.out.println("Scaffold:"+scaffold.exportToFormat("smiles")+".");
                    for(int i=0;i<molseq.size();i++){
                        System.out.print(molseq.get(i).exportToFormat("smiles")+".");
                    }
                }catch(Exception e){
                    e.printStackTrace();
                }
                System.out.println();
            }
            
            logger.warning("Scaffold cardinality is too small: "+subset);
            return null;
        }

        RGroupTable rgd = new RGroupTable (molseq);
        rgd.rgroupCount = 0;
        { int natt = 0;
            for (int r = 0; r < rgroup.length; ++r) {
                rgd.rgroupCount += rgroup[r];
                if (rgroup[r] > 0) {
                    ++natt;
                }
            }
            
            if (natt >= scaffold.getAtomCount()) {
                // number of attachment positions is at least the number of
                //  atoms, so this scaffold isn't very interesting...
                logger.warning("Number of attachments ("+natt
                               +") >= scaffold size ("
                               +scaffold.getAtomCount()+")!");
                return null;
            }
        }

        int[] labelMap;
        
        // allocate the number of r-groups to be the 
        //  maximum theoretical value
        {       int maxr = scaffold.getAtomCount()*3;
            if(numXgroups>0
               //&&false
               ){
                rgd.xgroupCount=numXgroups;
                rgd.xgroups = new Molecule[subset.cardinality()][numXgroups];
                rgd.xgroupLabels = new String[numXgroups];
                for(int i=0;i<rgd.xgroupCount;i++){
                    rgd.xgroupLabels[i]=((char)('X'+i))+"";
                }
            }
            rgd.rgroups = new Molecule[subset.cardinality()][maxr];
            rgd.rgroupLabels = new String[maxr];
            rgd.ligands = new Map[maxr];
            labelMap = new int[maxr];
        }
        rgd.core = scaffold; // the original input scaffold
        if (scaffold.getDim() < 2) {
            scaffold.clean(2, null);
        }
        rgd.scaffold = scaffold = rgd.core.cloneMolecule();
        rgd.rows = new int[subset.cardinality()];
        rgd.hits = new int[rgd.rows.length][];
        rgd.fp = queryFp;
        rgd.complexity = ChemUtil.complexity(rgd.core);
        
        if (debug) {
            System.out.println();
            System.out.print
                (scaffold.toFormat("cxsmarts:q") + "\t" 
                 + rgd.rgroupCount + "-group at positions");
            for (int i = 0; i < rgroup.length; ++i) {
                if (rgroup[i] > 0) {
                    System.out.print(" " + (i+1)+":"+rgroup[i]);
                }
            }
            System.out.println();
        }
        int lastLabel = 0;
        Map<String,Molecule> cachemap = new HashMap<String,Molecule>();
        
        for (int r = 0, g = 0; r < rgroup.length; ++r) {
            if (rgroup[r] == 0) {
                continue;
            }
        
            int[] bits = new int[fpFact.getNumInts()];
            for (int i = 0; i < bits.length; ++i) {
                bits[i] = ~0;
            }
        
            BitSet rmap = new BitSet ();
            //for each member, get the rgroup table
            //at that column
            for (int i = subset.nextSetBit(0), row = 0; i >= 0; i = subset
                     .nextSetBit(i + 1)) {
                Molecule[] ligands = rtable3.get(keyTranslate.get(i))[r];
                if(debug){
                    System.out.print(i+"\t");
                    System.out.print(r+"\t");
                    for(int m=0;m<ligands.length;m++){
                        try{
                            System.out.print(ligands[m].exportToFormat("cxsmarts:q")+".");
                        }catch(Exception e){
                            //System.out.println(m+"\tnull");
                        }
                    }
                    System.out.println();
                }
            }
            for (int i = subset.nextSetBit(0), row = 0; i >= 0; i = subset
                     .nextSetBit(i + 1)) {
                Molecule[] ligands = rtable3.get(keyTranslate.get(i))[r];
                for (int k = 0; k < rgroup[r]; ++k) {
                    Molecule ligand = ligands[k];

                    if (ligand != null) {

                        if (bits != null) {
                            int[] fp = fpFact.generate(ligand);
                            for (int j = 0; j < fp.length; ++j) {
                                bits[j] &= fp[j];
                            }
                        }

                        int[] rsizes = ligand.getSmallestRingSizeForIdx();
                        MolAtom[] atoms = ligand.getAtomArray();
                        List toRemove = new ArrayList();

                        StringBuilder sb = new StringBuilder();
                        for (int j = 0; j < atoms.length; ++j) {

                            MolAtom a = atoms[j];
                            if (debug) {
                                sb.append(j + ":" + a.getSymbol()
                                          + a.getBondCount() + "_");
                            }
                            if ((r + 1) == a.getAtomMap()) {
                                if (rsizes[j] > 0) {
                                    a.setAttach(MolAtom.ATTACH1);
                                    a.setAtomMap(g + k + 1);
                                    rmap.set(a.getAtomMap());
                                    for (int bi = 0; bi < a.getBondCount(); bi++) {
                                        MolAtom ma = a.getBond(bi)
                                            .getOtherAtom(a);
                                        if (ma.isPseudo()||ma.getAtno()==MolAtom.ANY) {
                                            for (int bi2 = 0; bi2 < ma
                                                     .getBondCount(); bi2++) {
                                                toRemove.add(ma.getBond(bi2));
                                            }
                                            toRemove.add(ma);
                                        }
                                    }
                                } else {
                                    a.setAtomMap(0);
                                    for (int bi = 0; bi < a.getBondCount(); bi++) {
                                        MolAtom ma = a.getBond(bi)
                                            .getOtherAtom(a);
                                        if (ma.isPseudo()||ma.getAtno()==MolAtom.ANY) {
                                            ma.setAtno(MolAtom.PSEUDO);
                                            ma.setAttach(MolAtom.ATTACH1);
                                        }
                                    }
                                }
                            }
                        }
                        if (debug) {
                            sb.append("\n");
                            for (MolBond mb : ligand.getBondArray()) {
                                sb.append(mb.toString());
                            }
                        }
                        for (Object o : toRemove) {
                            if (o instanceof MolBond) {
                                ligand.removeEdge((MolBond) o);
                            } else if (o instanceof MolAtom) {
                                ligand.removeNode((MolAtom) o);
                            }
                        }

                        Map ligs = rgd.ligands[g + k];
                        if (ligs == null) {
                            rgd.ligands[g + k] = ligs = new HashMap<String, Integer>();
                                                        
                        }
                        String lig = null;
                        try {
                            lig = ligand.toFormat("cxsmarts:q");
                        } catch (Exception e) {
                            e.printStackTrace();
                            return null;
                        }
                                                
                        if(cachemap.get(lig)==null){
                            cachemap.put(lig, ligand);
                        }else{
                            ligand=cachemap.get(lig);
                        }
                                                
                        if (debug) {
                            try {
                                System.out.println("adding group " + (g + k)+"\t" + lig);
                                //System.out.println(ligand.toFormat("sdf"));
                            } catch (Exception ex) {

                            }
                        }

                        Integer c = (Integer) ligs.get(lig);
                        ligs.put(lig, c == null ? 1 : (c + 1));//increment lig count.

                        resetEZ(ligand,0);
                        ligand.clean(2, null);
                        resetEZ(ligand,0);
                        rgd.rgroups[row][g + k] = ligand;
                                                
                    } else {
                        // mask all bits out
                        bits = null;
                    }
                }
 
                if (rgd.hits[row] == null) {
                    rgd.hits[row] = hits[keyTranslate.get(i)];
                }

                // mapping to the original table
                rgd.rows[row++] = i;
            }
        
            if (debug) {
                System.out.println((r+1)+": " + rmap);
            }
            for (int k = 0, j = g; k < rgroup[r]; ++k, ++j) {
                boolean newRgroup = true;
            
                //  now check to see if this r-group is degenerate (having
                //  the same ligands for all of its members).  Note that
                //  this is handled for for non-ring ligands.  We assume that
                //  the fragment enumeration will take care of the 
                //  degenerate ring ligands.
                //TODO:fix multiple attach
                
        
                Map ligands = rgd.ligands[j];
                if (ligands != null && ligands.size() == 1 //&& false
                    ) {
                        
                    // don't include this r-group; attach this ligand into
                    //  the scaffold
                    Map.Entry<String, Integer> entry = 
                        (Map.Entry)ligands.entrySet().iterator().next();
                    if (entry.getValue() == rgd.getRowCount()) {
                        if (debug) {
                            logger.info("Pruning R-group: R" 
                                        + (j+1)+" "+entry.getKey());
                        }
                        
                        try {
                            MolHandler mh = new MolHandler (entry.getKey());
                            Molecule ligand = mh.getMolecule();
                            // remove the attachment atoms from the ligand
                            Map<MolAtom, MolAtom> pseudo = new HashMap<MolAtom, MolAtom>();
                            Map<MolAtom, Integer> pseudoBondOrder =     new HashMap<MolAtom, Integer>();
                            int ring=0;
                            for (MolAtom a : ligand.getAtomArray()) {
                                if (a.getAtno() == MolAtom.ANY || 
                                    a.getAtno() == MolAtom.PSEUDO) {
                                    if (a.getBondCount() != 1) {
                                        logger.log
                                            (Level.WARNING, 
                                             "Expecting pseudo atom to "
                                             +"have one bond but instead "
                                             +"found "+a.getBondCount());
                                    }
                                    pseudo.put(a, a.getBond(0).getOtherAtom(a));
                                    pseudoBondOrder.put(a, a.getBond(0).getType());
                                }else{
                                    if(a.getAtomMap()!=0){
                                        ring++;
                                    }
                                }
                            }
                            
                            if (!pseudo.isEmpty()
                                //&&pseudo.size()<=1
                                && ring==0
                                //&& (Math.random()>.5)
                                        
                                ) {
                                
                                for (MolAtom a : pseudo.keySet()) {
                                    ligand.removeNode(a);
                                }
                                
                                scaffold.fuse(ligand);                      
                                for (MolAtom a : pseudo.keySet()) {
                                    // now add the attachment
                                    scaffold.add
                                        (new MolBond (pseudo.get(a), scaffold.getAtom(r),//MolBond.ANY));
                                                      pseudoBondOrder.get(a)));
                                }
                                
                                ligands.clear();
                                newRgroup = false;
                            }
                        }
                        catch (MolFormatException ex) {
                            logger.log(Level.SEVERE, "Bad ligand at scaffold "
                                       +"position "+r, ex);
                        }
                    }
                }
        
                if (newRgroup) {
                    int label = usedLabels.nextClearBit(1);
                    if (assignedGroups != null && assignedGroups[r] != null) {
                        for (int n = 0; n < assignedGroups[r].length; ++n) {
                            if (assignedGroups[r][n] > 0) {
                                label = assignedGroups[r][n];
                                assignedGroups[r][n] = 0;
                                break;
                            }
                        }
                    }
                    else {
                        lastLabel = label;
                    }
                    usedLabels.set(label);
                    labelMap[g] = label;
                    rgd.rgroupLabels[g] = "R"+label;
        
                    if (debug) {
                        System.out.println("## adding new group: " + g 
                                           + " => " + rgd.rgroupLabels[g]
                                           + " used labels " + usedLabels);
                    }
        
                    if (bits != null) { 
                        int c = 0;
                        for (int i = 0; i < bits.length; ++i) {
                            if (bits[i] != 0) {
                                ++c;
                            }
                        }
                        if (debug) {
                            System.out.println
                                ("## possible extension to R"+label+"? " 
                                 + (c>0));
                        }
        
                        /**
                         * it's possible to extend this scaffold at this
                         *   position further
                         */
                        if (c > 0) {
                            rgd.extension.set(g);
                        }
                    }
                    ++g;
                    //if adding rgroup to position, there can be no chirality, because it's meaningless
                    //TODO: handle such cases.
                    chirality[r]=0;
                    addRgroup (scaffold, label, r);
                }
                else {
                    // remove ligands for this r-group
                    for (int i = 0; i < rgd.getRowCount(); ++i) {
                        rgd.rgroups[i][j] = null;
                    }
                    String label = rgd.rgroupLabels[g];
                    if (label != null) {
                        // this label is available for use
                        usedLabels.clear
                            (Integer.parseInt(label.substring(1)));
                    }
                    labelMap[g] = 0;
                    if (debug) {
                        System.out.println
                            ("## remove degenerate R-group " + label);
                    }
                    --rgd.rgroupCount;
                }
            }
        }
        
        /*try {
          System.out.println(scaffold.exportToFormat("cxsmiles:q"));
          } catch (MolExportException e) {
          e.printStackTrace();
          }*/
        hits = null;
        for (int r = 0; r < rgd.rgroups.length; ++r) {
            Set<MolAtom> adjustLabel = new HashSet<MolAtom>();
        
            for (int g = 0; g < rgd.rgroupCount; ++g) {
                Molecule ligand = rgd.rgroups[r][g];
                if (ligand != null) {
                    MolAtom[] atoms = ligand.getAtomArray();
                    BitSet attachments = new BitSet (atoms.length);
                    for (int i = 0; i < atoms.length; ++i) {
                        MolAtom a = atoms[i];
                        if (a.getAtomMap() > 0) {
                            adjustLabel.add(a);
                            attachments.set(i);
                        }
                    }
        
                    if (attachments.cardinality() == 1) {
                        // don't adjust label for single attachment
                        adjustLabel.remove
                            (atoms[attachments.nextSetBit(0)]);
                    }
                }
            }
        
            for (MolAtom a : adjustLabel) {
                // make sure the atom map is proper
                int map = a.getAtomMap();
                if (map > 0) { // unnecessary check?
                    a.setAtomMap(labelMap[map-1]);
                    if (debug) {
                        System.out.println
                            (molseq.get(rgd.rows[r]).getName() 
                             + ": adjusting map " + map + " => "
                             + labelMap[map-1]);
                    }
                }
            }
        }
        
        // do another pass over (sigh) to handle attachment within scaffold
        for (int r = 0; r < rgd.rgroups.length; ++r) {
            for (int g = 0; g < rgd.rgroupCount; ++g) {
                Molecule ligand = rgd.rgroups[r][g];
                if (ligand == null) {
                    continue;
                }
        
                MolAtom[] atoms = ligand.getAtomArray();
                BitSet attachments = new BitSet (atoms.length);
                for (int i = 0; i < atoms.length; ++i) {
                    MolAtom a = atoms[i];
                    if (a.getAtomMap() > 0) {
                        attachments.set(i);
                    }
                }
        
                if (attachments.cardinality() == 1) {
                    int att = attachments.nextSetBit(0);
                    MolAtom a = atoms[att];
                    
                    int scaffoldAtom = -1;
                    int seq = a.getExtraLabelSetSeq() - 1;
                    // see if this atom belongs to the scaffold
                    for (int i = 0; 
                         i < rgd.hits[r].length && scaffoldAtom < 0; ++i) {
                        if (rgd.hits[r][i] == seq) {
                            /* this assumes that the atom index is the 
                             * same as the original core!!!!
                             */
                            scaffoldAtom = i;
                        }
                    }
                    //a.setExtraLabelSetSeq(-1);
                    
                    if (scaffoldAtom < 0) {
                        MolAtom ra = new MolAtom (MolAtom.PSEUDO);
                        ra.setAttach(MolAtom.ATTACH1);
                        ligand.add(ra);
                        ligand.add(new MolBond (ra, a));
                        
                        // only one attachment point, so no need for the
                        //  label
                        int hcount = a.getImplicitHcount();
                        if (hcount > 0) {
                            a.setImplicitHcount(hcount-1);
                        }
                    }
                    else {
                        // there should be as many R-groups as there
                        //  are bonds at this position in the scaffold
                        MolAtom sa = scaffold.getAtom(scaffoldAtom);
                        
                        BitSet nr = new BitSet (scaffold.getAtomCount());
                        for (int i = 0; i < sa.getBondCount(); ++i) {
                            MolAtom xa = sa.getBond(i).getOtherAtom(sa);
                            if (xa.getAtno() == MolAtom.RGROUP) {
                                nr.set(scaffold.indexOf(xa));
                            }
                        }
        
                        lastLabel = usedLabels.nextClearBit(1);
                        for (int i = nr.cardinality(); 
                             i < a.getBondCount(); ++i) {
                            int j = rgd.rgroupCount++; 
                            logger.info
                                ("Adding new R-group due to "
                                 +"scaffold atom "+(scaffoldAtom+1)
                                 +" matching ligand ring atom "
                                 +"at " + (seq+1) +" in "
                                 +molseq.get(rgd.rows[r]).getName());
                            rgd.rgroupLabels[j] = "R"+lastLabel;
                            MolAtom ra = addRgroup 
                                (scaffold, lastLabel, scaffoldAtom);
                            nr.set(scaffold.indexOf(ra));
        
                            usedLabels.set(lastLabel);
                            lastLabel = usedLabels.nextClearBit
                                (lastLabel+1);
                        }
        
                        /*
                          if (nr.cardinality() != a.getBondCount()) {
                          // scaffold has more R-group than the number
                          //  of bonds
                          Molecule core = scaffold.cloneMolecule();
                          core.getAtom(scaffoldAtom).setAtomMap
                          (scaffoldAtom+1);
        
                          logger.log(Level.WARNING, 
                          "There are " + nr.cardinality() 
                          +" R-rgroups at scaffold position " 
                          +(scaffoldAtom+1) +"\n"
                          +core.toFormat("cxsmarts:q")
                          + "\nbut there are "
                          +a.getBondCount() + " posible "
                          +"mappings at position "+(seq+1)
                          +" in molecule " 
                          + molseq.get(rgd.rows[r]).getName()
                          );
                          }
                        */
        
                        int j = nr.nextSetBit(0);
                        for (int i = 0; i < a.getBondCount(); ++i) {
                            MolAtom xa = a.getBond(i).getOtherAtom(a);
                            xa.setAttach(MolAtom.ATTACH1);
                            int rlabel = scaffold.getAtom(j).getRgroup();
                            xa.setAtomMap(rlabel);
                            for (int l = 0; 
                                 l < rgd.getRGroupCount(); ++l) {
                                if (rgd.getRGroupLabel(l).equals
                                    ("R"+rlabel) 
                                    && rgd.rgroups[r][l] == null) {
                                    rgd.rgroups[r][l] = ligand;
                                    Map ligs = rgd.ligands[l];
                                    if (ligs == null) {
                                        rgd.ligands[l] = ligs 
                                            = new HashMap<String, Integer>();
                                    }
                                    String lig = ligand.toFormat
                                        ("cxsmarts:q");
                                    Integer c = (Integer)ligs.get(lig);
                                    ligs.put(lig, c==null?1 : (c+1));
                                }
                            }
                            j = nr.nextSetBit(j+1);
                        }
                    } // endif attachment is part of scaffold
        
                    // reclean
                    a.setAttach(0);
                    a.setAtomMap(0);
                    
                    
                    ligand.clean(2, null);
                    ChemUtil.resetEZ(ligand);
                } // endif one attachment
                ligand.aromatize();
                ligand.dearomatize();
            } // endfor each group
        }
        
        //TODO: Determine if this is a good idea
        //this takes care of the x-groups as a secondary concern
        if(rgd.xgroupCount>0){
            int[] qatompos = new int[rgd.xgroupCount];
            int aindex=0;
            int qindex=0;
            for(MolAtom ma:scaffold.getAtomArray()){
                if(ma.getAtno()==MolAtom.LIST||ma.getAtno()==MolAtom.ANY){
                    //ma.set
                    ma.setAliasstr(rgd.xgroupLabels[qindex]);
                    ma.setExtraLabel(" ");
                    qatompos[qindex++]=aindex;
                                
                }else{
                    ma.setExtraLabel(" ");
                }
                aindex++;
            }
            Molecule[] mols=rgd.getMembers();
            for(int i=0;i<rgd.getMemberCount();i++){
                int x=0;
                for(int qpos : qatompos){
                    int matomindex = rgd.hits[i][qpos];
                    Molecule m= new Molecule();
                    MolAtom xatom=(MolAtom)mols[i].getAtom(matomindex).clone();
                    xatom.setImplicitHcount(0);
                    m.add(xatom);
                    if(false){                          
                        for(int j=0;j<mols[i].getAtom(matomindex).getBondCount();j++){
                            MolAtom matt = new MolAtom(MolAtom.ANY);
                            m.add(matt);
                            MolBond mb = new MolBond(xatom,matt,mols[i].getAtom(matomindex).getBond(j).getType());
                            m.add(mb);
                        }
                        m.clean(2, null);
                    }
                                
                                
                    rgd.xgroups[i][x++]=m;
                }
            }
        }
        
        //calculates score
        if (rgd.rgroupCount > 0 && rgd.getRowCount() > 1) {
            if(DOSETUP){
                rgd.setupColumns();
            }

            rgd.calculateScore();
            scaffold.hydrogenize(false);
            scaffold.clean(2, null);
            resetEZ (scaffold, 0);
            //fixes query atoms
            Map<MolAtom,MolAtom> realSmarts = new HashMap<MolAtom,MolAtom>();
            for(MolAtom ma:scaffold.getAtomArray()){
                if(ma.isQuery()){
                    realSmarts.put(ma, (MolAtom)ma.clone());
                    ma.setAtno(6);
                }
            }
            scaffold.dearomatize();
            //add chirality features
            for(int i=0;i<chirality.length;i++){
                if(chirality[i]>0){
                    try{
                        scaffold.setChirality(i,chirality[i]);
                    } catch(Exception e){
                        
                    }
                }
            }
            
            for(MolAtom ma:realSmarts.keySet()){
                ma.set(realSmarts.get(ma));
            }

            if (debug) {
                dump (rgd);
            }
        }
        else {
            rgd = null;
        }

        return rgd;
    } // solve ()

    /**
     * This method generates an R-group decomposition from the data already
     * provided in the molecules
     */
    public static RGroupTable solveLoad (int scaffNum, Molecule scaffold, 
                                         DataSeq<Molecule> molseq, 
                                         BitSet subset) {

        MolSearch msearch = new MolSearch ();
        msearch.setOrderSensitiveSearch(true);
        MolFpFactory fpFact = MolFpFactory.getInstance();
        scaffold.hydrogenize(false);
        scaffold.expandSgroups();
        List<Object> toRemove = new ArrayList<Object>();
        BitSet usedLabels = new BitSet ();
        int[][] assignedGroups = null;
        {
            Set<MolAtom> remove = new HashSet<MolAtom>();
            
            Map<MolAtom, BitSet> labeledRgroup = 
                new HashMap<MolAtom, BitSet>();
            
            MolAtom[] atoms = scaffold.getAtomArray();
            for (int i = 0; i < atoms.length; ++i) {
                MolAtom a = atoms[i];
                switch (a.getAtno()) {
                case MolAtom.PSEUDO:
                    // force this atom to be any
                    a.setAtno(MolAtom.ANY);
                    break;
                    
                case MolAtom.RGROUP:
                    if (a.isTerminalAtom()) {
                        // remove this atom... and make a note 
                        if (a.getBondCount() > 0) {
                           
                            MolAtom xa = a.getBond(0).getOtherAtom(a);
                            BitSet bs = labeledRgroup.get(xa);
                            if (bs == null) {
                                labeledRgroup.put(xa, bs = new BitSet ());
                            }
                            
                            bs.set(a.getRgroup());
                            usedLabels.set(a.getRgroup());
                            remove.add(a);
                            toRemove.add(a.getBond(0));
                            toRemove.add(a);
                            
                        }
                    }
                    else { // treat this as a query atom
                        a.setAtno(MolAtom.ANY);
                    }
                    break;
                    
                case MolAtom.ANY:
                case MolAtom.LIST: 
                    break;
                }
            }

            for (Object o : toRemove) {
                if(o instanceof MolBond){
                    scaffold.removeEdge((MolBond)o);
                }else{
                    scaffold.removeNode((MolAtom)o);
                }
            }

            if (!labeledRgroup.isEmpty()) {
                assignedGroups = new int[scaffold.getAtomCount()][];
                for (Map.Entry<MolAtom, BitSet> e : labeledRgroup.entrySet()) {
                    int pos = scaffold.indexOf(e.getKey());
                    BitSet bs = e.getValue();
                    assignedGroups[pos] = new int[bs.cardinality()];
                    for (int i = bs.nextSetBit(0), j = 0; i >= 0;
                         i = bs.nextSetBit(i+1), ++j) {
                        assignedGroups[pos][j] = i;
                    }
                }
            }
        }
        scaffold.dearomatize();
        scaffold.aromatize();
        int[] queryFp = fpFact.generate(scaffold);
        //scaffold.dearomatize();
        msearch.setQuery(scaffold);

        // R-group table for this scaffold
        Molecule[][][] rtable = new 
            Molecule[molseq.size()][scaffold.getAtomCount()][3];
        int[][] hits = new int[molseq.size()][];
        // r-group count at each position (there should be at most 3 groups 
        //  at a particular atom)
        int[] rgroup = new int[scaffold.getAtomCount()];

        boolean warnNonHit = true;
        if (subset == null) {
            // all bits
            subset = new BitSet (molseq.size());
            subset.flip(0, molseq.size());
            warnNonHit = false;
        }

        if (debug) {
            System.out.println("## scaffold decomposition");
            System.out.print
                (scaffold.toFormat("smarts:u") + "\t" + subset.cardinality());
        }

        // sort by complexity
        Comparator<Molecule> sorter = new Comparator<Molecule>() {
            public int compare (Molecule m1, Molecule m2) {
                return ChemUtil.complexity(m1) - ChemUtil.complexity(m2);
            }
        };

        /*
          The idea here is that there are few cases where we don't have 
          exact ligands matching at a particular position, but we do want
          "similar" ligands to align.
        */
        int[][] alignments = new int[scaffold.getAtomCount()][];
        for (int i = 0; i < alignments.length; ++i) {
            int[] fp = new int[fpFact.getNumInts()];
            for (int j = 0; j < fp.length; ++j) {
                fp[j] = ~0;
            }
            alignments[i] = fp;
        }
        int rgroupCount = scaffold.getAtomCount()*3; 
        int[] rows = new int[subset.cardinality()];
        hits = new int[rows.length][];
        int r=0;
        List<String> rGroupList = new ArrayList<String>();
        for (int i = subset.nextSetBit(0); i >= 0; i = subset.nextSetBit(i + 1)) {
            Molecule target = molseq.get(i);
            rows[r++] = i;

            StringBuilder sb = new StringBuilder();
            for (int j = 1; j <= rgroupCount; j++) {
                String cxsmiles = target.getProperty(RGroupTable.PREFIX_RGROUP + scaffNum + "R" + j);
                target.setProperty(RGroupTable.PREFIX_RGROUP + scaffNum + "R" + j, "");
                if (cxsmiles == null) {
                    rgroupCount = j - 1;
                    break;
                } else {
                    if (j > 1)
                        sb.append(";");
                    sb.append(cxsmiles);
                }
            }
            // System.out.println("R:"+sb.toString());
            rGroupList.add(sb.toString());

            boolean substructure = true;
            target.aromatize();
            int[] targetFp = fpFact.generate(target);
            for (int j = 0; j < queryFp.length && substructure; ++j) {
                substructure = (queryFp[j] & targetFp[j]) == queryFp[j];
            }

            if (!substructure) {
                subset.clear(i);
                continue;
            }

            msearch.setTarget(target);
            if (debug) {
                System.out.print(" " + target.getName() + "{");
            }
            try {

                int[][] allHits = msearch.findAll();
                if (allHits == null) {
                    System.out.println("Couldn't find:"
                                       + scaffold.toFormat("cxsmarts:u"));
                    System.out.println("Couldn't find:"
                                       + target.toFormat("cxsmiles:u"));
                } else {
                    hits[r - 1] = allHits[0];
                }
            } catch (chemaxon.sss.search.SearchException ex) {
                logger.log(Level.SEVERE, "substructure search error: query=\""
                           + scaffold.toFormat("cxsmarts") + "\" target=\""
                           + target.toFormat("smiles:q") + "\"", ex);
            }
            if (debug) {
                System.out.print("}");
            }
        }
        alignments = null;
        RGroupTable rgd = new RGroupTable (molseq);
        { 
            //System.out.println("Number of rgroups:" + rgroupCount);
            rgd.rgroupCount=rgroupCount;
            rgd.rgroups = new Molecule[rGroupList.size()][rgroupCount];
            rgd.ligands = new Map[rgroupCount];
            rgd.rgroupLabels = new String[rgroupCount];
            for(int i=0;i<rgroupCount;i++){
                rgd.rgroupLabels[i]="R"+(i+1);
            }
            for(int i=0;i<rGroupList.size();i++){
                int j=0;
                for(String smi:rGroupList.get(i).split(";")){
                    /*this gets the count of ligands like this one*/
                    Map ligs = rgd.ligands[j];
                    if(ligs==null)rgd.ligands[j]=ligs=new HashMap<String,Integer>();
                    Integer c=(Integer) ligs.get(smi);
                    ligs.put(smi, (c==null)?1:c+1);
                    Molecule m=null;
                    if(smi.length()>0 && !smi.equals(" ")){
                        try{
                            m=MolImporter.importMol(smi);
                            MolAtom[] atoms = m.getAtomArray();
                            for (MolAtom a:atoms) {
                                                
                                if(a.getAtno()==MolAtom.PSEUDO ||a.getAtno()==MolAtom.RGROUP || a.isPseudo()||a.isQuery()){
                                    //System.out.println(a.getAtno());
                                    //System.out.println("PSUDO");
                                    a.setAttach(MolAtom.ATTACH1);
                                    a.setAtno(MolAtom.PSEUDO);
                                    //a.setAliasstr("*");
                                }else{
                                    if(a.getAtomMap()>0){
                                        a.setAttach(MolAtom.ATTACH1);
                                    }
                                }
                            }
                        }catch(Exception e){
                            e.printStackTrace();
                        }
                    }
                    try{
                        rgd.rgroups[i][j]=m;
                    }catch(Exception e){
                        System.out.println("RgroupCount:"+ rgroupCount);
                        System.out.println("Rgd Size 1d:"+ rgd.rgroups.length);
                        System.out.println("Rgd Size 2d:"+ rgd.rgroups[i].length);
                        System.out.println("i:"+ i);
                        System.out.println("j:"+ j);
                        throw new ArrayIndexOutOfBoundsException();
                    }
                    j++;
                    if(j>=rgroupCount)break;
                }
            }
            //labelMap = new int[rgroupCount];
        }
        rgd.core = scaffold.cloneMolecule(); // the original input scaffold
        for (int i=toRemove.size()-1;i>=0;i--) {
            Object o=toRemove.get(i);
            if(o instanceof MolBond){
                scaffold.add((MolBond)o);
            }else{
                scaffold.add((MolAtom)o);
            }
        }
        /*
          if (scaffold.getDim() < 2) {
          scaffold.clean(2, null);
          }*/
        
        rgd.scaffold = scaffold.cloneMolecule();
        rgd.rows = rows;
        
        rgd.hits = hits;
        rgd.fp = queryFp;
        rgd.complexity = ChemUtil.complexity(rgd.core);
        
        if (rgd.rgroupCount > 0 && rgd.getRowCount() > 1) {
            //rgd.setupColumns();

            // calculate r-group score
            double sig = 0., size = rgd.core.getAtomCount();
            for (int i = 0; i < rgd.rows.length; ++i) {
                Molecule m = molseq.get(rgd.rows[i]);
                double x = m.getAtomCount() - size;
                sig += x*x;
            }
            
            if (sig > 0.) {
                /*
                  double snr = size
                  / Math.sqrt(sig/rgd.rows.length);
                  rgd.score = Math.sqrt(subset.cardinality()) 
                  * Math.log10(snr);//Math.log10(2);
                  */
                //rgd.score = Math.log10(rgd.rows.length)* Math.sqrt(snr);
                rgd.score = size*rgd.rows.length/Math.sqrt(sig);
            }

            scaffold.aromatize();
            scaffold.dearomatize();

            /*
              int[] hit = new int[rgd.core.getAtomCount()];
              for (int i = 0; i < hit.length; ++i) {
              hit[i] = i;
              }
              scaffold.partialClean(2, hit, null);
            */
           

            if (debug) {
                dump (rgd);
            }
        }
        else {
            rgd = null;
           
        }
        
        scaffold.clean(2, null);
        //MolAligner.align(rgd.core, scaffold, hit);
        resetEZ (scaffold, 0);
        

        return rgd;
    } // solve ()
    
    
    /**
     * Attempts to simplify existing r-group table to be more in line with expectations.
     * This just changes multiple substituent groups around a phenyl to be a variable attachment.
     * That is, it uses locants. 
     * 
     * @param rgroup  -- the origninal rgroup table
     * @return the rgroup simplified (and, for now, also locally modified)
     */
    public static RGroupTable markushAttach(RGroupTable rgroup){
        //TODO: fix ortho,meta,para designation to be 2,3,4,5,6
        int ATTACH = 1;
        //String[] ringPositionLocants = new String[]{null,"ORTHO","META","PARA", "META","ORTHO"};
        String[] ringPositionLocants = new String[]{null,"1","2","3","4", "5","6"};
        Set<MolAtom> remove = new HashSet<MolAtom>();
        Map<MolAtom, Integer> labeledRgroup = new HashMap<MolAtom, Integer>();
        Map<Integer,Integer> labeledRgroupAtomNumber = new HashMap<Integer, Integer>();
        //all I really want to know is whether there's a nonfused pheny group in the scaffold
        //that has 2 or more attached R-groups
        //So:
        // 1) remove R-groups from scaffold. Keep track of position
        // 2) search for smarts
        Molecule scaff=(Molecule)rgroup.getScaffold().cloneMolecule();
            
        Molecule query=null;
        MolHandler mh = new MolHandler();
        try {
            mh.setMolecule("*-[b,c,n,o,p,s:1]1[b,c,n,o,p,s;D2:2][b,c,n,o,p,s;D2:3][b,c,n,o,p,s;D2:4][b,c,n,o,p,s;D2:5][b,c,n,o,p,s;D2:6]1");
            query=mh.getMolecule();
            query.clean(2,null);
        } catch (MolFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
           
        MolAtom[] atoms = scaff.getAtomArray();
        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
                        
            if(a.getAtno()==MolAtom.RGROUP || a.getAtno()== MolAtom.LIST || a.getAtno()== MolAtom.PSEUDO || a.getAtno()== MolAtom.ANY){
                if (a.isTerminalAtom()) {
                    //System.out.println("OK...");
                    // remove this atom... and make a note 
                    if (a.getBondCount() > 0) {
                        MolAtom xa = a.getBond(0).getOtherAtom(a);
                        labeledRgroup.put(xa, a.getRgroup());
                        labeledRgroupAtomNumber.put(a.getRgroup(), i);
                        remove.add(a);
                    }
                }
            }
        }
        for(MolAtom rem:remove){
            scaff.removeNode(rem);
        }
        scaff.aromatize();
        if(debug){
            try{
                System.out.println("Testing"+ rgroup.getScaffold().exportToFormat("cxsmarts:q") +"\t" + scaff.exportToFormat("cxsmarts:q"));
            }catch(Exception e){
                        
            }
        }
        MolSearch ms = new MolSearch();
        ms.setQuery(query);
        ms.setSearchType(MolSearch.SUBSTRUCTURE);
     
        ms.setTarget(scaff);
        int hits[][]=null;
        try{
            hits=ms.findAll();
        }catch(Exception e){
            e.printStackTrace();
        }
        int changeCount=0;
        if(hits!=null){
            HashMap<Integer,Integer> groupsToRemove = new HashMap<Integer,Integer>();
            List<List<Molecule>> groupsToAdd = new ArrayList<List<Molecule>>();
            if(debug){
                System.out.println("\t" + "Looks possible...");
            }
            for(int[] inst:hits){
                Molecule ligEx = query.cloneMolecule();
                        
                //the mapping from R-group label to 2,3,4,5,6 positions
                Map<Integer,Integer> rgroupNum = new LinkedHashMap<Integer,Integer>();
                //the mapping from R-group label to ligex/query position
                Map<Integer,Integer> rgroupLig = new LinkedHashMap<Integer,Integer>();
                        
                int maxpos=0;
                int minpos=7;
                //Set<Object> ligObs = new LinkedHashSet<Object>();
                
                for(int m=0;m<inst.length;m++){
                    MolAtom scaffAtom=scaff.getAtom(inst[m]);
                    //ligObs.add(scaffAtom);
                    int atomType=query.getAtom(m).getAtomMap();
                    if(atomType!=0){
                        ligEx.getAtom(m).setFlags(0);
                        ligEx.getAtom(m).setAtno(6);
                        ligEx.getAtom(m).setAliasstr("");
                        ligEx.getAtom(m).setAtomMap(0);
                        ligEx.getAtom(m).set((MolAtom)scaffAtom.clone());
                        ligEx.getAtom(m).setSetSeq(3);
                    }else{
                        ligEx.getAtom(m).setSetSeq(3);
                        ligEx.getAtom(m).setAliasstr("*");
                    }
                    if(labeledRgroup.containsKey(scaffAtom)){
                        if(atomType!=ATTACH && atomType!=0){
                            rgroupLig.put(labeledRgroup.get(scaffAtom),m);
                            rgroupNum.put(labeledRgroup.get(scaffAtom),atomType);
                            if(atomType>=maxpos){
                                maxpos=atomType;
                            }
                            if(atomType<=minpos){
                                minpos=atomType;
                            }
                        }
                    }
                }
                if(rgroupNum.size()>1){ // interesting enough to be variable
                        
                    //minimize starting position
                    if(minpos>=4){
                        for(int rgroupPos:rgroupNum.keySet()){
                            int p=rgroupNum.get(rgroupPos);
                            rgroupNum.put(rgroupPos,8-p);
                        }
                    }
                         
                        
                    if(debug){
                        System.out.println("\t\t" + "Has enough rgroups..." + ligEx.toFormat("cxsmarts:q"));
                    }
                    if (new HashSet<Integer>(rgroupNum.values()).size() > 1) {
                        if(debug){
                            System.out.println("\t\t\t" + "More than one subst. position");
                        }
                        boolean terminalonly = true;
                        for (int r : rgroupNum.keySet()) {
                                                
                            // if the rgroups have any fused rings, don't
                            // simplify
                            for (int i = 0; i < rgroup.getMemberCount(); i++) {
                                int acount = 0;
                                Molecule lig = rgroup.getRGroup(i, r - 1);
                                if (lig != null) {
                                    for (MolAtom ma : lig.getAtomArray()) {
                                        if (ma.getAtno() == MolAtom.ANY
                                            || ma.getAtno() == MolAtom.PSEUDO
                                            || ma.getAtno() == MolAtom.RGROUP) {
                                            acount++;
                                        }
                                    }
                                    if (acount <= 0) {
                                        terminalonly = false;
                                        break;
                                    }
                                }
                                                                
                            }
                            if (!terminalonly) {
                                break;
                            }
                        }
                        if (terminalonly) { //let's do it
                            changeCount++;
                            for (int r : rgroupNum.keySet()) {
                                groupsToRemove.put(r,changeCount);
                            }
                            if(debug){
                                System.out.println("\t\t\t\t" + "And all ligands are terminal ...");
                            }
                            List<Molecule> newAttachments = new ArrayList<Molecule>();
                            for (int i = 0; i < rgroup.getMemberCount(); i++) {
                                Molecule ligadd = ligEx.cloneMolecule();
                                Map<Integer,MolAtom> attns = new HashMap<Integer,MolAtom>();
                                for (int r : rgroupNum.keySet()) {
                                    MolAtom attn = ligadd.getAtom(rgroupLig.get(r));
                                    attns.put(r,attn);
                                }
                                Map<Integer,String> subnames=new LinkedHashMap<Integer,String>();
                                int keepbits=0;
                                int changebits=0;
                                for (int r : rgroupNum.keySet()) {
                                    int SUB=rgroupNum.get(r);
                                    Molecule lig = rgroup.getRGroup(i, r - 1);
                                                                        
                                    if (lig != null) {
                                        keepbits=keepbits|1<<SUB;
                                        changebits=changebits|1<<(8-SUB);
                                                                                
                                        lig=lig.cloneMolecule();
                                        MolAtom attn = attns.get(r);
                                        MolAtom attStar = null;
                                        MolAtom att = null;

                                        for (MolAtom ma : lig.getAtomArray()) {
                                            if (ma.getAtno() == MolAtom.ANY
                                                || ma.getAtno() == MolAtom.PSEUDO
                                                || ma.getAtno() == MolAtom.RGROUP) {
                                                attStar = ma;
                                                att = ma.getBond(0)
                                                    .getOtherAtom(ma);
                                            }
                                        }
                                                                                
                                        if(att.getAtno()!=1){
                                            subnames.put(SUB,"-"+SimpleNamer.getName(lig));
                                            lig.removeNode(attStar);
                                            ligadd.fuse(lig);
                                            ligadd.add(new MolBond(att, attn));
                                        }
                                    }
                                }
                                //not exhaustive
                                if(keepbits>=changebits){
                                    Map<Integer,String> tmpnames = new LinkedHashMap<Integer,String>();
                                    for(int loc:subnames.keySet()){
                                        tmpnames.put(8-loc,subnames.get(loc));
                                    }
                                    subnames=tmpnames;
                                }
                                String name="<html>";
                                boolean first=true;
                                for(int n=0;n<7;n++){
                                    String sname = subnames.get(n);
                                    if(sname!=null){
                                        sname = n+sname;
                                        if(first)name+=sname;
                                        else name+="<br/>" +sname;
                                        first=false;
                                    }
                                }
                                name+="</html>";
                                ligadd.setName(name);
                                ligadd.clean(2,null);
                                newAttachments.add(ligadd);
                            }
                            groupsToAdd.add(newAttachments);
                                                        
                        }
                                                
                    }
                                        
                }
            }//all done
            if(!groupsToAdd.isEmpty() && !groupsToRemove.isEmpty()){
                Molecule[][] oldRgroups = rgroup.rgroups;
                Molecule[][] newRgroups = new Molecule[oldRgroups.length][rgroup.rgroupCount
                                                                          - groupsToRemove.size() + groupsToAdd.size()];
                Map<Integer, Integer> oldToNew = new HashMap<Integer, Integer>();
                int oldkeep=rgroup.rgroupCount
                    - groupsToRemove.size();
                int t = 0;
                for (int i = 0; i < rgroup.rgroupCount; i++) {
                    int nmap = -1;

                    if (!groupsToRemove.containsKey(i+1)) {
                        nmap = t;
                        t++;
                    } else {
                        nmap = 0 - groupsToRemove.get(i+1);
                    }
                    oldToNew.put(i, nmap);
                }
                for (int j = 0; j < rgroup.rgroupCount; j++) {
                    int newIndex = oldToNew.get(j);
                    if (newIndex >= 0) {
                        for (int i = 0; i < oldRgroups.length; i++) {
                            newRgroups[i][newIndex] = oldRgroups[i][j];
                        }
                    }
                }
                int start = rgroup.rgroupCount - groupsToRemove.size();
                for (int i = 0; i < groupsToAdd.size(); i++) {
                    for(int j:groupsToRemove.keySet()){
                        int c=groupsToRemove.get(j);
                        if(c==i+1){
                            oldToNew.put(j-1,i + start);
                        }
                    }
                                        
                    List<Molecule> add = groupsToAdd.get(i);
                    for (int j = 0; j < add.size(); j++) {
                        newRgroups[j][i + start] = add.get(j);
                    }
                }
                Map<Integer,MolBond> rbond =new HashMap<Integer,MolBond>(); 
                // update scaffold to look like what we have now:
                for (MolAtom ma : labeledRgroup.keySet()) {
                    int rindex=labeledRgroup.get(ma)-1;
                    int newMap  = oldToNew.get(rindex);
                    //int newMap = newmapset.iterator().next() + 1;
                    if (newMap <oldkeep) {
                        MolAtom ma2 = new MolAtom(MolAtom.RGROUP);
                        ma2.setRgroup(newMap+1);
                        scaff.add(ma2);
                        scaff.add(new MolBond(ma2, ma));
                    }else{
                        if(ma.getBondCount()<3){
                            if(ma.getBond(0).getOtherAtom(ma).getBondCount()<3)
                                rbond.put(newMap,ma.getBond(0));
                            else
                                rbond.put(newMap,ma.getBond(1));
                        }
                        //ma.setExtraLabel("R"+(newMap+1));
                        //ma.setAliasstr(" ");
                    }
                }
                                
                scaff.clean(2, null);
                try{
                    scaff.dearomatize();
                }catch(Exception e){
                    logger.warning("Problem dearomatizing markush.");
                }
                Molecule rattach = new Molecule();
                rattach.add(new MolAtom(MolAtom.ANY));
                rattach.add(new MolAtom(MolAtom.RGROUP));
                rattach.getAtom(0).setAliasstr(" ");
                rattach.add(new MolBond(rattach.getAtom(0),rattach.getAtom(1)));
                for(int i:rbond.keySet()){
                    Molecule r2 = rattach.cloneMolecule();
                    r2.getAtom(1).setAliasstr("R"+(i+1));
                    scaff.fuse(r2);
                                        
                    MolAtom ma1=rbond.get(i).getAtom1();
                    MolAtom ma2=rbond.get(i).getAtom2();
                    MolAtom ma3=null;
                    ma3=ma1.getBond(0).getOtherAtom(ma1);
                    if(ma3.equals(ma2)){
                        ma3=ma1.getBond(1).getOtherAtom(ma1);
                    }
                    //ma3.setAliasstr("t");
                    double x1=ma1.getX();
                    double x2=ma2.getX();
                    double x3=ma3.getX();
                    double y1=ma1.getY();
                    double y2=ma2.getY();
                    double y3=ma3.getY();
                    double dx=x2-x1;
                    double dy=y2-y1;
                    double cx=(x2+x1)/2;
                    double cy=(y2+y1)/2;
                    double nx1=cx-dy;
                    double ny1=cy+dx;
                    double nx2=cx+dy;
                    double ny2=cy-dx;
                    double dist1=Math.sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3));
                    double dist2=Math.sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
                    r2.getAtom((dist1>dist2)?1:0).setXY(nx2, ny2);
                    r2.getAtom((dist1>dist2)?0:1).setXY(nx1, ny1);
                }
                rgroup.rgroupCount = newRgroups[0].length;
                rgroup.rgroups = newRgroups;
                rgroup.scaffold = scaff;
                rgroup.calculateScore();
                                
            }
                
                
        }else{
            return null;
        }
        return null;
    }
    
    /**
     * perform extension on the given rgroup scaffold.  if no extension is
     * possible, then null is returned.
     */
    public static RGroupTable extend (RGroupTable rgroup) {
        if (!rgroup.isExtensible()) {
            return null;
        }

        DefaultAtomComparator comparator = new DefaultAtomComparator (true);
        MCSMaxClique mcs = new MCSMaxClique ();
        // treat query atoms as literals
        mcs.setAtomComparator(comparator);

        MCSMaxClique mcs2 = new MCSMaxClique ();
        mcs2.setAtomComparator(new DefaultAtomComparator (true));

        MolFpFactory fpFact = MolFpFactory.getInstance();
        int ng = rgroup.getRGroupCount();

        RGroupTable gtab = null;
        Molecule scaffold = rgroup.getScaffold().cloneMolecule();
        logger.info("Extending scaffold: " + scaffold.toFormat("cxsmarts:q"));

        // extension at each R-group atom
        Map<MolAtom, Molecule> extensions = new HashMap<MolAtom, Molecule>();

        for (int i = 0; i < ng; ++i) {
            System.out.print(rgroup.getRGroupLabel(i) + ": ");

            if (!rgroup.extension.get(i)) {
                System.out.println();
                continue;
            }
            
            // R-group Id
            int R = Integer.parseInt
                (rgroup.getRGroupLabel(i).substring(1));

            // now try to extend this rgroup
            Molecule[] ligands = {};
            {
                Map<String, Molecule> ligvec = new HashMap<String, Molecule>();
                for (int r = 0; r < rgroup.getRowCount(); ++r) {
                    Molecule lig = rgroup.getRGroup(r, i);
                    if (lig != null) {
                        String key = lig.toFormat("cxsmarts:q");
                        if (!ligvec.containsKey(key)) {
                            ligvec.put(key, lig);
                        }
                    }
                }
                ligands = ligvec.values().toArray(new Molecule[0]);
            }

            rgroup.extension.clear(i);
            if (ligands.length < 2) {
                continue;
            }
            
            Molecule core = null;
            int[] fpCore = null;
            boolean extensible = true;
            int qatt, tatt; // number of query and target attachments (resp)

            /**
             * perform all pairwise mcs and keep the smallest mcs.
             * it's possible for the ring ligand to not have explicit 
             * attachment points.  handle this special case by treating
             * the pseudo atom as equivalent to MolAtom.ANY in the mcs
             * matching.
             */
            for (int j = 0; j < ligands.length && extensible; ++j) {
                mcs.setQuery(ligands[j]);

                MolAtom[] atoms = ligands[j].getAtomArray();
                qatt = 0;
                for (int k = 0; k < atoms.length; ++k) {
                    if (atoms[k].getAttach() != 0) {
                        int map = atoms[k].getAtomMap();
                        if (map == 0 || map == R) {
                            //mcs.addQueryConstraint(k);
                        }
                        ++qatt;
                    }
                }

                for (int k = j+1; k < ligands.length; ++k) {
                    mcs.setTarget(ligands[k]);

                    atoms = ligands[k].getAtomArray();
                    tatt = 0;
                    for (int n = 0; n < atoms.length; ++n) {
                        if (atoms[n].getAttach() != 0) {
                            int map = atoms[n].getAtomMap();
                            if (map == 0 || map == R) {
                                //mcs.addTargetConstraint(n);
                            }
                            ++tatt;
                        }
                    }

                    // if there is more than one attachment points, then
                    //   we turn off literal matching
                    comparator.setLiteral(qatt == 1 && tatt == 1);
                    
                    if (mcs.search()) {
                        // the core must contain the attachment point
                        Molecule mos = mcs.getResultAsMolecule();
                        int[] fp = fpFact.generate(mos);
                        if (core == null) {
                            core = mos;
                            fpCore = fp;
                            mcs2.setQuery(core);
                        }
                        else { // see if the core overlap with the mos
                            int sup = 0, sub = 0;
                            for (int n = 0; n < fp.length; ++n) {
                                int mask = fp[n] & fpCore[n];
                                if (mask == fp[n]) ++sub;
                                if (mask == fpCore[n]) ++sup;
                            }
                            if (sup == fp.length || sub == fp.length) {
                                // keep the smaller core
                                if (sub == fp.length) {
                                    core = mos;
                                    fpCore = fp;
                                    mcs2.setQuery(core);
                                }
                                else { // core is already smaller
                                }
                            }
                            else {
                                // last resort... 
                                mcs2.setTarget(mos);
                                if (mcs2.search()) {
                                    core = mcs2.getResultAsMolecule();
                                    fpCore = fpFact.generate(core);
                                    mcs2.setQuery(core);
                                }
                                else {
                                    System.out.println
                                        ("no mcs found between "
                                         + mcs2.getQuery()
                                         .toFormat("cxsmarts:q") 
                                         + " and " 
                                         + mcs2.getTarget()
                                         .toFormat("cxsmarts:q"));
                                    
                                    // bail out... we can't extend this group
                                    //  any further
                                    extensible = false;
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        System.out.println
                            ("no mcs found between ligands "
                             + mcs.getQuery().toFormat("cxsmarts:q") 
                             + " and " 
                             + mcs.getTarget().toFormat("cxsmarts:q")
                             + "... bailing out");

                        // bail out...
                        extensible = false;
                        break;
                    }
                }
            }

            if (extensible && core != null) {
                System.out.println(core.toFormat("cxsmarts:q"));
                // now extend the ligand
                for (MolAtom a : scaffold.getAtomArray()) {
                    if (a.getAtno() == MolAtom.RGROUP) {
                        if (a.getRgroup() == R) {
                            extensions.put(a, core);
                        }
                    }
                }
                if (extensions.isEmpty()) {
                    logger.log(Level.WARNING, 
                               "Unable to extend R"+R+ " with " 
                               + core.toFormat("cxsmarts:q"));
                }
            }
        } // endfor each R-group

        if (!extensions.isEmpty()) {
            for (Map.Entry<MolAtom, Molecule> e : extensions.entrySet()) {
                MolAtom ra = e.getKey();
                MolAtom xa = ra.getBond(0).getOtherAtom(ra);
                scaffold.removeNode(ra);
                
                // now find the attachment point in the core
                Molecule core = e.getValue();
                MolAtom rb = null;
                MolAtom[] atoms = core.getAtomArray();
                for (int j = 0; j < atoms.length; ++j) {
                    if (atoms[j].getAtno() == MolAtom.PSEUDO) {
                        if (rb != null) {
                            logger.log
                                (Level.WARNING, "Extension core " 
                                 + core.toFormat("cxsmarts:q") 
                                 + " contains more than one attachments!");
                        }
                        rb = atoms[j];
                    }
                    /*
                      atoms[j].setSetSeq(0);
                      atoms[j].setExtraLabelSetSeq(0);
                    */
                }
                
                if (rb == null) {
                    logger.log(Level.WARNING, "Extension core " 
                               + core.toFormat("cxsmarts:q") 
                               + " has no attachment point!");
                }
                else {
                    MolAtom xb = rb.getBond(0).getOtherAtom(rb);
                    core.removeNode(rb);
                    scaffold.fuse(core);
                    scaffold.add(new MolBond (xa, xb));
                }
            }

            // now remove all r-group node
            for (MolAtom a : scaffold.getAtomArray()) {
                if (a.getAtno() == MolAtom.RGROUP) {
                    scaffold.removeNode(a);
                }
            }

            BitSet subset = new BitSet ();
            for (int j = 0; j < rgroup.rows.length; ++j) {
                subset.set(rgroup.rows[j]);
            }

            scaffold.clean(2, null);
            scaffold.aromatize();
            resetEZ (scaffold, MolBond.TRANS|MolBond.CIS);
            //ChemUtil.resetEZ(scaffold);

            logger.info("Solving for new scaffold "
                        + scaffold.toFormat("cxsmarts:q"));

            try {
                gtab = solve (scaffold, rgroup.molseq, subset);
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        return gtab;
    }
    static MolAtom addRgroup (Molecule scaffold, int g, int atom) {
        MolAtom ra = new MolAtom (MolAtom.RGROUP);
        ra.setRgroup(g);
        scaffold.add(ra);

        MolAtom xa = scaffold.getAtom(atom);
        int charge = xa.getCharge();
        scaffold.add(new MolBond (ra, xa));

        if (xa.getCharge() != charge) {
            xa.setCharge(charge);
        }

        resetEZ (scaffold, MolBond.TRANS|MolBond.CIS);
        //ChemUtil.resetEZ(scaffold);
        return ra;
    }

    static void resetEZ (Molecule ligand, int flags) {
        MolBond[] bonds = ligand.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond b = bonds[i];
            if (/*ligand.isRingBond(i) ||*/ b.getType() == 2) {
                MolAtom a1 = b.getCTAtom1();
                MolAtom a4 = b.getCTAtom4();
                if (a1 != null && a4 != null) {
                    b.setStereo2Flags(a1, a4, flags);
                    
                }
            }
        }
    }

    protected static void dump (RGroupTable rgd) {
        try {
        System.out.println(" " + rgd.scaffold.toFormat("cxsmarts"));
        } catch (Exception ex) {ex.printStackTrace(); }
        
                
        for (int i = 0; i < rgd.getRowCount(); ++i) {
            //Molecule m = (Molecule)rgd.getValueAt(i, 1);
            String name = (String)rgd.getValueAt(i, 0);
            System.out.println(name+":");
            for (int r = 0; r < rgd.getRGroupCount(); ++r) {
                Molecule lig = rgd.getRGroup(i, r);
                if (lig != null) {
                    System.out.println
                        (String.format("%1$5s: ", 
                                       rgd.getRGroupLabel(r))
                         +lig.toFormat("cxsmarts"));
                }
            }
        }
        
        System.out.println("## unique ligands");
        for (int i = 0; i < rgd.getRGroupCount(); ++i) {
            System.out.print(rgd.getRGroupLabel(i)+":");
            final Map ligands = rgd.ligands[i];
            if (ligands != null) {
                Vector sorted = new Vector(ligands.keySet()); 
                Collections.sort(sorted, new Comparator () {
                        public int compare (Object o1, 
                                            Object o2) {
                            Integer n1 = (Integer)ligands.get(o1);
                            Integer n2 = (Integer)ligands.get(o2);
                            return n2-n1;
                        }
                    });
                
                for (Object s : sorted) {
                    System.out.print(" {" + s+","+ligands.get(s)+"}");
                }
                System.out.println();
            }
        }
        System.out.println();
    }
    protected boolean getBestHitsAndLigands(Molecule target, int[][] amorph,int[][] allHits, int[][] alignments,Molecule[] ligands, int[] hit, int[] chirality){
        if (target.getDim() < 2)
            target.clean(2,null);
        resetEZ(target,0);
        
        if (allHits != null && allHits.length > 0) {
            //int[] hit = null;
            //Molecule[] ligands = null;
            int[] targetchiral = new int[allHits[0].length];
        
            int mincard = Integer.MAX_VALUE;
            int besta=0;
            for (int j = 0; j < allHits.length; ++j) {
                int[] ttargetchiral = new int[allHits[0].length];
                Molecule[] ligs = getLigands (target, allHits[j],ttargetchiral);
                StringBuilder sb = new StringBuilder();
                for(int a=0;a<targetchiral.length;a++){
                    if(targetchiral[a]!=0)
                        sb.append(targetchiral[a]+";");
                }
                try{
                    //System.out.println(target.exportToFormat("cxsmiles:q")+ "\t"+ sb.toString());
                }catch(Exception e){
                                        
                }
                for (int a = 0; a < amorph.length; ++a) {
                    //hit goes from query to target, so permute on automorphism
                    //That is, hit[i] becomes hit[amorph[a][i]]
                    int[] amorphhit= new int[allHits[j].length];
                    Molecule[] amorphlig= new Molecule[ligs.length];
                    for(int k=0;k<allHits[j].length;k++){
                        amorphhit[k]=allHits[j][amorph[a][k]];
                        amorphlig[k]=ligs[amorph[a][k]];
                    }
                    int card = 0;
                    for (int k = 0; k < allHits[j].length; ++k) {
                        if (amorphlig[k] != null) {
                            int[] fp = fpFact.generate(amorphlig[k]);
                            boolean hyd=true;
                            for (int p = 0; p < fp.length; ++p) {
                                if(fp[p]!=hFP[p]){hyd=false;break;}
                            }
                            if(hyd){fp[0]=23556;}//why does this work?
                            boolean different=false;
                            for (int p = 0; p < fp.length; ++p) {
                                card += ChemUtil.countBits
                                    (alignments[k][p] & ~fp[p]);
                                if(alignments[k][p] !=fp[p])different=true;
                            }
                            if(!different)card--;
                        }
                    }
                    if (card < mincard) {
                        besta=a;
                        //hit = amorphhit;
                        for(int i=0;i<hit.length;i++){
                            hit[i]=amorphhit[i];
                        }
                        for(int i=0;i<ligands.length;i++){
                            ligands[i]=amorphlig[i];
                        }
                        mincard = card;
                        targetchiral=ttargetchiral;
                    }
                }
            }
                
            if (hit == null || ligands == null) {
                logger.log
                    (Level.SEVERE, "Hits or ligands is null; "
                     +"this should never happened!");
            }
                
            for (int j = 0; j < ligands.length; ++j) {
                if (ligands[j] != null) {
                    for (MolAtom ma : ligands[j].getAtomArray()) {
                        int amap = ma.getAtomMap() - 1;
                        if (amap >= 0)
                            ma.setAtomMap(amorph[besta][amap] + 1);
                    }
                    int[] fp = fpFact.generate(ligands[j]);
                    boolean hyd = true;
                    for (int p = 0; p < fp.length; ++p) {
                        if (fp[p] != hFP[p]) {
                            hyd = false;
                            break;
                        }
                    }
                    if (hyd) {
                        fp[0] = 23556;
                    }
                    for (int k = 0; k < fp.length; ++k) {
                        alignments[j][k] &= fp[k];
                    }
                }

            }
            //targetchiral=getStereo(target,hit);
            for (int j = 0; j < chirality.length; j++) {
                int chi = chirality[amorph[besta][j]];
                int stereo = targetchiral[j];
                if (chi != 0 ) {                // is chiral/unassigned
                    if (chi == -1) {    // unassigned
                        chi = stereo;
                    } else { // assigned chiral
                        if (stereo != chi) { // disagreement, default
                            // achiral
                            chi = 0;
                            if (debug) {
                                String old = "";
                                String newo = "";
                                switch (chi) {
                                case MolAtom.CHIRALITY_R:
                                    old = "R";
                                    break;
                                case MolAtom.CHIRALITY_S:
                                    old = "S";
                                    break;
                                default:
                                    old = "OTHER";
                                    break;
                                }
                                switch (stereo) {
                                case MolAtom.CHIRALITY_R:
                                    newo = "R";
                                    break;
                                case MolAtom.CHIRALITY_S:
                                    newo = "S";
                                    break;
                                default:
                                    newo = "OTHER";
                                    break;
                                }
                            }
                        }
                    }
                    chirality[amorph[besta][j]] = chi;
                }
            }
            if (debug) {
                System.out.print(mincard);
            }
            return true;
        }
        
        return false;

    }

    protected static Molecule[] getLigands (Molecule m, int[] hit,int[] chiral) {
        
        Molecule mol = m.cloneMolecule();
        mol.dearomatize();
        //mol.hydrogenize(true);

        MolAtom[] atoms = mol.getAtomArray();

        int[] rmap = new int[atoms.length];     // all set to -1
        for (int i = 0; i < rmap.length; ++i) {
            rmap[i] = -1;
            atoms[i].setAtomMap(0);
            atoms[i].setExtraLabelSetSeq(i+1);
        }
        //hits go from atom[query#]=target#
        
        //rmap set from each target atom to query rmap[target#]=query#
        for (int i = 0; i < hit.length; ++i) {
            rmap[hit[i]] = i;
        }

        BitSet[] rings = new BitSet[mol.getAtomCount()];
        
        int[][] sssr = mol.getSSSR();
        for (int i = 0; i < rings.length; ++i) {
            rings[i] = new BitSet (sssr.length);
        }
        
        for (int j = 0; j <  sssr.length; ++j) {
            int[] r = sssr[j];
            for (int i = 0; i < r.length; ++i) {
                rings[r[i]].set(j);
            }
        }
        final int REMOVE_TYPE_OUTER = 1;
        final int REMOVE_TYPE_INNER = 2;
        
        Map<MolBond,Integer> removeType = new HashMap<MolBond,Integer>();
        //Set<MolBond> remove = new HashSet<MolBond>();

        int[] rsizes = mol.getSmallestRingSizeForIdx();
        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < bonds.length; ++i) {
            MolBond bond = bonds[i];
            int a1 = mol.indexOf(bond.getAtom1());
            int a2 = mol.indexOf(bond.getAtom2());
            if (rmap[a1] >= 0 && rmap[a2] >= 0) { //both atoms part of scaffold
            }
            else if (rmap[a1] < 0 && rmap[a2] >= 0) {
                // attachment point
                if (mol.isRingBond(i)) {
                    MolAtom a = bond.getAtom2();
                    for (int j = 0; j < a.getBondCount(); ++j) {
                        MolBond b = a.getBond(j);
                        MolAtom xa = b.getOtherAtom(a);
                        if (!rings[a1].intersects(rings[mol.indexOf(xa)])) {
                            removeType.put(b, (rmap[mol.indexOf(xa)]>=0)?
                                           REMOVE_TYPE_INNER:
                                           REMOVE_TYPE_OUTER
                                           );
                            //remove.add(b);
                        }
                    }
                    atoms[a2].setAtomMap(rmap[a2]+1);
                }
                else {
                    atoms[a1].setAtomMap(rmap[a2]+1);
                    removeType.put(bond,REMOVE_TYPE_OUTER);
                    //remove.add(bond);
                }
            }
            else if (rmap[a1] >= 0 && rmap[a2] < 0) {
                // attachment point
                if (mol.isRingBond(i)) {
                    MolAtom a = bond.getAtom1();
                    for (int j = 0; j < a.getBondCount(); ++j) {
                        MolBond b = a.getBond(j);
                        MolAtom xa = b.getOtherAtom(a);
                        if (!rings[a2].intersects(rings[mol.indexOf(xa)])) {
                            removeType.put(b, (rmap[mol.indexOf(xa)]>=0)?
                                           REMOVE_TYPE_INNER:
                                           REMOVE_TYPE_OUTER
                                           );
                            //remove.add(b);
                        }
                    }
                    atoms[a1].setAtomMap(rmap[a1]+1);
                }
                else {
                    atoms[a2].setAtomMap(rmap[a1]+1);
                    removeType.put(bond, REMOVE_TYPE_OUTER);
                    //remove.add(bond);
                }
            }
        }

        // now do another pass over to remove all mapped bonds that are
        //  not connected to any mapped atoms
        for (int i = 0; i < bonds.length; ++i) {
            MolBond bond = bonds[i];
            int a1 = mol.indexOf(bond.getAtom1());
            int a2 = mol.indexOf(bond.getAtom2());
            
            boolean del = true;
            if (mol.isRingBond(i)) {
                BitSet bs = (BitSet)rings[a1].clone();
                bs.and(rings[a2]);
                for (int j = bs.nextSetBit(0); j >= 0 && del; 
                     j = bs.nextSetBit(j+1)) {
                    int[] r = sssr[j];
                    for (int k = 0; k < r.length && del; ++k) {
                        // if this ring contains at least one atom that 
                        // doesn't match the scaffold, then keep this
                        // ring in tact
                        if (rmap[r[k]] < 0) {
                            del = false;
                        }
                    }
                }
            }

            if (del) {
                if (rmap[a1] >= 0 && rmap[a2] >= 0 
                    && atoms[a1].getAtomMap() == 0 
                    && atoms[a2].getAtomMap() == 0) { //basically, remove all scaffold bonds.
                    //Unless it's part of a
                    removeType.put(bond, REMOVE_TYPE_INNER);
                    //remove.add(bond);
                }
            }
        }
        rings = null;

        //remove each outer bond 
        for (MolBond bond : removeType.keySet()) {
            switch(removeType.get(bond)){
            case REMOVE_TYPE_OUTER:
                mol.removeEdge(bond);
                break;
            case REMOVE_TYPE_INNER:
            default:
            }
        }
        for(int i=0;i<chiral.length;i++){
            chiral[i]=mol.getChirality(hit[i]);
        }
        
        // add back bonds one at a time, calculate stereo flags, then re-remove
        // this is only useful if the r-group attached is degenerate. Otherwise,
        // R,S are meaningless.
        // Doesn't work right now ... not sure why. Apparently removing bonds
        // from a molecule seems to affect the chirality gen
        /*try{
          System.out.println("removed:" + mol.exportToFormat("cxsmarts:q"));
          }catch(Exception e){
                
          }*/
        
        for(MolBond bond:removeType.keySet()){
            if(removeType.get(bond)==REMOVE_TYPE_OUTER){
                mol.add(bond);
                int ma1=mol.indexOf(bond.getAtom1());
                int ma2=mol.indexOf(bond.getAtom2());
                if(rmap[ma1]>=0){ //atom1 is scaffold
                    chiral[rmap[ma1]]=mol.getChirality(ma1);
                }else if(rmap[ma2]>=0){
                    chiral[rmap[ma2]]=mol.getChirality(ma2);
                }
                mol.removeEdge(bond);
            }
        }
        
        //add psuedo atoms
        for (MolBond bond : removeType.keySet()) {
            if(removeType.get(bond)==REMOVE_TYPE_INNER){
                mol.removeEdge(bond);
            }
            MolAtom ma1 = bond.getAtom1();
            MolAtom ma2 = bond.getAtom2();
            MolAtom mtemp1=new MolAtom(MolAtom.PSEUDO);
            MolAtom mtemp2=new MolAtom(MolAtom.PSEUDO);
            mol.add(mtemp1);
            mol.add(mtemp2);
            mol.add(new MolBond(mtemp1,ma1,bond.getType()));
            mol.add(new MolBond(mtemp2,ma2,bond.getType()));
         
        }
        

        Molecule[] ligands = new Molecule[hit.length];
        // now we should have ligands and one fragment corresponding
        //   to the core scaffold!
        Molecule[] frags = mol.convertToFrags();
        for (int i = 0; i < frags.length; ++i) {
            Molecule f = frags[i];
            for (MolAtom a : f.getAtomArray()) {
                int map = a.getAtomMap() - 1;
                if (map < 0) {
                }
                else { // attachment point...
                    if (ligands[map] != null) {
                        // multiple attachments at this position
                        try{
                            //System.out.println("FUSE! : " + ligands[map].exportToFormat("cxsmarts:q") + "," + f.exportToFormat("cxsmarts:q"));
                        }catch(Exception e){
                                
                        }
                        ligands[map].fuse(f);
                    }
                    else {
                        ligands[map] = f;
                    }
                }
            }
        }
        
        
        for (int i = 0; i < ligands.length; ++i) {
            if (ligands[i] != null) {
                ligands[i].hydrogenize(false);
            }
        }
        /*
          for (int i = 0; i < ligands.length; ++i) {
          if (ligands[i] != null) {
          try {
          String smarts = ligands[i].toFormat("smiles:q");
          ligands[i] = MolImporter.importMol(smarts);
          for (MolAtom a : ligands[i].getAtomArray()) { 
          a.setRadical(0);
          a.setCharge(0);
          }
          }
          catch (Exception ex) {
          logger.log(Level.SEVERE, 
          "Can't clean up ligand: "+ligands[i], ex);
          }
          }
          }
        */
        
        return ligands;
    }
    
    private int[] getStereo(Molecule m2,int hit[]){
        Molecule m= m2.cloneMolecule();
        int[] chiral = new int[hit.length];
        Set<Integer> keep = new HashSet<Integer>();
        List<MolAtom> toRemove =new ArrayList<MolAtom>();
        for(int i=0;i<hit.length;i++){
            keep.add(hit[i]);
            m.getAtom(hit[i]).setAtomMap(i+1);
        }
        for(int i=0;i<m.getAtomCount();i++){
            if(!keep.contains(i)){
                toRemove.add(m.getAtom(i));
            }
        }
        for(MolAtom ma:toRemove){
            m.removeNode(ma);
        }
        try{
            //System.out.println(m.exportToFormat("cxsmiles:q"));
        }catch(Exception e){
                
        }
        for(int i=0;i<chiral.length;i++){
            MolAtom ma=m.getAtom(i);
            chiral[ma.getAtomMap()-1]=m.getChirality(i);
        }
        return chiral;
    }
    
   
    
    public static void main(String[] args) throws Exception {
        String ms2="O=C(NC(CCN(CC1=CC=CS1)SC2=CC=CC=C2)CC3=CC=CC=C3)C4CN(C(=O)O4)C5=CC=CC=C5";
        //String ms="COC1=CC(=CC=C1)S(=O)(=O)N(C[C@@H](O)[C@H](CC2=CC=CC=C2)NC(=O)[C@@H]3CN(C(=O)O3)C4=CC(F)=CC=C4)CC5=CC=CS5";
        MolImporter mi = new MolImporter();
        mi.setFileName("test.sdf");
        final ArrayList<Molecule> mols = new ArrayList<Molecule>();
        
        DataSeq<Molecule> myseq=new DataSeq<Molecule>(){
            @Override
            public int size() {
                // TODO Auto-generated method stub
                return mols.size();
            }
            @Override
            public Molecule get(int index) {
                return mols.get(index);
            }
                
        };
        for(Molecule mread;(mread=mi.read())!=null;){
            mols.add(mread);
        }
        mi.close();
        Molecule q=MolImporter.importMol(ms2);
        q.aromatize();
        RGroupSolver rgs = new RGroupSolver(myseq);
        RGroupTable rgt=rgs.solve(q,null);
        for(int j=0;j<rgt.getRowCount();j++){
            System.out.print(j+"\t");
            for(int i=0;i<rgt.rgroupCount;i++){
                try{
                    System.out.print("\t"+rgt.rgroups[j][i].exportToFormat("cxsmarts:q"));
                }catch(Exception e){
                    System.out.print("\t");
                }
            }
            System.out.println();
        }
        System.out.println("Done");
        System.exit(0);
    }
}
