// $Id: RGroupGenerator.java 3864 2009-12-18 22:09:45Z nguyenda $

package gov.nih.ncgc.rgroup;

import gov.nih.ncgc.algo.graph.VFLib2;
import gov.nih.ncgc.descriptor.MolecularFramework;
import gov.nih.ncgc.model.DataSeq;
import gov.nih.ncgc.util.MolFpFactory;
import gov.nih.ncgc.util.MolStandardizer;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.logging.Logger;


import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.util.MolExportException;
import chemaxon.sss.search.MolSearch;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;


public class RGroupGenerator implements DataSeq<Molecule> {
    static final String MOLECULE_CACHE = "molecule-cache";
    protected ConcurrentMap<Integer, Molecule> cache =
        new ConcurrentHashMap<Integer,Molecule>();;
    static int cachePuts=0;
    final int cacheMax = 8000;
    final int cacheTrims = 1000;
        

    static private boolean debug = false;
    
    static {
        try {
            debug = Boolean.getBoolean("rgroup.debug");
        }
        catch (Exception ex) {
        }
    }
    static private final Logger logger = Logger.getLogger
        (RGroupGenerator.class.getName());

    private static final int MIN_FRAG_SIZE = 3;
    private static final Molecule DONE = new Molecule ();
    private static final Fragment DONEFrag = new Fragment("null","null",null);
    private static final List<Fragment> DONEFragList = new ArrayList<Fragment>();
    
    
    //private ArrayList<Molecule> molvec = new ArrayList<Molecule>();
    private ArrayList<String> molCompressVec = new ArrayList<String>();   
    
    private ConcurrentMap<String, Fragment> fragments = new ConcurrentHashMap<String, Fragment>();

    private MolStandardizer standardizer = new MolStandardizer ();
    //private RGroupSolver solver = new RGroupSolver(this,false);

    private static final int DEFAULT_maxVariableAtoms = 1;
    // max molecule size; anything else is probably not useful
    static final int DEFAULT_MOLSIZE = 50; 
    
    private boolean markushlike=true;
    private int maxVariableAtoms;
    private int minFragMemberCount = 2; 
    private int minFragmentSize = MIN_FRAG_SIZE;
    private BitSet singletons = new BitSet ();
    private int maxMolSize = DEFAULT_MOLSIZE;
    private String name;

    private ArrayList<RGroupTable> scaffolds = new ArrayList<RGroupTable>();
   
    
    private final PropertyChangeSupport propChange = 
        new PropertyChangeSupport (this);

    private BlockingQueue<Molecule> stagingQueue = 
        new ArrayBlockingQueue<Molecule>(1000);
    private BlockingQueue<Fragment> fragQueue =
        new ArrayBlockingQueue<Fragment>(10);
    private BlockingQueue<List<Fragment>> fragVariableQueue = 
        new ArrayBlockingQueue<List<Fragment>>(1000);
    
    private ExecutorService threadPool;
    final ReentrantLock lock = new ReentrantLock ();
    
    private Future[] threads;
    private boolean fullyLoaded = false;
    private int remainingMoleculesForProcessing=-1;
    
    static class Fragment{
        String _can_smiles;
        String _id;
        String _topoForm;
        MolStandardizer standardizer = new MolStandardizer ();        
        BitSet _bs;
        private int _preload_id;
        boolean generic=false;
        public Fragment(String can_smiles,String id, BitSet bs){
            this(can_smiles,id,bs,null,false);
                
        }
        public Fragment(String can_smiles,String id, BitSet bs, String topoform, boolean generic){
            _id=id;
            _can_smiles=can_smiles;
            _bs=bs;
            if(_id.startsWith(RGroupTable.PREFIX_SCAFFOLD)){
                _preload_id = Integer.parseInt(_id.replace(
                                                           RGroupTable.PREFIX_SCAFFOLD, ""));
            }else{
                _preload_id = -1;
            }
            this.generic=generic;
            _topoForm = topoform;
        }
        public void setTopoForm(String topo){
            _topoForm=topo;
        }
        public String getTopoForm(){
            if(!generic){
                if(_topoForm == null){
                    Molecule topoform = null;
                    try {
                        topoform = getTopological(_can_smiles,true,false);
                        standardizer.standardize(topoform);
                        String topHash = MolStandardizer.hashKey(topoform);
                        _topoForm = topHash;
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                                        
                }
            }
            return _topoForm;
        }
        public boolean isPreLoaded(){
            return _preload_id!=-1;
        }
        public int getPreloadID(){
            return _preload_id;
        }
        public String getSmarts(){
            return this._can_smiles;
        }
        public String getID(){
            return this._id;
        }
        public BitSet getBitSet(){
            return this._bs;
        }
        public boolean isGeneric(){
            return generic;
        }
    }
    class FragmentWorker implements Runnable {
        MolecularFramework framework =
            MolecularFramework.createMurckoInstance();
        String name;

        FragmentWorker (String name) {
            this.name = name;
            //framework.setMinFragmentSize(minFragmentSize);
            framework.setAllowBenzene(false);
            framework.setNumThreads(1);
        }

        public void run () {
            Thread.currentThread().setName(name);
            logger.info("Thread "+name+" is running...");
            try {
                MolHandler mh = new MolHandler ();
                for (Molecule mol; (mol = stagingQueue.take()) != DONE; ) {
                    boolean readError = false;
                    List<String> frags = new ArrayList<String>();
                    String prop = mol.getProperty(RGroupTable.PREFIX_SCAFFOLD);
                    if (prop != null) {
                        for (String line : prop.split("\n")) {
                            if(!line.equals(RGroupTable.SINGLETON_VALUE)){
                                String[] toks = line.split("\t");
                                if (toks.length >=2) {
                                    frags.add(toks[1]+"\t" + RGroupTable.PREFIX_SCAFFOLD+toks[0]);
                                }
                            }
                        }
                    }else{
                        try{
                            //this is a long process to standardize the molecule
                            //not sure why its necessary
                            Molecule m2=mol.cloneMolecule();
                            mh.setMolecule(m2.exportToFormat("smiles"));
                            m2 = mh.getMolecule();
                                
                            framework.setMolecule(m2, false);
                            framework.run();
                                
                            //System.out.println("Fragmenting... DONE!");
                            for (Enumeration<Molecule> en = 
                                     framework.getFragments();
                                 en.hasMoreElements(); ) {
                                Molecule m = en.nextElement();
                                m.aromatize();                      
                                //System.out.println(m2.exportToFormat("cxsmarts:q") + ":" +m.exportToFormat("cxsmarts:q"));
                                try {
                                    String toadd=m.exportToFormat("cxsmarts:q") +"\t"+ m.getName();
                                    frags.add(toadd);
                                } catch (MolExportException e) {
                                    e.printStackTrace();
                                }                   
                                    
                            }
                        }catch(Exception ex){
                            ex.printStackTrace();
                            readError=true;
                        }
                    }
                    if(!readError){
                        updateFragmentsSmiles (mol, frags);
                    }else{
                        System.out.println("Read error on Molecule");
                    }
                    if(fullyLoaded){
                        updateUnprocessedProgress();
                    }
                }
                
                logger.info("Thread "+name+" is done!");
            }
            catch (InterruptedException ex) {
                logger.info("Thread "+name+" is interrupted!");
            }
        }
    }
    class RGroupWorker implements Runnable {
        String name;
        RGroupWorker (String name) {
            this.name = name;
        }
        public void run () {
            Thread.currentThread().setName(name);
            logger.info("Thread "+name+" is running...");
            try {
                for (Fragment fragEnt; (fragEnt = fragQueue.take())
                         != DONEFrag;) {
                    /*
                    logger.info(name+": processing fragment "
                                +fragEnt._id+" "
                                +fragEnt._can_smiles+"...");
                    */
                    RGroupTable rgd = generate(fragEnt);
                    if (rgd != null) {
                        lock.lock();
                        try {
                            singletons.andNot(rgd.getMemberBitset());
                        }
                        finally {
                            lock.unlock();
                        }
                    }
                    else {
                        logger.warning("Unable to generate R-group table "
                                       +"for fragment "+fragEnt._id+" "
                                       +fragEnt._can_smiles+"!");
                        //fragments.remove(fragEnt.getID(), fragEnt);
                        fragments.remove(fragEnt.getID());
                    }
                }
            } catch (InterruptedException ex) {
                logger.info("Thread " + name + " is interrupted!");
            }
            logger.info("Thread " + name + " is done!");
        }
    }

    class VariableScaffoldWorker implements Runnable {
        String name;
        int diffAtoms=0;
        VariableScaffoldWorker (String name, int diffAtoms) {
            this.name = name;
            this.diffAtoms=diffAtoms;
        }
        public void run () {
            Thread.currentThread().setName(name);
            logger.info("Thread "+name+" is running...");
            try {
                for (List<Fragment> fragList; (fragList = fragVariableQueue.take()) != DONEFragList;) {
                    //System.out.println("got");
                    findVariableScaffolds(fragList,diffAtoms);
                }
            } catch (InterruptedException ex) {
                logger.info("Thread " + name + " is interrupted!");
            }
            logger.info("Thread " + name + " is done!");
        }
    }

    class SubsetDataSeq implements DataSeq<Molecule> {
        final Molecule[] mols;
        SubsetDataSeq (BitSet subset) {
            mols = new Molecule[subset.cardinality()];
            for (int i = subset.nextSetBit(0), j = 0;
                 i >= 0; i = subset.nextSetBit(i+1)) {
                mols[j++] = RGroupGenerator.this.get(i);
            }
        }
        public int size () { return mols.length; }
        public Molecule get (int pos) { return mols[pos]; }
    }
    
    public RGroupGenerator () {
        this (Executors.newCachedThreadPool(), 1);
    }
    public RGroupGenerator(int nthreads){
        this (Executors.newCachedThreadPool(), nthreads);
    }

    public RGroupGenerator (ExecutorService threadPool, int nthreads) {
        this(threadPool,nthreads,DEFAULT_maxVariableAtoms);
    }
    public RGroupGenerator (ExecutorService threadPool, int nthreads, int variableAtomSize) {
        this.threadPool = threadPool;
        threads = new Future[nthreads];
        maxVariableAtoms=variableAtomSize;
    }

    public int getMaxMolSize () { return maxMolSize; }
    public void setMaxMolSize (int maxMolSize) { 
        this.maxMolSize = maxMolSize; 
    }
    public int getMinFragmentSize () { 
        return minFragmentSize;
    }
    public void setMinFragmentSize (int size) { 
        this.minFragmentSize = size;
    }
    public void setMaxVaribaleAtoms (int size) { 
        this.maxVariableAtoms = size;
    }
    public void setMarushLike(boolean mlike){
        markushlike=mlike;
    }

    public void setName (String name) { this.name = name; }
    public String getName () { return name; }

    public void setNumThreads (int nthreads) {
        if (nthreads < 1) {
            throw new IllegalArgumentException
                ("Invalid number of threads: "+nthreads);
        }

        for (int i = 0; i < threads.length; ++i) {
            if (!(threads[i].isDone() || threads[i].isCancelled())) {
                throw new IllegalStateException 
                    ("Can't update thread count while there are "
                     +"running threads");
            }
        }

        // otherwise, allocate new threads
        threads = new Future[nthreads];
    }

    public void shutdown () {
        threadPool.shutdown();
    }

    //  this allows user to specify the scaffolds of interest in addition
    //  to automatically generated ones.
    synchronized public Fragment addScaffold (Molecule scaffold) {
        return addScaffold(scaffold,null,false);
    }
    synchronized public Fragment addScaffold (Molecule scaffold,BitSet bs, boolean generic) {
        String frag = scaffold.toFormat("cxsmarts:q");
        String topo = null;
        if(!generic){
            Molecule m = getTopological(scaffold, true, false);
            topo = MolStandardizer.hashKey(m);
        }
        String id=MolStandardizer.hashKey(scaffold);
        return addScaffold(id,frag,topo,bs,generic);
    }
    synchronized public Fragment addScaffold (String ID, String rep,String topo,BitSet bs, boolean generic) {
        Fragment f = new Fragment(rep,ID,bs,topo,generic);
        fragments.putIfAbsent(f.getID(),f);
        return f;
    }
    synchronized public Fragment addScaffold (Fragment f) {
        fragments.putIfAbsent(f.getID(),f);
        return f;
    }
    synchronized private RGroupTable addScaffold (RGroupTable rtab, boolean prune){
        //logger.info(Thread.currentThread()+": adding scaffold "+rtab);
        if (!prune) {
            scaffolds.add(rtab);
            return rtab;
        }

        
        // now compare this scaffold with the ones we've discovered
        //  so far to see if it's a sub- or supersubstructure of
        //  another scaffold.  if so, we remove the smaller scaffold
        Set<RGroupTable> remove = new HashSet<RGroupTable>();
        
        
        for (RGroupTable tab : scaffolds) {
            int sup = 0, sub = 0;
            for (int i = 0; i < rtab.fp.length; ++i) {
                if ((rtab.fp[i] & tab.fp[i]) == rtab.fp[i]) {
                    ++sub;
                }
                if ((rtab.fp[i] & tab.fp[i]) == tab.fp[i]) {
                    ++sup;
                }
            }
            
            if (sub == rtab.fp.length) {
                // new scaffold might be a substructure of this
                //  scaffold.  to really verify this we should
                //  be doing MolSearch on the two scaffolds.
                if (rtab.rows.length <= tab.rows.length) {
                    int k = 0;
                    for (; k < rtab.rows.length 
                             && rtab.rows[k] == tab.rows[k]; ++k)
                        ;
                    if (k == rtab.rows.length) {
                        // ok, we have with high probability that
                        //  rtab is a substructure of tab
                        // return null to indicate that we already 
                        //  have this scaffold
                        return null;
                    }
                }
            }
            else if (sup == rtab.fp.length) {
                // new scaffold might be a superstructure of this
                //  scaffold
                if (tab.rows.length <= rtab.rows.length) {
                    int k = 0;
                    for (; k < tab.rows.length 
                             && rtab.rows[k] == tab.rows[k]; ++k)
                        ;
                    if (k == tab.rows.length) {
                        // tab is a substructure of rgd, so remove it
                        remove.add(tab);
                    }
                }
            }
        }

        // remove any that was a substructure of the current scaffold
        for (RGroupTable tab : remove) {
            scaffolds.remove(tab);
        }
        scaffolds.add(rtab);

        return rtab;
    }

    public void add (String mol) {
        add (mol, true);
    }

    public void add (String mol, boolean standardize) {
        Molecule m = null;
        try {
            MolHandler mh = new MolHandler (mol);
            add (m = mh.getMolecule(), standardize);
        }
        catch (MolFormatException ex) {
            logger.log(Level.WARNING, "Not a molecule format " + mol, ex);
        }
    }

    public void add (Molecule mol) {
        add (mol, true);
    }

    public void add (Molecule mol, boolean standardize) {
        if (mol.getAtomCount() <= maxMolSize || maxMolSize <= 0) {
            if (standardize) {
                try {
                    standardizer.standardize(mol);
                }
                catch (Exception ex) {
                    logger.log(Level.WARNING, 
                               "Can't standardize " +mol.getName(), ex);
                }
            }
            
            // check to see if the threads are still running
            for (int i = 0; i < threads.length; ++i) {
                if (threads[i] == null 
                    || threads[i].isDone() 
                    || threads[i].isCancelled()) {
                    
                    // start a new thread
                    threads[i] = threadPool.submit
                        (new FragmentWorker ("Fragmenter-"+(i+1)));
                }
            }
            
            // block if the queue becomes too large (i.e., fragmentation falls
            //  behind)
            try {
                stagingQueue.put(mol);
            }
            catch (InterruptedException ex) {
                logger.log(Level.WARNING, "Staging molecule for fragmentation "
                           +"interrupted", ex);
            }
        }
        else {
            // singletons
            addMolecule (mol);
        }
    }

    /*
      synchronized protected void updateFragments 
      (Molecule mol, Collection<Molecule> frags) {
      int pos = molvec.size();
      for (Molecule frag : frags) {
      frag.aromatize();
      String fsmi = frag.toFormat("cxsmarts:q");
      BitSet member = fragments.get(fsmi);
      if (member == null) {
      fragments.put(fsmi, member = new BitSet ());
      try{
      Integer.parseInt(frag.getName());
      fragmentID.put(frag.getName(),fsmi);
      }catch(Exception e){
                        
      }
                
      }
      member.set(pos);
      }
      molvec.add(mol);
      }*/
    //ConcurrentHashMap<String, String> fragmentID2 = new ConcurrentHashMap<String,String>();

        
    synchronized protected void updateFragmentsSmiles 
        (Molecule mol, Collection<String> frags) {
        //int pos = molvec.size();
        int pos = getMoleculeCount();
        for (String frag : frags) {
            //frag.aromatize();
            String fsmi = frag.split("\t")[0];
            String fID = frag.split("\t")[1];
            Fragment f = fragments.get(fID);
            if (f == null) {
                f = new Fragment(fsmi,fID,new BitSet());
                fragments.put(fID, f);
            }
            f.getBitSet().set(pos);
        }
        //molvec.add(mol);
        try {
            //addMolecule(mol);
            cache.put(pos, mol);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    synchronized protected void addMolecule(Molecule m) /*throws Exception*/{
        molCompressVec.add(m.toFormat("cssdf"));
    }
        
    public int getFragmentCount () { return fragments.size(); }
    public Enumeration<Molecule> getMolecules () {
        Vector<Molecule> mols = new Vector<Molecule>();
        for (int i = 0;i<getMoleculeCount();i++) {
            mols.add(getMolecule(i));
        }
        return mols.elements();
    }
    public Enumeration<Molecule> getSingletons () {
        Vector<Molecule> mols = new Vector<Molecule>();
        for (int i = singletons.nextSetBit(0); 
             i >= 0; i = singletons.nextSetBit(i+1)) {
            mols.add(getMolecule(i));
        }
        return mols.elements();
    }
    
    private Molecule getMoleculeCache (int pos){
        return cache.get(pos);
    }
    private synchronized void putMoleculeCache(int pos, Molecule m){
        
        cache.put(pos,m);
        cachePuts++;
        if(cachePuts>cacheMax){
            cachePuts-=cacheTrims;
            int[] toRemove = new int[cacheTrims];
            int i=0;
            for(int r:cache.keySet()){
                toRemove[i++]=r;
                if(i>cacheTrims)break;
            }
            for(int r:toRemove){
                cache.remove(r);
            }                    
            System.out.println("cache flush");
        }
        
    }
    public synchronized Molecule getMolecule (int pos) {
        Molecule m = getMoleculeCache(pos);
        
        if(m!=null){
            return m;
        }else{
            try{
                m =MolImporter.importMol(molCompressVec.get(pos));
                putMoleculeCache(pos,m);
                return m;
            }catch(Exception ex){
                ex.printStackTrace();
                return null;
            }
        }
    }
    public int getMoleculeCount () { 
        //return molCompressVec.size();
        return cache.size();
    }
    public int getSingletonCount () { return singletons.cardinality(); }

    public int getScaffoldCount () { return scaffolds.size(); }
    public RGroupTable getRGroupTable (int scaffold) { 
        //System.out.println("Getting Rtable row ..." + scaffold);
        return scaffolds.get(scaffold);
    }
    public RGroupTable[] getRGroupTables () {
        return scaffolds.toArray(new RGroupTable[0]);
    }

    /** 
     * DataSeq interface
     */
    public int size () { return getMoleculeCount(); }
    public Molecule get (int index) { return getMolecule(index); }
    // return all scaffolds that contain the specified member index
    public RGroupTable[] getRGroupTablesFor (int index) {
        List<RGroupTable> rgs = new ArrayList<RGroupTable>();
        for (RGroupTable rg : scaffolds) {
            for (int i = 0; i < rg.rows.length; ++i) {
                if (rg.rows[i] == index) {
                    rgs.add(rg);
                    break;
                }
            }
        }
        return rgs.toArray(new RGroupTable[0]);
    }

    public int getRemainingUnprocessedMoleculeCount(){
        int l=0;
        for (int i = 0; i < threads.length; ++i) {
            l+=stagingQueue.size();
        }
        return l;
    }
    private int oldVal=-1; 
    public synchronized void updateUnprocessedProgress(){
        int pct = 100-100*getRemainingUnprocessedMoleculeCount()/remainingMoleculesForProcessing;
        propChange.firePropertyChange("progress", oldVal, pct);
        oldVal = pct;
    }
    public void closeFragmentationThreads(){
        if(fullyLoaded)return;
        try {
                
            for (int i = 0; i < threads.length; ++i) {
                stagingQueue.put(DONE);
            }
            fullyLoaded=true;
            remainingMoleculesForProcessing=getRemainingUnprocessedMoleculeCount();

            for (int i = 0; i < threads.length; ++i) {
                // block until this thread finishes
                logger.info("Waiting for thread to finish..."+(i+1));
                threads[i].get();
            }
        }
        catch (Exception ex) {
            logger.log(Level.WARNING, "Thread interrupted while waiting "
                       +"for fragmentation", ex);
            ex.printStackTrace();
        }
    }
    public void closeRGroupThreads(){
        try {
            for (int i = 0; i < threads.length; ++i) {
                fragQueue.put(DONEFrag);
            }
            for (int i = 0; i < threads.length; ++i) {
                // block until this thread finishes
                logger.info("Waiting for thread to finish..."+(i+1));
                threads[i].get();
            }
        }
        catch (Exception ex) {
            logger.log(Level.WARNING, "Thread interrupted while waiting "
                       +"for Rgroup Decomposition", ex);
            ex.printStackTrace();
        }
    }
    public void closeVariableScaffoldThreads(){
        try {
            for(int i=0;i<this.threads.length;i++){
                fragVariableQueue.put(DONEFragList);
            }
            for (int i = 0; i < threads.length; ++i) {
                logger.info("Waiting for thread to finish..."+(i+1));
                threads[i].get();
            }
        }
        catch (Exception ex) {
            logger.log(Level.WARNING, "Thread interrupted while waiting "
                       +"for Rgroup Decomposition", ex);
            ex.printStackTrace();
        }
    }
    
    private void processFragmentList(List<Fragment> l, int numAtoms){
        
        for (int i = 0; i < threads.length; ++i) {
            if (threads[i] == null 
                || threads[i].isDone() 
                || threads[i].isCancelled()) {

                // start a new thread
                threads[i] = threadPool.submit
                    (new VariableScaffoldWorker("VariableFragmentWorker-"+(i+1),numAtoms));
            }
        }
        try {
            fragVariableQueue.put(l);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }       
    }
    private static Molecule getTopological(Molecule m, boolean atoms, boolean bonds){
        Molecule r = m.cloneMolecule();
        r.hydrogenize(false);
        r.expandSgroups();
        if(atoms){
            for (MolAtom a : r.getAtomArray()) {
                a.setAtno(6);
                a.setRadical(0);
                a.setCharge(0);
                a.setFlags(0);
                a.setMassno(0);
            }
        }

        /* may or may not be a good idea */
        if(bonds){
            for(MolBond b : r.getBondArray()) {
                b.setFlags(0);
                b.setType(1);
            }
        }
                
        return r;
    }
    private static Molecule getTopological(String smiles, boolean atoms, boolean bonds) throws MolFormatException{
        MolHandler mh = new MolHandler();
        mh.setMolecule(smiles);
        return getTopological(mh.getMolecule(),atoms,bonds);
        
    }
    private void findVariableScaffolds(int diffAtoms){
        if(diffAtoms==0)return;
                
        //This is a clustering of the fragments based on a Topological invariant.
        HashMap<String,List<Fragment>> topoHashMap = new HashMap<String,List<Fragment>>();
        for (String ID : fragments.keySet()) {
            Fragment f = fragments.get(ID);
            if(!f.isPreLoaded() && !f.isGeneric()){
                String topHash=ID.split("-")[0];
                //System.out.println(topHash);
                List<Fragment> members=topoHashMap.get(topHash);
                if(members==null){
                    members=new ArrayList<Fragment>();
                    topoHashMap.put(topHash,members);
                }
                members.add(f);                         
            }
        }
                
        int k=0;
        double pctDone=0;
        for(String topHash:topoHashMap.keySet()){
                        
            List<Fragment> l = topoHashMap.get(topHash);
                        
            if(l.size()>1){
                HashMap<String,List<Fragment>> stopoHashMap = new HashMap<String,List<Fragment>>();
                //System.out.println("test It");
                for(Fragment f: l){
                    String stopHash=f.getTopoForm();
                    List<Fragment> members=stopoHashMap.get(stopHash);
                    if(members==null){
                        members=new ArrayList<Fragment>();
                        stopoHashMap.put(stopHash,members);
                    }
                    members.add(f);     
                }
                for(String stopHash:stopoHashMap.keySet()){
                    List<Fragment> m = stopoHashMap.get(stopHash);
                    if(m.size()>1){
                        //System.out.println(topHash + "\t" + l.size());
                        processFragmentList(m,diffAtoms);
                        int pct = (int) (100. * (k + 1) / topoHashMap.size() + 0.5);
                        propChange.firePropertyChange("progress", pctDone, pct);
                        pctDone = pct;
                    }
                }
                                
            }
            k++;
                        
        }
                
                
        closeVariableScaffoldThreads();
                
    }
    private void findVariableScaffolds(List<Fragment> l,int diffAtoms){
        MolHandler mh = new MolHandler();
        MolSearch ms = new MolSearch();
        ms.setSubgraphSearch(true);
        int[][] amorph;
        try {
            Set<String> columns = new HashSet<String>();
            Fragment rep = l.get(0);
            String topHash= rep.getTopoForm();
                        
            mh.setMolecule(rep.getSmarts());
            Molecule mrep = getTopological(mh.getMolecule(),true,false);                        
            if(debug)System.out.println("START:" + mrep.exportToFormat("cxsmarts:q"));
            VFLib2 vf = VFLib2.automorphism(mrep);
            amorph = vf.findAll();
                        
            if(debug){
                try {
                    System.out.println();
                    System.out.println(topHash+ "\t"+ "Automorph:" + amorph.length);
                    System.out.println(topHash+ "\t"+ mrep.exportToFormat("smiles"));
                } catch (MolExportException e1) {
                    // TODO Auto-generated catch block
                    e1.printStackTrace();
                }
            }
            ms.setQuery(mrep);
            HashMap<Integer,List<String>> atomLists = new HashMap<Integer,List<String>>(mrep.getAtomCount());
            List<MolAtom[]> dumbThingAboutBonds = new ArrayList<MolAtom[]>();
            List<Integer> dumbThingAboutBonds2 = new ArrayList<Integer>();
            List<Integer> dumbThingAboutBonds3 = new ArrayList<Integer>();
            List<Fragment> listOriginator = new ArrayList<Fragment>(); 
            for(int i=0;i<mrep.getAtomCount();i++){
                atomLists.put(i,new ArrayList<String>());
            }
            for(MolBond mb:mrep.getBondArray()){
                MolAtom[] maarr = new MolAtom[]{mb.getAtom1(),mb.getAtom2()};
                dumbThingAboutBonds.add(maarr);
                dumbThingAboutBonds2.add(mb.getType());
                dumbThingAboutBonds3.add(mb.getFlags());
            }
                        
            for (Fragment f : l) {
                try {
                    mh.setMolecule(f._can_smiles);
                    Molecule m = mh.getMolecule();
                    Molecule minv = getTopological(m,true,false);
                                        
                    ms.setTarget(minv);
                                        
                    int[][] match = ms.findAll();
                    MolAtom[] maArr = m.getAtomArray();
                    if(match!=null){
                        for(int[] amorphi:amorph){
                            for(int[] matches:match){
                                String[] aCol = new String[maArr.length];
                                String uniq="";
                                for (int i = 0; i < maArr.length; i++) {
                                    aCol[amorphi[i]]=maArr[matches[i]].getSymbol();
                                    //"#"+maArr[matches[i]].getAtno();
                                }
                                for (int i = 0; i < maArr.length; i++) {
                                    uniq+=aCol[i];
                                }
                                if(!columns.contains(uniq)){
                                    columns.add(uniq);
                                    for(int i=0;i<maArr.length;i++){
                                        atomLists.get(i).add(aCol[i]);
                                    }
                                    listOriginator.add(f);
                                }
                            }
                        }
                        //System.out.println(topHash + "\t" + f._can_smiles);
                                                
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            //Remove all columns that are constant
            //And set all common atoms to their correct values
            MolAtom[] marr=mrep.getAtomArray();
            for (int i = 0; i < mrep.getAtomCount(); i++) {
                /*
                Set<String> atomSet=new HashSet<String> (atomLists.get(i));
                marr[i].setSMARTS(atomSet.toString());
                */
                List<String> atoms = atomLists.get(i);
                BitSet unique = new BitSet ();
                for (String a : atoms) {
                    unique.set(MolAtom.numOf(a));
                }

                if(unique.cardinality() <= 1) {
                    atomLists.remove(i);
                    marr[i].setAtno(unique.nextSetBit(0));
                }
                else {
                    marr[i].setAtno(MolAtom.LIST);
                    int[] list = new int[unique.cardinality()];
                    for (int j = unique.nextSetBit(0), k = 0; j>=0;
                         j = unique.nextSetBit(j+1)) {
                        list[k++] = j;
                    }
                    marr[i].setList(list);
                }
            }
                        
            //Make per-fragment instance atom list Matrix
            List<ReducedFragmentSequence> sequencesQuery = new ArrayList<ReducedFragmentSequence>();
            int[] index_to_atom_index=new int[atomLists.size()];
            int k=0;
            for(int j=0;j<listOriginator.size();j++){
                ReducedFragmentSequence rfs = new ReducedFragmentSequence();
                rfs.addParentFragments(listOriginator.get(j));
                sequencesQuery.add(rfs);
            }
            for (int i:atomLists.keySet()) {
                List<String> aSymbols = atomLists.get(i);
                for(int j=0;j<listOriginator.size();j++){
                    sequencesQuery.get(j).addAtomSet(k,aSymbols.get(j));
                }
                index_to_atom_index[k]=i;
                k++;
            }
            List<ReducedFragmentSequence> sequenceTargets = sequencesQuery;
                        
            /*
             * Set Target list to Query List.
             * 
             * For each member of Query fragment sequence list, compare to Target members of fragment sequence list.
             *          For each pair, if the atom disagreement is less than the 1:
             *                  Generate the consensus string, with "X" for the disagreements.
             *                  Generate a Map from that string, to a list of sets of atoms.
             *                  e.g.
             *                          CCCN
             *                          CCCC
             * 
             *                          CCCX ->  0:{C}
             *                                           1:{C}
             *                                           2:{C}
             *                                           3:{C,N}
             *                  If that map already exists, combine the lists of sets of atoms.
             *                  e.g.
             *                          CCCO
             *                          CCCC
             * 
             *                          CCCX ->  0:{C}
             *                                           1:{C}
             *                                           2:{C}
             *                                           3:{C,O,N}
             * 
             *          For each consensus variable sequence, generate a molecule, with the 
             *          template atoms replaced by the set, if the smarts is unique:
             *                  Add it as a new fragment.
             *                  Have the bitset inherit all members from the parent fragments.
             * 
             *  Replace Query list with the consensus values generated. Repeat n times.
             * 
             * 
             */
                        
                        
            for (int u = 0; u < diffAtoms; u++) {
                Map<String, ReducedFragmentSequence> newSequences = new HashMap<String, ReducedFragmentSequence>();
                Set<String> alreadyCompared = new HashSet<String>();
                for (ReducedFragmentSequence parent : sequencesQuery) {
                    for (ReducedFragmentSequence comp : sequenceTargets) {
                        if (!alreadyCompared.contains(parent
                                                      .getAtomSequence()
                                                      + "_"
                                                      + comp.getAtomSequence())) {
                            ReducedFragmentSequence common = parent
                                .getCommonFragmentSequence(comp);
                            if (common != null) {
                                ReducedFragmentSequence prevc = newSequences
                                    .get(common.getAtomSequence());
                                if (prevc != null) {
                                    prevc.addParentsAndAtoms(common);
                                } else {
                                    newSequences.put(
                                                     common.getAtomSequence(),
                                                     common);
                                }
                            }
                            alreadyCompared.add(comp.getAtomSequence()
                                                + "_" + parent.getAtomSequence());
                        }
                    }
                }
                sequencesQuery = new ArrayList<ReducedFragmentSequence>(newSequences.values());

                Map<String, Set<Fragment>> newScaffolds = new HashMap<String, Set<Fragment>>();
                Map<String, String> newScaffoldsRep = new HashMap<String, String>();
                for (String aDiff : newSequences.keySet()) {
                    ReducedFragmentSequence newFrag = newSequences
                        .get(aDiff);
                    int j=0;
                    for (int i = 0; i < newFrag.getVariableAtomCount(); i++) {
                        if (newFrag.getAtomSet(i).size() > 1) {
                            marr[index_to_atom_index[i]]
                                .setAliasstr(((char)('X' + j++))+"");
                           // marr[index_to_atom_index[i]].setSMARTS(newFrag.getAtomSet(i).toString());
                            marr[index_to_atom_index[i]].setAtno(MolAtom.LIST);
                            int[] list = new int[newFrag.getAtomSet(i).size()];
                            int n=0;
                            for (String atSym:newFrag.getAtomSet(i)) {
                                list[n++] = MolAtom.numOf(atSym);
                            }
                            marr[index_to_atom_index[i]].setList(list);
                            marr[index_to_atom_index[i]].setQueryAromaticity(0);
                        } else {
                            marr[index_to_atom_index[i]]
                                .setAliasstr("");
                            marr[index_to_atom_index[i]].setAtno(MolAtom.numOf(newFrag.getAtomSet(i).iterator().next()));
                            marr[index_to_atom_index[i]].setQueryAromaticity(0);
                        }
                    }
                    int inc=0;                  
                    String smarts = mrep.toFormat("cxsmarts:q");
                    Set<Fragment> parentFragList = newScaffolds.get(smarts);
                    if (parentFragList == null) {
                        parentFragList = new HashSet<Fragment>();
                        newScaffolds.put(smarts, parentFragList);
                        newScaffoldsRep.put(smarts, mrep.toFormat("sdf"));
                    }
                    parentFragList.addAll(newFrag.parentFragments);
                }
                for (String scaffSmiles : newScaffolds.keySet()) {
                    Set<Fragment> parentFragList = newScaffolds
                        .get(scaffSmiles);
                    BitSet bs = new BitSet();
                    bs.clear();
                    for (Fragment parentFrag : parentFragList) {
                        bs.or(parentFrag.getBitSet());
                    }
                    Fragment newFrag =new Fragment(scaffSmiles,newScaffoldsRep.get(scaffSmiles),bs,null,true); 
                    //Do not add a new fragment if y
                    //addScaffold(scaffSmiles,newScaffoldsRep.get(scaffSmiles),null,bs,true);
                    //newFrag._bs = bs;
                }
                if (debug) {
                    System.out.println();
                    for (Fragment f : l) {
                        System.out.print(f.getSmarts() + ".");
                    }
                    for (String s : newScaffolds.keySet()) {
                        System.out.print(s + ".");
                    }
                    System.out.println();
                    for (int i : atomLists.keySet()) {
                        System.out.println(topHash + "\t"
                                           + atomLists.get(i));
                    }
                    System.out.print(topHash + "\t"
                                     + mrep.toFormat("cxsmarts"));
                }
            }
        } catch (Exception e2) {
            e2.printStackTrace();
        }
    }
    public static class ReducedFragmentSequence{
        Set<Fragment> parentFragments=null;
        Map<Integer,Set<String>> allowedAtoms = new HashMap<Integer,Set<String>>();
        String atomSequence=null;
        
        public String getAtomVar(int pos){
            Set<String> atoms = allowedAtoms.get(pos);
            if(atoms!=null){
                if(atoms.size()>1){
                    return "X";
                }else{
                    return (new ArrayList<String>(atoms)).get(0);
                }
            }
            return null;
        }
        public int getVariableAtomCount(){
            return allowedAtoms.size();
        }
        public Set<String> getAtomSet(int pos){
            return allowedAtoms.get(pos);
        }
        public void addAtomSet(int pos,String atomSet){
            addAtomSet(pos,new HashSet<String>(Arrays.asList(atomSet)));
        }
        public void addAtomSet(int pos,Set<String> atomSet){
                
            Set<String> catomSet = allowedAtoms.get(pos);
            if(catomSet==null){
                catomSet = new HashSet<String>();
                allowedAtoms.put(pos, catomSet);
            }else{
                if(catomSet.size()==1){
                    catomSet.addAll(atomSet);
                    if(catomSet.size()!=1){
                        atomSequence=null;
                    }
                    return;
                }
            }
            catomSet.addAll(atomSet);
        }
        public void addParentsAndAtoms(ReducedFragmentSequence rfs){
            addParentFragments(rfs.parentFragments);
            for(int i=0;i<allowedAtoms.size();i++){
                addAtomSet(i,rfs.getAtomSet(i));
            }
        }
        public void addParentFragments(Set<Fragment> parents){
            if(parentFragments==null){
                parentFragments=new HashSet<Fragment>();
            }
            parentFragments.addAll(parents);
        }
        public void addParentFragments(Fragment parent){
            if(parentFragments==null){
                parentFragments=new HashSet<Fragment>();
            }
            parentFragments.add(parent);
        }
        public String getAtomSequence(){
            if(atomSequence==null){
                StringBuilder sb = new StringBuilder();
                for(int i=0;i<allowedAtoms.size();i++){
                    sb.append(getAtomVar(i));
                }
                atomSequence=sb.toString();
            }
            return atomSequence;
        }
        public ReducedFragmentSequence getCommonFragmentSequence(ReducedFragmentSequence rfs){
            int cutoff =1;
            if(parentFragments.containsAll(rfs.parentFragments))return null;
            int diff = 0;
            ReducedFragmentSequence commonFrag = new ReducedFragmentSequence();
                        
            //List<Set<String>> allowableAtoms = new ArrayList<Set<String>>();
            for (int k = 0; k < allowedAtoms.size(); k++) {
                if(!this.getAtomVar(k).equals(rfs.getAtomVar(k))) {
                    if(!this.getAtomVar(k).equals("X") && !rfs.getAtomVar(k).equals("X")){
                        diff++;
                    }
                }
                if (diff > cutoff) {
                    break;
                }
                commonFrag.addAtomSet(k, this.getAtomSet(k));
                commonFrag.addAtomSet(k, rfs.getAtomSet(k));
            }
            if (diff <= cutoff) {
                commonFrag.addParentFragments(this.parentFragments);
                commonFrag.addParentFragments(rfs.parentFragments);
                return commonFrag;
            }else{
                return null;
            }
        }
        
    }
    
    
    /**
     * This does the R-group generation
     * 
     */
    public void run() {
        propChange.firePropertyChange("status", null, "Preprocessing ... ");
        closeFragmentationThreads();
        propChange.firePropertyChange("status", null, "Looking for generic structures ... ");
        findVariableScaffolds(this.maxVariableAtoms);
        propChange.firePropertyChange("status", null, "Generating R-groups ... ");
                
        MolHandler mh = new MolHandler();
        int pctDone = 0;
        singletons.clear();
        for (int i = 0; i < getMoleculeCount(); ++i) {
            singletons.set(i);
        }
        scaffolds.clear();
                
        mh.setQueryMode(true);

        int fragsize = fragments.size();
        if (fragsize > 0) {
            for (int i = 0; i < threads.length; ++i) {
                // start a new thread
                threads[i] = threadPool.submit
                    (new RGroupWorker ("RGroupWorker-"+(i+1)));
            }
        }
        
        int k = 0;
        Set<String> remove = new HashSet<String>();
        for (String ID : fragments.keySet()) {
            Fragment f = fragments.get(ID);
            BitSet mb = f.getBitSet();
            if (((mb!=null)?mb.cardinality():1000) < minFragMemberCount) {
                if(debug){
                    System.out.println("Removing scaffold:"+f._can_smiles + "\t" + f._id);
                }
                //fragments.remove(f.getID(), f);
                remove.add(f.getID());
            } else {
                try {
                    //logger.info("Queuing fragment "+ID+"...");
                    fragQueue.put(f);
                }
                catch (InterruptedException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                } 
            }
            int pct = (int) (100. * (k++ + 1) / fragsize + 0.5);
            propChange.firePropertyChange("progress", pctDone, pct);
            pctDone = pct;
        }

        for (String id: remove)
            fragments.remove(id);
        
        if (debug) {
            if (!singletons.isEmpty()) {
                System.out.print("** singletons");
                for (int i = singletons.nextSetBit(0); i >= 0; i = singletons
                         .nextSetBit(i + 1)) {
                    System.out.print(this.getMolecule(i).getName());
                }
                System.out.println();
            }
        }

        /*
        try {
            // now try extend singletons
            propChange.firePropertyChange
                ("status", null, "Extending singletons...");
            Fragment f = addScaffold (new MolHandler ("c1ccccc1", true)
                                      .getMolecule(), singletons, false);
            
            fragQueue.put(f);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        */
        closeRGroupThreads();
        
        propChange.firePropertyChange("progress", pctDone, 100);
        if(markushlike){
            propChange.firePropertyChange("status", null, "Generating variable attachments ... ");
            k=0;
            for(RGroupTable rt:scaffolds){
                RGroupSolver.markushAttach(rt);
                int pct = (int) (100. * (k++ + 1) / scaffolds.size() + 0.5);
                propChange.firePropertyChange("progress", pctDone, pct);
                pctDone = pct;
            }
            propChange.firePropertyChange("progress", pctDone, 100);
        }
    }

    public RGroupTable generate (String scaffold) {
        try {
            MolHandler mh = new MolHandler (scaffold, true);
            return generate (mh.getMolecule());
        }
        catch (MolFormatException ex) {
            throw new IllegalArgumentException
                ("Bogus scaffold specified: " + scaffold);
        }
    }

    public RGroupTable generate (Fragment f) {
        try {
            MolHandler mh = new MolHandler (f.getSmarts(), true);
            Molecule m = mh.getMolecule();
            return generate (m, f.getBitSet(), true,f.getPreloadID(),(false)&f.isGeneric());
        } catch (MolFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;
    }
    public RGroupTable generate (Molecule scaffold) {
        return generate (scaffold, null, true,-1,false);
    }

    public RGroupTable generate (Molecule scaffold, boolean prune) {
        return generate (scaffold, null, prune,-1,false);
    }
    public RGroupTable generate (Molecule scaffold, BitSet subset) {
        return generate (scaffold, subset, true,-1,false);
    }
    public RGroupTable generate (Molecule scaffold, 
                                 BitSet subset, 
                                 boolean prune, int scaffNum) {
        return generate (scaffold, subset, prune,scaffNum,false);
    }
    public RGroupTable generate (Molecule scaffold, 
                                 BitSet subset, 
                                 boolean prune, int scaffNum, boolean force) {
        RGroupTable rtab = null;
        try {
            rtab = (scaffNum==-1)?
                RGroupSolver.solve(scaffold, this, subset, force)
                : RGroupSolver.solveLoad(scaffNum, scaffold, this, subset);
            
            if(scaffNum!=-1)
                prune = false;
            
            if (rtab != null)
                rtab = addScaffold (rtab,prune);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return rtab;
    }

    public RGroupTable extendScaffold (RGroupTable tab) {       
        RGroupTable xtab = RGroupSolver.extend(tab);
        if (xtab != null) {
            scaffolds.add(xtab);
        }
        return xtab;
    }

    public void sort (final Molecule ref) {
        // sort the scaffolds using ref as the anchor
        final MolFpFactory fpFact = MolFpFactory.getInstance();
        Collections.sort(scaffolds, new Comparator<RGroupTable>() {
                             public int compare 
                                 (RGroupTable a, RGroupTable b) {
                                 double simA = fpFact.tanimotoSim
                                     (ref, a.getScaffold());
                                 double simB = fpFact.tanimotoSim
                                     (ref, b.getScaffold());
                                 int d = 0;
                                 if (simA > simB) d = -1;
                                 else if (simA < simB) d = 1;

                                 if (d == 0) {
                                     d = b.getRowCount() - a.getRowCount();
                                 }
                                 if (d == 0) {
                                     d = b.getScaffoldComplexity() 
                                         - a.getScaffoldComplexity();
                                 }
                                 return d;
                             }
                         });
        propChange.firePropertyChange("reordering", ref, ref);
    }

    public void addPropertyChangeListener (PropertyChangeListener l) {
        propChange.addPropertyChangeListener(l);
    }
    public void removePropertyChangeListener (PropertyChangeListener l) {
        propChange.removePropertyChangeListener(l);
    }
    public PropertyChangeSupport getPropertyChangeSupport () { 
        return propChange; 
    }

    public static void main (String[] argv) throws Exception {
        
        
        RGroupGenerator rgroup = new RGroupGenerator ();
        for (int i = 0; i < argv.length; ++i) {
            MolImporter molimp = new MolImporter (argv[i]);
            for (Molecule mol; (mol = molimp.read()) != null; ) {
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

        // print each scaffold and its members
        for (int i = 0; i < rgroup.getScaffoldCount(); ++i) {
            RGroupTable rtab = rgroup.getRGroupTable(i);
            System.out.print(rtab.getScaffold().toFormat("cxsmarts"));
            for (int r= 0; r < rtab.getRowCount(); ++r) {
                String id = (String)rtab.getValueAt(r, 0); 
                Molecule m = (Molecule)rtab.getValueAt(r, 1);
                System.out.print(" " + id + ":");
                // now the rest of the columns are r-group ligands
                for (int g = 0; g < rtab.getRGroupCount(); ++g) {
                    String rg = rtab.getRGroupLabel(g);
                    Molecule ligand = rtab.getRGroup(r, g);
                    if (ligand != null) {
                        System.out.print
                            ("<"+rg+">"+ligand.toFormat("cxsmarts")
                             +"</"+rg+">");
                    }
                }
            }
            System.out.println();
        }

        rgroup.shutdown();
    }
}

