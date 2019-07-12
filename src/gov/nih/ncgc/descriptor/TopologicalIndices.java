// $Id: TopologicalIndices.java 4028 2010-02-01 21:32:31Z nguyenda $

package gov.nih.ncgc.descriptor;

import java.util.Vector;
import java.util.Arrays;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;
import chemaxon.struc.MolBond;
import chemaxon.struc.MolAtom;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.plugin.PluginException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.calculations.ChargePlugin;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Algebra;


public class TopologicalIndices {
    static final double ONE_OVER_LOG2 = 1.44269504088896340737;

    static double log2 (double x) {
	return Math.log(x)*ONE_OVER_LOG2;
    }

    // table lookup of hybridization & ratios r / r_C(sp3)
    static double[][] KappaTable = new double[255][5];
    static {
	KappaTable[6][MolAtom.HS_SP2] = -.13;
	KappaTable[6][MolAtom.HS_SP]  = -.22;
	KappaTable[7][MolAtom.HS_SP3] = -.04;
	KappaTable[7][MolAtom.HS_SP2] = -.20;
	KappaTable[7][MolAtom.HS_SP]  = -.29;
	KappaTable[8][MolAtom.HS_SP3] = -.04;
	KappaTable[8][MolAtom.HS_SP2] = -.20;
	KappaTable[15][MolAtom.HS_SP3] = .43;
	KappaTable[15][MolAtom.HS_SP2] = .30;
	KappaTable[16][MolAtom.HS_SP3] = .35;
	KappaTable[16][MolAtom.HS_SP2] = .22;
	for (int i = 0; i < 5; ++i) {
	    KappaTable[9][i] = -.07; // F
	    KappaTable[17][i] = .29; // Cl
	    KappaTable[35][i] = .48; // Br
	    KappaTable[53][i] = .73; // I
	}
    }

    static String[] _indices = null;
    public static String[] getNames () {
	if (_indices == null) {
	    try {
		Method[] methods = TopologicalIndices.class.getMethods();
		Vector<String> names = new Vector<String>();
		for (Method m : methods ) {
		    if (m.getName().startsWith("TI")) {
			names.add(m.getName());
		    }
		}
		_indices = names.toArray(new String[0]);
		Arrays.sort(_indices);
	    }
	    catch (Exception ex) {
		ex.printStackTrace();
	    }
	}
	return _indices;
    }

    private double[][] _adjmat; // adjacency matrix
    private double[][] _dismat; // distance matrix
    private double[][] _lapmat; // laplace matrix
    private double[][] _burden; // burden matrix
    private double[][] _charge; // charge matrix
    private double[][] _dist3d; // distance matrix 3D

    private MatrixSolver _solverD; // solver for distance matrix
    private MatrixSolver _solverA; // solver for adjacency matrix
    private MatrixSolver _solverL; // solver for laplacian matrix
    private MatrixSolver _solverB; // solver for burden matrix
    private MatrixSolver _solverC; // solver for charge matrix

    private Molecule mgraph;

    static class MatrixSolver {
	DoubleMatrix2D M;
	EigenvalueDecomposition _evd;
	Algebra alg = new Algebra ();

	public MatrixSolver (double[][] matrix) {
	    M = DoubleFactory2D.dense.make(matrix);
	}

	protected synchronized EigenvalueDecomposition getEvd () {
	    if (_evd == null) {
		_evd = new EigenvalueDecomposition (M);
	    }
	    return _evd;
	}

	public double getFiedlerValue () {
	    DoubleMatrix1D eval = getEvd().getRealEigenvalues();
	    // return the first non-zero eigenvalue starting at the
	    //  second eigenvalue.  if the matrix contains disconnected
	    //  components, then there will as many zero eigenvalues.
	    for (int i = 1; i < eval.size(); ++i) {
		double x = eval.getQuick(i);
		if (x > 1e-6) {
		    return x;
		}
	    }
	    return 0.;
	}

	public double[] getFiedlerVector () {
	    DoubleMatrix1D eval = getEvd().getRealEigenvalues();
	    // return the first non-zero eigenvalue starting at the
	    //  second eigenvalue.  if the matrix contains disconnected
	    //  components, then there will as many zero eigenvalues.
	    for (int i = 1; i < eval.size(); ++i) {
		double x = eval.getQuick(i);
		if (x > 1e-6) {
		    return getEvd().getV().viewColumn(i).toArray();
		}
	    }
	    return new double[]{};
	}

	public double[] eigenvalues () { 
	    return getEvd().getRealEigenvalues().toArray(); 
	}
	public double eigenvalue (int which) {
	    return getEvd().getRealEigenvalues().getQuick(which);
	}
	public double[][] eigenvectors () { 
	    return getEvd().getV().toArray();
	}
	public double[] eigenvector (int which) {
	    return getEvd().getV().viewColumn(which).toArray();
	}

	// solve Az = v, return vector z
	public double[] solve (double[] a) {
	    if (a.length != M.rows()) {
		throw new IllegalArgumentException 
		    ("Input vector has wrong dimension!");
	    }

	    DoubleMatrix2D V = DoubleFactory2D.dense.make(a.length, 1);
	    DoubleMatrix2D A = M.copy();
	    for (int i = 0; i < a.length; ++i) {
		double s = 0.;
		for (int j = 0; j < a.length; ++j) {
		    s += A.getQuick(i, j);
		}
		V.setQuick(i, 0, s);
		A.setQuick(i, i, a[i]);
	    }

	    DoubleMatrix2D z = null;
	    try {
		z = alg.solve(A, V);
	    }
	    catch (Exception ex) {
		ex.printStackTrace();
		System.err.println("A = " + A);
		System.err.println("V = " + V);
	    }
	    return z != null ? z.viewColumn(0).toArray() : new double[]{};
	}

	public double get (int r, int c) {
	    return M.getQuick(r, c);
	}
    }

    public TopologicalIndices (Molecule mgraph) {
	if (mgraph == null) {
	    throw new IllegalArgumentException 
		("Input molecule graph is null");
	}
	// make sure we have a hydrogen depleted graph
	mgraph.implicitizeHydrogens(MolAtom.ALL_H);
	mgraph.aromatize();
	mgraph.calcHybridization();
	this.mgraph = mgraph;
    }

    public synchronized double[][] getAdjacencyMatrix () {
	if (_adjmat == null) {
	    int[][] ctab = mgraph.getCtab();
	    _adjmat = new double[ctab.length][ctab.length];
	    for (int i = 0; i < ctab.length; ++i) {
		for (int j = 0; j < ctab[i].length; ++j) {
		    _adjmat[i][ctab[i][j]] = 1.;
		}
	    }
	}
	return _adjmat;
    }

    public synchronized double[][] getDistanceMatrix () {
	if (_dismat == null) {
	    int atomCount = mgraph.getAtomCount();
	    _dismat = new double[atomCount][atomCount];
	    int[][] tab = mgraph.getBtab();

	    for (int i = 0; i < atomCount; ++i) {
		_dismat[i][i] = 0.;
		for (int j = i+1; j < atomCount; ++j) {
		    _dismat[j][i] = _dismat[i][j] 
			= tab[i][j] < 0 ? atomCount : 1.;
		}
	    }
	    tab = null;

	    /* Floyd's all-pairs shortest path algorithm */
	    for (int k = 0; k < atomCount; ++k)
		for (int i = 0; i < atomCount; ++i) 
		    for (int j = 0; j < atomCount; ++j)
			_dismat[i][j] = Math.min
			    (_dismat[i][j], _dismat[i][k] + _dismat[k][j]);
	    /*
	    System.out.println("distance matrix...");
	    for (int i = 0; i < atomCount; ++i) {
		for (int j = 0; j < atomCount; ++j) {
		    System.out.printf("%1$3.0f", _dismat[i][j]);
		}
		System.out.println();
	    }
	    */
	}
	return _dismat;
    }

    public synchronized double[][] getDistanceMatrix3d () {
	if (_dist3d == null) {
	    if (mgraph.getDim() < 3) {
		mgraph.clean(3, null);
	    }
	    double x, y, z;
	    MolAtom[] atoms = mgraph.getAtomArray();
	    _dist3d = new double[atoms.length][atoms.length];
	    for (int i = 0; i < atoms.length; ++i) {
		for (int j = i+1; j < atoms.length; ++j) {
		    x = atoms[i].getX() - atoms[j].getX();
		    y = atoms[i].getY() - atoms[j].getY();
		    z = atoms[i].getZ() - atoms[i].getZ();
		    _dist3d[i][j] = _dist3d[j][i] = Math.sqrt(x*x + y*y + z*z);
		}
	    }
	}
	return _dist3d;
    }

    public synchronized double[][] getBurdenMatrix () {
	if (_burden == null) {
	    int atomCount = mgraph.getAtomCount();

	    _burden = new double[atomCount][atomCount];
	    // setup default values..
	    for (int i = 0; i < atomCount; ++i) {
		for (int j = i+1; j < atomCount; ++j) {
		    _burden[i][j] = _burden[j][i] = 0.001;
		}
		_burden[i][i] = mgraph.getAtom(i).getAtno();
	    }

	    for (MolBond b : mgraph.getBondArray()) {
		MolAtom a1 = b.getAtom1();
		MolAtom a2 = b.getAtom2();
		int idx1 = mgraph.indexOf(a1);
		int idx2 = mgraph.indexOf(a2);

		switch (b.getType()) {
		case 1: 
		    _burden[idx1][idx2] = 0.1;
		    _burden[idx2][idx1] = 0.1;
		    break;
		case 2:
		    _burden[idx1][idx2] = 0.2;
		    _burden[idx2][idx1] = 0.2;
		    break;
		case 3:
		    _burden[idx1][idx2] = 0.3;
		    _burden[idx2][idx1] = 0.3;
		    break;
		case MolBond.AROMATIC:
		    _burden[idx1][idx2] = .15;
		    _burden[idx2][idx1] = .15;
		}

		if (a1.isTerminalAtom() || a2.isTerminalAtom()) {
		    _burden[idx1][idx2] += 0.01;
		    _burden[idx2][idx1] += 0.01;
		}
	    }
	}
	return _burden;
    }

    public synchronized double[][] getChargeMatrix () {
	if (_charge == null) {
	    double[][] A = getAdjacencyMatrix ();
	    double[][] D = getDistanceMatrix ();
	    double[][] M = new double[D.length][D.length];
	    double x, d;
	    for (int i = 0; i < D.length; ++i) {
		for (int j = 0; j < D.length; ++j) {
		    x = 0.;
		    for (int k = 0; k < D.length; ++k) {
			d = k != j ? (1./(D[k][j]*D[k][j])) : 0.;
			x += A[i][k] * d;
		    }
		    M[i][j] = x;
		}
	    }
	    _charge = M;
	}
	return _charge;
    }

    public static double[][] makeRawLaplacian (MoleculeGraph g) {
	int[][] tab = g.getBtab();
	int size = g.getAtomCount();
	double[][] laplace = new double[size][size];

	for (int i = 0; i < size; ++i) {
	    double sum = 0.;
	    for (int j = 0; j < i; ++j) {
		sum += laplace[i][j];
	    }

	    for (int j = i + 1; j < size; ++j) {
		double v = 0.;
		if (tab[i][j] >= 0) {
		    MolBond b = g.getBond(tab[i][j]);
		    v = b.getType();

		    sum += v;
		}
		laplace[i][j] = v;
		laplace[j][i] = v;
	    }
	    
	    laplace[i][i] = sum;
	}
	tab = null;

	for (int i = 0; i < size; ++i) {
	    for (int j = 0; j < size; ++j) {
		if (i != j) {
		    laplace[i][j] *= -1.;
		}
	    }
	}
	
	return laplace;
    }

    public static double[][] makeNormalizedLaplacian 
	(MoleculeGraph g, boolean signless) {

	int[][] tab = g.getBtab();
	int size = g.getAtomCount();

	double[][] laplace = new double[size][size];
	double[][] sim = new double[size][size]; // similarity matrix
	double vol[] = new double[size]; // volume

	for (int i = 0; i < vol.length; ++i) {
	    vol[i] = 0.;
	    for (int j = 0; j < i; ++j) {
		vol[i] += sim[i][j];
	    }

	    for (int j = i + 1; j < vol.length; ++j) {
		double v = 0.;
		if (tab[i][j] >= 0) {
		    //v = m.getBond(tab[i][j]).getType();
		    v = 1.;
		}

		sim[i][j] = v;
		sim[j][i] = v;

		vol[i] += v;
	    }

	    if (vol[i] > 0.) {
		vol[i] = 1./Math.sqrt(vol[i]);
	    }
	}
	tab = null;

	// construct a normalized laplacian matrix: I - L
	for (int i = 0; i < vol.length; ++i) {
	    for (int j = i + 1; j < vol.length; ++j) {
		double x = vol[i] * vol[j] * sim[i][j];
		if (x != 0.) {
		    double v = signless ? x : -x;
		    laplace[i][j] = v;
		    laplace[j][i] = v;
		}
	    }
	    laplace[i][i] = 1.;
	}

	vol = null;
	sim = null;

	return laplace;
    }

    public double[][] getLaplacianMatrix () {
	return getLaplacianMatrix (true);
    }

    public synchronized double[][] getLaplacianMatrix (boolean normalized) {
	if (_lapmat == null) {
	    _lapmat = normalized ? makeNormalizedLaplacian (mgraph, false)
		: makeRawLaplacian (mgraph);
	}
	return _lapmat;
    }

    public double TIbalabanJ () {
	return balabanJ (getDistanceMatrix ());
    }
    public double TIbalabanJ3d () {
	return balabanJ (getDistanceMatrix3d ());
    }

    public double balabanJ (double[][] D) {

	double J = 0.;
	for (MolBond e : mgraph.getBondArray()) {
	    double[] arow = D[mgraph.indexOf(e.getAtom1())];
	    double[] brow = D[mgraph.indexOf(e.getAtom2())];
	    double sa = 0., sb = 0.;
	    for (int i = 0; i < arow.length; ++i) {
		sa += arow[i];
		sb += brow[i];
	    }

	    if (sa > 0. && sb > 0.) {
		J += 1./Math.sqrt(sa * sb);
	    }
	}
	return J * mgraph.getBondCount() / (mgraph.getSSSR().length + 1);
    }

    public double wiener (double[][] D) {
	double w = 0.;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		w += D[i][j];
	    }
	}
	return w;
    }

    public double TIwiener () {
	return wiener (getDistanceMatrix ());
    }

    public double TIwiener3d () {
	return wiener (getDistanceMatrix3d ());
    }

    public double TImeanWiener () {
	double size = mgraph.getAtomCount();
	return 2.* TIwiener() / (size * (size - 1));
    }

    public double TIquasiWiener () {
	double[] eval = getLaplacianSolver().eigenvalues();
	double W = 0.;
	for (int i = 1; i < eval.length; ++i) {
	    W += 1/eval[i]; // we better not have disconnected components
	}
	W *= eval.length;
	return W;
    }

    public double TIevenWiener () {
	return evenWiener (1);
    }
    public double evenWiener (double p) {
	double[][] D = getDistanceMatrix ();
	double e = 0.;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		int d = (int)D[i][j];
		if (d % 2 == 0) { // even
		    e += Math.pow(D[i][j], p);
		}
	    }
	}
	return e;
    }

    public double TIoddWiener () {
	return oddWiener (1);
    }
    public double oddWiener (double q) {
	double[][] D = getDistanceMatrix ();
	double o = 0.;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		int d = (int)D[i][j];
		if (d % 2 != 0) { // odd
		    o += Math.pow(D[i][j], q);
		}
	    }
	}
	return o;
    }

    public double TIratioWiener () {
	return ratioWiener (1, 1);
    }
    public double ratioWiener (double p, double q) {
	return evenWiener (p) / oddWiener (q);
    }

    public double TImohar1 () {
	return 2.*Math.log((double)mgraph.getBondCount()
			   / mgraph.getAtomCount()) * TIquasiWiener();
    }

    public double TImohar2 () {
	double[] eval = getLaplacianSolver().eigenvalues();
	return 4./(eval.length * eval[1]);
    }

    // J. Chem. Inf: Comput. Sci. 1993, 33, 630-634
    public double TIvaa1 () {
	double[] eval = getAdjacencySolver().eigenvalues();
	double vaa = 0.;
	for (int i = 0; i < eval.length; ++i) {
	    if (eval[i] > 0.) {
		vaa += eval[i];
	    }
	}
	return vaa/eval.length;
    }

    public double TIvad1 () {
	double[] eval = getDistanceSolver().eigenvalues();
	double vad1 = 0.;
	for (int i = 0; i < eval.length; ++i) {
	    if (eval[i] < 0.) {
		vad1 += Math.abs(eval[i]);
	    }
	}
	return vad1/eval.length;
    }

    public double TIvel1 () {
	double[] evec = getLaplacianSolver().eigenvector(1);
	double v = 0.;
	for (int i = 0; i < evec.length; ++i) {
	    if (evec[i] < 0.) {
		v += Math.abs(evec[i]);
	    }
	}
	return v;
    }
    public double TIvel2 () {
	double[] evec = getLaplacianSolver().eigenvector(1);
	double v = 0.;
	for (int i = 0; i < evec.length; ++i) {
	    if (evec[i] > 0.) {
		v += Math.abs(evec[i]);
	    }
	}
	return v;
    }
    public double TIvel3 () {
	return TIvel1() / TIvel2();
    }

    public double TIgel1 () {
	double[] evec = getLaplacianSolver().eigenvector(1);
	Arrays.sort(evec);
	double gap = 0., d;
	for (int i = 1; i < evec.length; ++i) {
	    d = evec[i] - evec[i-1];
	    if (d > gap) {
		gap = d;
	    }
	}
	return gap;
    }

    public double TIdel1 () {
	double[] evec = getLaplacianSolver().eigenvector(1);
	Arrays.sort(evec);
	double s = 0.;
	for (int i = 1; i < evec.length; ++i) {
	    s += evec[i] - evec[i-1];
	}
	return s;
    }

    public double TIvea1 () {
	// first eigenvector
	double[] evec = getAdjacencySolver().eigenvector(0); 
	double vea1 = 0.;
	for (int i = 0; i < evec.length; ++i) {
	    vea1 += evec[i];
	}
	return vea1;
    }

    public double TIved1 () {
	double[] evec = getDistanceSolver().eigenvector(0);
	double ved1 = 0.;
	for (int i = 0; i < evec.length; ++i) {
	    ved1 += evec[i];
	}
	return ved1;
    }

    /*
    public double TIvra1 () {
	double[] eval = getAdjacencySolver().eigenvector(0);
	double vra1 = 0., d;
	for (MolBond b : mgraph.getBondArray()) {
	    int a1 = mgraph.indexOf(b.getAtom1());
	    int a2 = mgraph.indexOf(b.getAtom2());
	    d = Math.sqrt(eval[a1] * eval[a2]);
	    if (d > 0.) {
		vra1 += 1./d;
	    }
	}
	return vra1;
    }

    public double TIvrd1 () {
	double[] eval = getDistanceSolver().eigenvector(0);
	double vrd1 = 0., d = 0.;
	for (MolBond b : mgraph.getBondArray()) {
	    int a1 = mgraph.indexOf(b.getAtom1());
	    int a2 = mgraph.indexOf(b.getAtom2());
	    d = Math.sqrt(eval[a1] * eval[a2]);
	    if (d > 0.) {
		vrd1 += 1./d;
	    }
	}
	return vrd1;
    }
    */

    public double TIevD () {
	double[] eval = getDistanceSolver().eigenvalues();
	return norm (eval, 2);
    }

    public double TIbcut1 () {
	double[] eval = getBurdenSolver().eigenvalues();
	return eval[0];
    }

    /*
    public double TIspanningTree () {
	// FIXME: don't use this on normalized laplacian matrix!
	double[] eval = getLaplacianSolver().eigenvalues();
	double T = 1.;
	for (int i = 1; i < eval.length; ++i) {
	    T *= eval[i];
	}
	return T/eval.length;
    }
    */

    public double TIrandic () {
	double r = 0.;
	for (MolBond b : mgraph.getBondArray()) {
	    MolAtom a1 = b.getAtom1();
	    MolAtom a2 = b.getAtom2();
	    r += 1./Math.sqrt(a1.getBondCount() * a2.getBondCount());
	}
	return r;
    }

    // molecular topological index  (mti)
    /*
     * mti = sum_i [v(A + D)]_i
     */
    public double TIschultz () {
	double[][] A = getAdjacencyMatrix ();
	double[][] D = getDistanceMatrix ();

	int N = mgraph.getAtomCount();
	double mti = 0;

	double[] v = new double[N];
	for (int i = 0; i < N; ++i) {
	    v[i] = mgraph.getAtom(i).getBondCount();
	}

	for (int i = 0; i < N; ++i) {
	    for (int j = 0; j < N; ++j) {
		mti += v[j]*(A[j][i] + D[j][i]);
	    }
	}
	v = null;

	return mti;
    }

    public double harary (double[][] D) {
	double h = 0., x;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		x = D[i][j];
		h += 1./(x*x);
	    }
	}
	return h;
    }

    public double TIharary () {
	return harary (getDistanceMatrix ());
    }
    public double TIharary3d () {
	return harary (getDistanceMatrix3d ());
    }

    public double TIchi0 () {
	double x = 0.;
	for (MolAtom a : mgraph.getAtomArray()) {
	    x += 1./Math.sqrt(a.getBondCount());
	}
	return x;
    }

    public double TIchi0v () {
	double x = 0., v;
	for (MolAtom a : mgraph.getAtomArray()) {
	    v = a.getValence();
	    x += 1./Math.sqrt((v - a.getImplicitHcount())
			      / (a.getAtno() - v - 1));
	}
	return x;
    }

    public double TIchi1 () {
	double x = 0.;
	for (MolBond b : mgraph.getBondArray()) {
	    x += 1./Math.sqrt((double)(b.getAtom1().getBondCount()
				       + b.getAtom2().getBondCount()));
	}
	return x;
    }

    public double TIchi1v () {
	double x = 0., v1, v2;
	MolAtom a1, a2;
	for (MolBond b : mgraph.getBondArray()) {
	    a1 = b.getAtom1();
	    a2 = b.getAtom2();
	    v1 = (double)(a1.getValence() - a1.getImplicitHcount())
		/ (a1.getAtno() - a1.getValence() - 1);
	    v2 = (double)(a2.getValence() - a2.getImplicitHcount())
		/ (a2.getAtno() - a2.getValence() - 1);
	    x += 1./Math.sqrt(v1*v2);
	}
	return x;
    }

    // Bonchev-Trinasjstic's molecular branching:
    //   BT = n log2 n - sum_I n_I log2 n_I
    public double TIcomplexity1 () {
	double[][] D = getDistanceMatrix ();
	int[] nI = new int[D.length];
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		++nI[(int)D[i][j]];
	    }
	}
	double se = 0.;
	int n = 0;
	for (int i = 0; i < nI.length; ++i) {
	    if (nI[i] > 0) {
		n += nI[i];
		se += nI[i] * log2 (nI[i]);
	    }
	}
	return n*log2 (n) - se;
    }

    public double TIcomplexity2 () {
	double[][] D = getDistanceMatrix ();
	int[] nI = new int[D.length];
	int n = 0, d;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		d = (int)D[i][j];
		n += d;
		++nI[d];
	    }
	}
	double se = 0., pI;
	for (int i = 0; i < nI.length; ++i) {
	    if (nI[i] > 0) {
		pI = (double)nI[i]/n;
		se +=  pI* log2 (pI);
	    }
	}
	return -se;
    }

    protected int countDistance (int order) {
	double[][] D = getDistanceMatrix ();
	int d = 0;

	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		if ((int)D[i][j] == order) {
		    ++d;
		}
	    }
	}

	return d;
    }

    protected double kappaAlpha () {
	double alpha = 0.;
	for (MolAtom a : mgraph.getAtomArray()) {
	    int hyb = a.getHybridizationState();
	    alpha += KappaTable[a.getAtno()][hyb];
	}
	return alpha;
    }

    // kappa: J. Chem. Inf: Comput. Sci. 1993, 33, 630-634
    public double TIkappa1 () {
	double alpha = kappaAlpha ();
	double d = countDistance (1) + alpha;
	double k = 0.;
	if (d > 0) {
	    double n = mgraph.getAtomCount() + alpha;
	    k = (double)(n * (n-1) * (n-1))/(d*d);
	}
	return k;
    }

    public double TIkappa2 () {
	double alpha = kappaAlpha ();
	double d = countDistance (2) + alpha;
	double k = 0.;
	if (d > 0) {
	    double n = mgraph.getAtomCount() +alpha;
	    k = (double)((n-1) * (n-2) * (n-2))/(d*d);
	}
	return k;
    }

    public double TIkappa3 () {
	double alpha = kappaAlpha ();
	double d = countDistance (3) + alpha;
	double k = 0.;
	if (d > 0) {
	    double n = mgraph.getAtomCount() + alpha;
	    if (n % 2 == 0) { // even
		k = (double)((n-2) * (n-2) * (n-3))/(d*d);
	    }
	    else { // odd
		k = (double)((n-1) * (n-3) * (n-3))/(d*d);
	    }
	}
	return k;
    }

    public double TIkierflex () {
	return TIkappa1() * TIkappa2() / mgraph.getAtomCount();
    }

    protected synchronized MatrixSolver getLaplacianSolver () {
	if (_solverL == null) {
	    _solverL = new MatrixSolver (getLaplacianMatrix ());
	}
	return _solverL;
    }
    
    protected synchronized MatrixSolver getAdjacencySolver () {
	if (_solverA == null) {
	    _solverA = new MatrixSolver (getAdjacencyMatrix ());
	}
	return _solverA;
    }

    protected synchronized MatrixSolver getDistanceSolver () {
	if (_solverD == null) {
	    _solverD = new MatrixSolver (getDistanceMatrix ());
	}
	return _solverD;
    }

    protected synchronized MatrixSolver getBurdenSolver () {
	if (_solverB == null) {
	    _solverB = new MatrixSolver (getBurdenMatrix ());
	}
	return _solverB;
    }

    protected synchronized MatrixSolver getChargeSolver () {
	if (_solverC == null) {
	    _solverC = new MatrixSolver (getChargeMatrix ());
	}
	return _solverC;
    }

    // or algebraic connectivity index
    public double TIfiedler () {
	return getLaplacianSolver().getFiedlerValue();
    }

    // v is a vector of length equals to number of atom
    public double[] lovi (double[] a) {
	return getAdjacencySolver().solve(a);
	//return getDistanceSolver().solve(a);
	//return getLaplacianSolver().solve(a);
    }

    public double[] loviConstant (double x) {
	double[] c = new double[mgraph.getAtomCount()];
	for (int i = 0; i < c.length; ++i) {
	    c[i] = x;
	}
	return getDistanceSolver().solve(c);
    }

    public double[] loviCharge () {
	return loviCharge (7.4);
    }

    public double[] loviCharge (double pH) {
	double[] z = null;
	try {
	    ChargePlugin charge = new ChargePlugin ();
	    charge.setpH(pH);
	    charge.setMolecule(mgraph);
	    charge.run();

	    double[] v = new double[mgraph.getAtomCount()];
	    for (int i = 0; i < v.length; ++i) {
		v[i] = charge.getTotalCharge(i);
		if (Double.isNaN(v[i])) {
		    v[i] = 0.;
		}
	    }
	    z = lovi (v);
	    v = null;
	}
	catch (PluginException ex) {
	    ex.printStackTrace();
	}
	return z;
    }

    public double TIlovi1 () {
	return norm (loviCharge ());
    }

    public double[] loviVanderWaalsRadii () {
	double[] radii = new double[mgraph.getAtomCount()];
	for (int i = 0; i < radii.length; ++i) {
	    int z = mgraph.getAtom(i).getAtno();
	    switch (z) {
	    case 6: radii[i] = 1.70; break; // C
	    case 7: radii[i] = 1.55; break; // N
	    case 8: radii[i] = 1.52; break; // O
	    case 9: radii[i] = 1.47; break; // F
	    case 15: radii[i] = 1.80; break; // P
	    case 16: radii[i] = 1.80; break; // S
	    case 17: radii[i] = 1.75; break; // Cl
	    default: radii[i] = 1.40;
	    }
	}
	double[] z = lovi (radii);
	radii = null;

	return z;
    }
    
    public double TIlovi2 () {
	return norm (loviVanderWaalsRadii ());
    }

    public double[] loviCovalentRadii () {
	double[] r = new double[mgraph.getAtomCount()];
	for (int i = 0; i < r.length; ++i) {
	    int z = mgraph.getAtom(i).getAtno();
	    switch (z) {
	    case 6: r[i] = 0.76; break; // C.sp3
	    case 7: r[i] = 0.71; break; // N
	    case 8: r[i] = 0.66; break; // O
	    case 9: r[i] = 0.57; break; // F
	    case 15: r[i] = 1.07; break; // P
	    case 16: r[i] = 1.05; break; // S
	    case 17: r[i] = 1.02; break; // Cl
	    case 53: r[i] = 1.39; break; // I
	    default: r[i] = 1.5;
	    }
	}
	double[] z = lovi (r);
	r = null;
	return z;
    }

    public double TIlovi3 () {
	return norm (loviCovalentRadii ());
    }

    //J. Chem. Inf: Comput. Sci. 1993, 33, 630-634
    /* [...] evaluate the charge transfers between pairs of 
       atoms, and therefore the lobal charge transfers in the 
       molecule.
    */
    public double chargeIndex (int order) {
	MatrixSolver M = getChargeSolver ();
	double[][] D = getDistanceMatrix ();
	double index = 0.;
	for (int i = 0; i < D.length; ++i) {
	    for (int j = i+1; j < D.length; ++j) {
		if (order == (int)D[i][j]) {
		    index += Math.abs(M.get(i, j) - M.get(j, i)) * D[i][j];
		}
	    }
	}
	return index / (D.length-1);
    }

    public double TIchargeDist1 () {
	return chargeIndex (1);
    }
    public double TIchargeDist2 () {
	return chargeIndex (2);
    }
    public double TIchargeDist3 () {
	return chargeIndex (3);
    }

    // Kier-Hall intrinsic value
    protected static double intrinsicValue (MolAtom a) {
	double I = 0.;

	// single, double, triple, and aromatic bonds
	int sb = 0, db = 0, tb = 0, ab = 0;
	for (int i = 0; i < a.getBondCount(); ++i) {
	    MolBond b = a.getBond(i);
	    switch (b.getType()) {
	    case 1: ++sb; break;
	    case 2: ++db; break;
	    case 3: ++tb; break;
	    case MolBond.AROMATIC: ++ab; break;
	    }
	}
	int h = a.getImplicitHcount();

	switch (a.getAtno()) {
	case 6: // C
	    {
		if (sb == 1 && h == 3) {
		    I = 2.; // -CH3
		}
		else if (sb == 2 && h == 2) {
		    I = 1.5; // -CH2-
		}
		else if (sb == 3 && h == 1) {
		    I = 1.33; // >CH-
		}
		else if (sb == 4) {
		    I = 1.25; // >C<
		}
		else if (ab > 0 && sb > 0) {
		    I = 1.667; // 
		}
		else if (ab > 0 && sb == 1) {
		    I = 2.;
		}
		else if (db == 1 && h == 2) {
		    I = 3.; // =CH2
		}
		else if (db == 1 && sb == 1 && h == 1) {
		    I = 2.; // =CH-
		}
		else if (db == 2) {
		    I = 2.5; // =C=
		}
		else if (db == 1 && sb == 2) {
		    I = 1.667; // =C<
		}
		else if (tb == 1 && h == 1) {
		    I = 4.; // #CH
		}
		else if (tb == 1 && sb == 1) {
		    I = 2.5; // #C-
		}
	    }
	    break;

	case 7:
	    {
		if (sb == 1 && h == 2) {
		    I = 4.; // -NH2
		}
		else if (sb == 2 && h == 1) {
		    I = 2.5; // -NH-
		}
		else if (sb == 3 && db == 0) {
		    I = 2.0; // >N-
		}
		else if (ab > 0 && h == 0) {
		    I = 2.5; // nH
		}
		else if (ab > 0 && h == 1) {
		    I = 3.0; // n
		}
		else if (db == 1 && h == 1) {
		    I = 5.0; // =NH
		}
		else if (db == 1 && sb == 1) {
		    I = 3.0; // =N-
		}
		else if (ab > 0 && db == 1) {
		    I = 2.0; // n=
		}
		else if (db == 2 && sb == 1) {
		    I = 2.0; // -N=(=)
		}
		else if (tb == 3) {
		    I = 6.0; // #N
		}
	    }
	    break;

	case 8:
	    {
		if (sb == 1 && h == 1) {
		    I = 6.0; // -OH
		}
		else if (sb == 2) {
		    I = 3.5; // -O-
		}
		else if (ab > 0) {
		    I = 3.5; // o
		}
		else if (db == 1) {
		    I = 7.0; // =O
		}
	    }
	    break;

	case 9:
	    {
		if (sb == 1) {
		    I = 8.0; // -F
		}
	    }
	    break;

	case 15:
	    {
		if (sb == 3 && db == 0) {
		    I = 1.074;
		}
		else if (sb == 3 && db == 1) {
		    I = 0.806;
		}
	    }
	    break;

	case 16:
	    {
		if (sb == 1 && h == 1) {
		    I = 3.222; // -SH
		}
		else if (sb == 2 && db == 0) {
		    I = 1.833; // -S-
		}
		else if (ab > 0) {
		    I = 1.833; // s
		}
		else if (db == 1 && sb == 0) {
		    I = 3.667; // =S
		}
		else if (sb == 2 && db == 1) {
		    I = 1.222; // =S<
		}
		else if (sb == 2 && db == 2) {
		    I = 0.917;
		}
	    }
	    break;

	case 17:
	    {
		if (sb == 1) {
		    I = 4.111; // -Cl
		}
	    }
	    break;

	case 35:
	    {
		if (sb == 1) {
		    I = 2.75; // -Br
		}
	    }
	    break;

	case 53:
	    {
		if (sb == 1) {
		    I = 2.12; // -I
		}
	    }
	    break;
	}
	return I;
    }

    public double[] eState (double order) {
	MolAtom[] atoms = mgraph.getAtomArray();
	/*
	 * NOTE: because only the relative intrinsic values are being
	 * considered, it doesn't matter here whether the topological 
	 * distance or 3d distance matrix is used here!
	 */
	double[][] D = getDistanceMatrix ();
	double[] S = new double[atoms.length];
	double dI, I;
	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom a = atoms[i];
	    // compute intrinsic value for this atom
	    I = intrinsicValue (a);
	    dI = 0.;
	    for (int j = 0; j < a.getBondCount(); ++j) {
		MolAtom b = a.getBond(j).getOtherAtom(a);
		dI += (I - intrinsicValue (b)) 
		    / Math.pow(D[i][mgraph.indexOf(b)]+1, order);
	    }
	    S[i] = I + dI;
	}
	return S;
    }

    public double TIeState () {
	return sum (eState (2.));
    }

    public static double sum (double[] v) {
	double s = 0.;
	for (int i = 0; i < v.length; ++i) {
	    s += v[i];
	}
	return s;
    }

    public static double norm (double[] v) {
	return norm (v, 2); // euclidean norm
    }

    public static double norm (double[] v, int order) {
	double lv = -1;
	if (order < 0) {
	    throw new IllegalArgumentException 
		("Invalid norm specified: " + order);
	}

	switch (order) {
	case 0: // inf-norm
	    {
		lv = Double.MIN_VALUE;
		double x;
		for (int i = 0; i < v.length; ++i) {
		    x = Math.abs(v[i]);
		    if (x > lv) {
			lv = x;
		    }
		}
	    }
	    break;

	case 1:
	    {
		lv = 0.;
		for (int i = 0; i < v.length; ++i) {
		    lv += Math.abs(v[i]);
		}
	    }
	    break;

	case 2:
	    {
		lv = 0.;
		for (int i = 0; i < v.length; ++i) {
		    lv += v[i]*v[i];
		}
		lv = Math.sqrt(lv);
	    }
	    break;

	default:
	    {
		lv = 0.;
		double x;
		for (int i = 0; i < v.length; ++i) {
		    x = Math.abs(v[i]);
		    for (int j = 1; j < order; ++j) {
			x *= x;
		    }
		    lv += x;
		}
		lv = Math.pow(lv, 1./order);
	    }
	}

	return lv;
    }

    public static double dotKernel (double[] u, double[] v) {
	return dotKernel (u, v, null);
    }

    public static double dotKernel (double[] u, double[] v, int[] uvmap) {
	if (uvmap == null && u.length != v.length) {
	    throw new IllegalArgumentException
		("Input dimensions for input vectors don't match");
	}
	double k = 0.;
	if (uvmap == null) {
	    for (int i = 0; i < u.length; ++i) {
		k += u[i]*v[i];
	    }
	}
	else {
	    for (int i = 0; i < uvmap.length; ++i) {
		int j = uvmap[i];
		if (j >= 0) {
		    k += u[i] * v[j];
		}
	    }
	}
	return k;
    }

    public static double rbfKernel (double[] u, double[] v, int[] uvmap) {
	double gamma = 0.;
	if (uvmap == null) {
	    gamma = 1./u.length;
	}
	else {
	    for (int i = 0; i < uvmap.length; ++i) {
		if (uvmap[i] > 0) {
		    ++gamma;
		}
	    }
	    if (gamma > 0.) {
		gamma = 1./gamma;
	    }
	}
	return rbfKernel (gamma, u, v, uvmap);
    }

    public static double rbfKernel (double gamma, double[] u, 
				    double[] v, int[] uvmap) {
	if (uvmap == null && u.length != v.length) {
	    throw new IllegalArgumentException
		("Input dimensions for input vectors don't match");
	}
	double k = 0., d;
	if (uvmap == null) {
	    for (int i = 0; i < u.length; ++i) {
		d = u[i] - v[i];
		k += d*d;
	    }
	}
	else {
	    for (int i = 0; i < uvmap.length; ++i) {
		int j = uvmap[i];
		if (j >= 0) {
		    d = u[i] - v[j];
		    k += d*d;
		}
	    }
	}
	return Math.exp(-gamma*k);
    }

    public static void main (String[] argv) throws Exception {
	MolImporter mi;
	if (argv.length == 0) {
	    System.err.println("** reading from STDIN...");
	    mi = new MolImporter (System.in);
	}
	else {
	    mi = new MolImporter (argv[0]);
	}

	String[] desc = TopologicalIndices.getNames();
	/*
	System.out.print("Compound");
	for (int i = 0; i < desc.length; ++i) {
	    System.out.print("," + desc[i]);
	}
	System.out.println();
	*/
	java.util.Vector<Molecule> mols = new java.util.Vector<Molecule>();
	java.util.Vector<TopologicalIndices> indices = 
	    new java.util.Vector<TopologicalIndices>();
	gov.nih.ncgc.util.MolStandardizer molstd = 
	    new gov.nih.ncgc.util.MolStandardizer();

	double[] mins = new double[desc.length];
	java.util.Arrays.fill(mins, Double.MAX_VALUE);
	double[] maxs = new double[desc.length];
	java.util.Arrays.fill(maxs, Double.MIN_VALUE);

	java.util.Map<String, double[]> data0 = 
	    new java.util.HashMap<String, double[]>();
	java.util.Map<String, double[]> data1 = 
	    new java.util.HashMap<String, double[]>();

	for (Molecule mol; (mol = mi.read()) != null; ) {
	    if (mol.getAtomCount() < 3) {
		continue;
	    }

	    boolean active = false;
	    String prop = mol.getProperty("herg-trans-p2/set1-CurveClass");
	    if (prop != null) {
		double c = Double.parseDouble(prop);
		active = c == -1.1 || c == -2.1;
	    }
	    else {
		prop = mol.getProperty("herg-trans-p2/set2-CurveClass");
		if (prop != null) {
		    double c = Double.parseDouble(prop);
		    active = c == -1.1 || c == -2.1;
		}
		else {
		    prop = mol.getProperty("herg-trans-p1/set1-CurveClass");
		    if (prop != null) {
			double c = Double.parseDouble(prop);
			active = c == -1.1 || c == -2.1;
		    }
		}
	    }


	    mol.calcHybridization();
	    System.out.println(mol.getName());
	    /*
	    for (int i = 0; i < mol.getAtomCount(); ++i) {
		MolAtom a = mol.getAtom(i);
		System.out.println((i+1) + ": hyb=" 
				   + a.getHybridizationState() + " v=" 
				   + a.getValence() + " z="
				   + a.getAtno() + " h="
				   + a.getImplicitHcount() + " lp="
				   + mol.getLonePairCount(i));
	    }
	    */

	    molstd.standardize(mol);
	    mol = molstd.getLargestFragment();
	    if (mol.getAtomCount() < 3) {
		continue;
	    }

	    double[] row = new double[desc.length];
	    TopologicalIndices ti = new TopologicalIndices (mol);
	    for (int i = 0; i < desc.length; ++i) {
		double v = (Double)ti.getClass().getMethod(desc[i]).invoke(ti);
		row[i] = v;
		if (Double.isNaN(v) || Double.isInfinite(v)) {
		    System.err.println("** warning: " + desc[i] + " for " 
				       + mol.getName() + " is " + v);
		}
		else {
		    if (v > maxs[i]) maxs[i] = v;
		    if (v < mins[i]) mins[i] = v;
		}

		//System.out.printf(",%1$.4f", v);
		//System.out.printf("\n%1$2d%2$20s %3$.4f", i+1,desc[i], v);
	    }

	    if (active) {
		data1.put(mol.getName(), row);
	    }
	    else {
		data0.put(mol.getName(), row);
	    }

	    /*
	    double[] lovi = ti.loviCharge();
	    for (int i = 0; i < lovi.length; ++i) {
		System.out.print("," + lovi[i]);
	    }
	    */
	    //System.out.println();

	    /*
	    double[] charge = ti.loviCharge();
	    double[] vwaals = ti.loviVanderWaalsRadii();
	    double[] covalent = ti.loviCovalentRadii();
	    System.out.println("        Charge     vWaals     Covalent");
	    for (int i = 0; i < charge.length; ++i) {
		System.out.printf("%1$2d: %2$10.4f %3$10.4f %4$10.4f\n",
				  i+1,charge[i], vwaals[i], covalent[i]);
	    }
	    */
	    mols.add(mol);
	    indices.add(ti);
	}
	mi.close();

	System.out.println(data1.size() + " active sample(s)!");
	System.out.println(data0.size() + " inactive sample(s)!");
	System.out.println("now generatin prob density functions...");
	gov.nih.ncgc.math.stat.DensityEstimator pdf1[] = new
	    gov.nih.ncgc.math.stat.DensityEstimator[desc.length];
	gov.nih.ncgc.math.stat.DensityEstimator pdf0[] = new
	    gov.nih.ncgc.math.stat.DensityEstimator[desc.length];
	gov.nih.ncgc.math.stat.DensityEstimator pdf;
	double lower = -5, upper = 5;

	for (double[] v : data1.values()) {
	    for (int i = 0; i < v.length; ++i) {
		pdf = pdf1[i];
		double d = maxs[i] - mins[i];
		if (pdf == null) {
		    pdf1[i] = pdf = new gov.nih.ncgc.math
			.stat.ContinuousDensityEstimator(0.2, 100);
		    pdf.setRangeUniform(lower, upper);
		}
		pdf.increment(lower + (upper-lower)*(v[i]-mins[i])/d);
	    }
	}

	for (double[] v : data0.values()) {
	    for (int i = 0; i < v.length; ++i) {
		pdf = pdf0[i];
		double d = maxs[i] - mins[i];
		if (pdf == null) {
		    pdf0[i] = pdf = new gov.nih.ncgc.math
			.stat.ContinuousDensityEstimator(0.2, 100);
		    pdf.setRangeUniform(lower, upper);
		}
		pdf.increment(lower + (upper-lower)*(v[i]-mins[i])/d);
	    }
	}

	double[] x = new double[500];
	java.util.Random rand = new java.util.Random();
	for (int i = 0; i < desc.length; ++i) {
	    pdf1[i].estimate();
	    pdf0[i].estimate();
	    for (int j = 0; j < x.length; ++j) {
		x[j] = lower + (upper-lower)*rand.nextDouble();
	    }
	    java.util.Arrays.sort(x);
	    java.io.PrintStream ps1 = new java.io.PrintStream
		(new java.io.FileOutputStream (desc[i] + "_1.txt"));
	    java.io.PrintStream ps0 = new java.io.PrintStream
		(new java.io.FileOutputStream (desc[i] + "_0.txt"));

	    for (int j = 0; j < x.length; ++j) {
		ps1.println(x[j] + " " + pdf1[i].probability(x[j]));
		ps0.println(x[j] + " " + pdf0[i].probability(x[j]));
	    }
	    ps1.close();
	    ps0.close();
	}
    }
}
