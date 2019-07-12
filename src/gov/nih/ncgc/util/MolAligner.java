// $Id: MolAligner.java 3451 2009-10-24 00:55:23Z nguyenda $
package gov.nih.ncgc.util;

import java.io.*;

import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.*;

import cern.colt.matrix.*;
import cern.colt.matrix.linalg.*;
import cern.colt.function.DoubleFunction;

public class MolAligner {
    private MolAligner () {
    }

    /**
     * Perform aligment of the target molecule given the query molecule
     * and the atom mappings.  The map array is indexed by the atom index
     * of ref, i.e., map[i] = j implies that the ith atom of query corresponds
     * to the jth atom of target where j >= 0, otherwise there is no mapping.
     */
    public static void align (Molecule query, Molecule target, int[] map) {
	if (query.getAtomCount() != map.length) {
	    throw new IllegalArgumentException 
		("Input mapping does match reference molecule size");
	}

	if (query.getDim() < 2) {
	    /*
	    throw new IllegalArgumentException 
		("Query molecule doesn't have coordinates");
	    */
	    query.clean(2, null); // default is to generate 2d coordinates
	}

	int size = 0;
	double qx = 0., qy = 0., qz = 0.;
	double tx = 0., ty = 0., tz = 0.;

	MolAtom[] qatoms = query.getAtomArray();
	MolAtom[] tatoms = target.getAtomArray();
	for (int i = 0; i < map.length; ++i) {
	    if (map[i] < 0) {
	    }
	    else {
		MolAtom q = qatoms[i];
		MolAtom t = tatoms[map[i]];
		qx += q.getX();
		qy += q.getY();
		qz += q.getZ();
		tx += t.getX();
		ty += t.getY();
		tz += t.getZ();
		++size;
	    }
	}

	qx /= size;
	qy /= size;
	qz /= size;
	tx /= size;
	ty /= size;
	tz /= size;

	// now center the vectors
	DoubleMatrix2D Y = DoubleFactory2D.dense.make(3, size);
	DoubleMatrix2D X = DoubleFactory2D.dense.make(3, size);
	for (int i = 0, j = 0; i < map.length; ++i) {
	    if (map[i] < 0) {
	    }
	    else {
		MolAtom q = qatoms[i];
		MolAtom t = tatoms[map[i]];
		
		X.setQuick(0, j, q.getX() - qx);
		X.setQuick(1, j, q.getY() - qy);
		X.setQuick(2, j, q.getZ() - qz);
		Y.setQuick(0, j, t.getX() - tx);
		Y.setQuick(1, j, t.getY() - ty);
		Y.setQuick(2, j, t.getZ() - tz);
		++j;
	    }
	}

	//System.out.println("Y = " + Y);
	//System.out.println("X = " + X);

	Algebra alg = new Algebra ();

	// now compute A = YX'
	DoubleMatrix2D A = alg.mult(Y, alg.transpose(X));	
	//System.out.println("A = " + A);

	SingularValueDecomposition svd = new SingularValueDecomposition (A);
	DoubleMatrix2D U = svd.getU();
	DoubleMatrix2D V = svd.getV();
	DoubleMatrix2D R = alg.mult(U, alg.transpose(V));
	//System.out.println("R = " + R);
	R = alg.transpose(R);

	/*
	// if A is non-singular, then solve for the rotation matrix R
	//   as R = A(A'A)^{-1/2}
	DoubleMatrix2D invA = alg.inverse(alg.mult(alg.transpose(A), A));
	//System.out.println("inv(A) = " + invA);
	invA.assign(new DoubleFunction () {
		public double apply (double x) {
		    return Math.sqrt(x);
		}
	    });

	// the rotation matrix R; if det(R) = -1, then we have rotoinversion
	//  (rotation and reflection).  det(R) = 1, only rotation.
	DoubleMatrix2D R = alg.mult(A, invA);
	System.out.println("R = " + R);
	double det = alg.det(R);
	if (det < 0.) {
	    System.out.println("** R is rotoinversion");
	}
	else {
	    R = alg.transpose(R);
	}
	System.out.println("det R = " + det);
	*/

	// compute translation... 
	/*
	double zx = tx - (R.getQuick(0, 0)*qx + R.getQuick(0, 1)*qy);
	double zy = ty - (R.getQuick(1, 0)*qx + R.getQuick(1, 1)*qy);
	System.out.println("translation: zx = " + zx+" zy = " + zy);
	*/
	
	
	// rotate the target atoms..
	for (MolAtom a : tatoms) {
	    double x = a.getX();
	    double y = a.getY();
	    double z = a.getZ();
	    double rx = x*R.getQuick(0, 0) 
		+ y*R.getQuick(0, 1) + z*R.getQuick(0, 2);
	    double ry = x*R.getQuick(1, 0) 
		+ y*R.getQuick(1, 1) + z*R.getQuick(1, 2);
	    double rz = x*R.getQuick(2, 0) 
		+ y*R.getQuick(2, 1) + z*R.getQuick(2, 2);
	    a.setXYZ(rx, ry, rz);
	}
		//if rotoinversion, invert dashes and wedges
		double det = R.getQuick(0, 0)*R.getQuick(1, 1)-R.getQuick(0, 1)*R.getQuick(1, 0);
		if(det<0){
			for(MolBond mb:target.getBondArray()){
				int ster=mb.getStereo1(mb.getAtom1());
				if(ster==MolBond.UP){
					mb.stepWedge();
				}else if(ster==MolBond.DOWN){
					mb.stepWedge();
					mb.stepWedge();
				}
			}
		}
    }

    public static void main (String[] argv) throws Exception {
	MolHandler mh = new MolHandler ();

	mh.setMolecule(
"NCGC00028940-01\n"+
"  Marvin  10140910092D          \n"+
"\n"+
" 19 21  0  0  0  0            999 V2000\n"+
"    2.3155    1.2734    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.6011    0.8609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.8866    1.2734    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.1721    0.8609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -0.6125    1.1159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -0.9480    1.8696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -1.7685    1.9558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -2.2534    1.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -1.9179    0.5347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -1.0974    0.4484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -0.6125   -0.2190    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.1721    0.0359    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.8866   -0.3766    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.6011    0.0359    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.8866   -1.2016    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.1721   -1.6141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.1721   -2.4391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    0.8866   -2.8516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -0.5423   -2.8516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  1  0  0  0  0\n"+
"  2  3  1  0  0  0  0\n"+
"  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n"+
"  5  6  2  0  0  0  0\n"+
"  6  7  1  0  0  0  0\n"+
"  7  8  2  0  0  0  0\n"+
"  8  9  1  0  0  0  0\n"+
"  9 10  2  0  0  0  0\n"+
"  5 10  1  0  0  0  0\n"+
" 10 11  1  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n"+
"  4 12  1  0  0  0  0\n"+
" 12 13  2  0  0  0  0\n"+
" 13 14  1  0  0  0  0\n"+
"  2 14  2  0  0  0  0\n"+
" 13 15  1  0  0  0  0\n"+
" 15 16  1  0  0  0  0\n"+
" 16 17  1  0  0  0  0\n"+
" 17 18  2  0  0  0  0\n"+
" 17 19  1  0  0  0  0\n"+
"M  END\n");

	Molecule target = mh.getMolecule();
	target.aromatize();

	mh.setMolecule(
"\n"+
"  Marvin  10220916572D          \n"+
"\n"+
" 13 15  0  0  0  0            999 V2000\n"+
"    1.9749   -1.8326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.7200   -1.0479    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.2720   -0.4349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.0790   -0.6064    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.7464   -0.1215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4139   -0.6064    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2208   -0.4349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.7729   -1.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.5179   -1.8326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7110   -2.0041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.1589   -1.3910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3339   -1.3910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.7819   -2.0041    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  2  0  0  0  0\n"+
"  2  3  1  0  0  0  0\n"+
"  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n"+
"  5  6  1  0  0  0  0\n"+
"  6  7  2  0  0  0  0\n"+
"  7  8  1  0  0  0  0\n"+
"  8  9  2  0  0  0  0\n"+
"  9 10  1  0  0  0  0\n"+
" 10 11  2  0  0  0  0\n"+
"  6 11  1  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n"+
"  4 12  1  0  0  0  0\n"+
" 12 13  2  0  0  0  0\n"+
"  1 13  1  0  0  0  0\n"+
"M  END\n");
	Molecule query = mh.getMolecule();
	query.setName("Query1");
	query.aromatize();

	mh.setMolecule(
"\n"+
"  Marvin  10220921442D          \n"+
"\n"+
" 13 15  0  0  0  0            999 V2000\n"+
"    5.5169   -1.8352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.7730   -1.0509    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2219   -0.4370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4146   -0.6074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.7479   -0.1215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.0798   -0.6054    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.2730   -0.4327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.7201   -1.0450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.9739   -1.8300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.7806   -2.0027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3336   -1.3904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.1586   -1.3916    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7097   -2.0055    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  2  0  0  0  0\n"+
"  2  3  1  0  0  0  0\n"+
"  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n"+
"  5  6  1  0  0  0  0\n"+
"  6  7  2  0  0  0  0\n"+
"  7  8  1  0  0  0  0\n"+
"  8  9  2  0  0  0  0\n"+
"  9 10  1  0  0  0  0\n"+
" 10 11  2  0  0  0  0\n"+
"  6 11  1  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n"+
"  4 12  1  0  0  0  0\n"+
" 12 13  2  0  0  0  0\n"+
"  1 13  1  0  0  0  0\n"+
"M  END\n");
	Molecule query2 = mh.getMolecule();
	query2.setName("Query2");

	mh.setMolecule(
"\n"+
"  Marvin  10220923552D          \n"+
"\n"+
" 13 15  0  0  0  0            999 V2000\n"+
"    5.5169   -0.4349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.7730   -1.2191    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2219   -1.8330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4146   -1.6627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.7479   -2.1486    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.0798   -1.6646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.2730   -1.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.7201   -1.2250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.9739   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.7806   -0.2673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3336   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.1586   -0.8784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7097   -0.2645    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  2  0  0  0  0\n"+
"  2  3  1  0  0  0  0\n"+
"  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n"+
"  5  6  1  0  0  0  0\n"+
"  6  7  2  0  0  0  0\n"+
"  7  8  1  0  0  0  0\n"+
"  8  9  2  0  0  0  0\n"+
"  9 10  1  0  0  0  0\n"+
" 10 11  2  0  0  0  0\n"+
"  6 11  1  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n"+
"  4 12  1  0  0  0  0\n"+
" 12 13  2  0  0  0  0\n"+
"  1 13  1  0  0  0  0\n"+
"M  END\n"
);
	Molecule query3 = mh.getMolecule();
	query3.setName("Query3");

	mh.setMolecule(
"\n"+
"  Marvin  10230900082D          \n"+
"\n"+
" 13 15  0  0  0  0            999 V2000\n"+
"    1.9759   -0.4349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    1.7199   -1.2191    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.2710   -1.8330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.0782   -1.6627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.7450   -2.1486    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4131   -1.6646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2198   -1.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.7727   -1.2250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.5190   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7122   -0.2673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.1593   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3343   -0.8784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.7832   -0.2645    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  2  0  0  0  0\n"+
"  2  3  1  0  0  0  0\n"+
"  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n"+
"  5  6  1  0  0  0  0\n"+
"  6  7  2  0  0  0  0\n"+
"  7  8  1  0  0  0  0\n"+
"  8  9  2  0  0  0  0\n"+
"  9 10  1  0  0  0  0\n"+
" 10 11  2  0  0  0  0\n"+
"  6 11  1  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n"+
"  4 12  1  0  0  0  0\n"+
" 12 13  2  0  0  0  0\n"+
"  1 13  1  0  0  0  0\n"+
"M  END\n");
	Molecule query4 = mh.getMolecule();
	query4.setName("Query4");

	MolSearch ms = new MolSearch ();
	ms.setQuery(query);
	ms.setTarget(target);

	int[] hit = ms.findFirst();

	PrintStream ps = new PrintStream 
	    (new FileOutputStream ("align-out.sdf"));
	ps.print(target.toFormat("sdf"));

	Molecule m0 = target.cloneMolecule();
	ps.print(query.toFormat("sdf"));
	MolAligner.align(query, m0, hit);
	ps.print(m0.toFormat("sdf"));

	Molecule m1 = target.cloneMolecule();
	ps.print(query2.toFormat("sdf"));
	MolAligner.align(query2, m1, hit);
	ps.print(m1.toFormat("sdf"));
	MolAligner.align(query2, m0, hit);
	ps.print(m0.toFormat("sdf"));

	Molecule m2 = target.cloneMolecule();
	ps.print(query3.toFormat("sdf"));
	MolAligner.align(query3, m2, hit);
	ps.print(m2.toFormat("sdf"));
	MolAligner.align(query3, m1, hit);
	ps.print(m1.toFormat("sdf"));
	MolAligner.align(query3, m0, hit);
	ps.print(m0.toFormat("sdf"));

	Molecule m3 = target.cloneMolecule();
	ps.print(query4.toFormat("sdf"));
	MolAligner.align(query4, m3, hit);
	ps.print(m3.toFormat("sdf"));
	MolAligner.align(query4, m2, hit);
	ps.print(m2.toFormat("sdf"));
	MolAligner.align(query4, m1, hit);
	ps.print(m1.toFormat("sdf"));
	MolAligner.align(query4, m0, hit);
	ps.print(m0.toFormat("sdf"));

	ps.close();
    }
}
