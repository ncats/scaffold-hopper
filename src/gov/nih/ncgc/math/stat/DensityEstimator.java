// $Id: DensityEstimator.java 3490 2009-10-29 15:51:13Z nguyenda $

package gov.nih.ncgc.math.stat;

import gov.nih.ncgc.math.Histogram;

/*
 * probability density estimator.
 */

public abstract class DensityEstimator extends Histogram {
    static final double ONE_LOG2 = 1.44269504088896340737;

    protected DensityEstimator (int nbins) {
	super (nbins);
    }
    protected DensityEstimator (int nbins, double min, double max) {
	super (nbins, min, max);
    }

    public abstract double probability (double x);
    public abstract void estimate ();

    // compute (symmetric) kullback-leibler distance
    public double klDistance (DensityEstimator pdf) {
	if (pdf.size() != bin.length) {
	    throw new IllegalArgumentException 
		("Input density is not compatible");
	}
	double d12 = 0., d21 = 0.;
	for (int i = 0; i < bin.length; ++i) {
	    double x1 = bin[i];
	    double x2 = pdf.bin[i];
	    d12 += x1 * Math.log(x1/x2);
	    d21 += x2 * Math.log(x2/x1);
	}
	return d12 + d21;
    }

    // compute differential shannon entropy
    public double entropyDiff (DensityEstimator pdf) {
	if (pdf.size() != bin.length) {
	    throw new IllegalArgumentException 
		("Input density is not compatible");
	}
	
	double s1 = 0, s2 = 0, s12 = 0.;
	for (int i = 0; i < bin.length; ++i) {
	    double x1 = bin[i], x2 = pdf.bin[i];
	    s1 += x1 * ONE_LOG2 * Math.log(x1);
	    s2 += x2 * ONE_LOG2 * Math.log(x2);
	    s12 += (x1+x2)*ONE_LOG2 * Math.log(x1+x2);
	}
	return (s1+s2)/2. - s12;
    }
}
