package gov.nih.ncgc.rgroup;

import chemaxon.struc.Molecule;

public class ScaffoldData {
    double mean, std;

    public ScaffoldData (String prop, RGroupTable tab) {
        mean = 0.;
        int c = 0;
        Molecule[] mb = tab.getMembers();
        double[] vx = new double[mb.length];
        for (Molecule m : mb) {
            try {
                String val = m.getProperty(prop);
                if (val != null) {
                    double x = Double.parseDouble(val);
                    mean += x;
                    vx[c++] = x;
                }
            }
            catch (NumberFormatException ex) {
            }
        }

        if (c > 0) {
            mean /= c;
        }

        if (c > 1) {
            std = 0.;
            for (int i = 0; i < c ; ++i) {
                double x = (vx[i] - mean);
                std += x*x;
            }
            std /= (c-1);
            std = Math.sqrt(std);
        }
    }

    public double getMean () { return mean; }
    public double getStd () { return std; }
}
