package com.akasiyanik.problems;

/**
 * @author akasiyanik
 */
public class PoissonProblem {

    private int n;

    private double[][] J;

    public PoissonProblem(int n) {
        this.n = n;
        this.J = calculateJ();
    }

    private double[][] calculateJ() {
        double[][] J = new double[n][n];
        double tmp = (n + 1) * (n + 1);
        for (int i = 0; i < n; i++) {
            J[i][i] = - 2 * tmp;
            if (i < n - 1) {
                J[i][i + 1] = tmp;
            }
            if (i > 0) {
                J[i][i - 1] = tmp;
            }
        }
        return J;
    }

    public double[] getY0() {
        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0 + i;
        }
        return y0;
    }

    public double[][] getJ() {
        return J;
    }

    public int getN() {
        return n;
    }
}
