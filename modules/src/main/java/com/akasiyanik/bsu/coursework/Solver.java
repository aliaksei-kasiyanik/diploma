package com.akasiyanik.bsu.coursework;

/**
 * @author: akasiyanik
 */
public interface Solver {

    public double[] solve(double t0, double[] y0, double tau, double eps);

    public int getIterationCount();

    public int getIterationCountWithClarifying();

}
