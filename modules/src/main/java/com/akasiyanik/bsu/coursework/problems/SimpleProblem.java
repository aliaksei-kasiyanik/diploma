package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

/**
 * @author: akasiyanik
 */
public class SimpleProblem extends SteadyingEquation {
    public SimpleProblem(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public SimpleProblem(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    @Override
    public double[] f(double t, double[] arg) {

        double f1 = -Math.pow(arg[0], 2) - arg[1];
        double f2 = 2.0 * arg[0] - Math.pow(arg[1], 3);

        double[] res = new double[]{f1, f2};
        return res;
    }

    @Override
    public double[][] JacobiMatrix(double[] arg) {
        double df1_darg0 = -2.0 * arg[0];
        double df1_darg1 = -1.0;
        double df2_darg0 = 2.0;
        double df2_darg1 = -3.0 * Math.pow(arg[1], 2);

        double[][] matr = {
                {df1_darg0, df1_darg1},
                {df2_darg0, df2_darg1}
        };
        return matr;
    }

    @Override
    public double getJacobiMatrMaxEigenvalue(double t, double[] arg) {// < 0
        //temporary!!!!!!!!!!!
        return -2.5;
    }

    @Override
    public double getAMatrMaxEigenvalue() {  // > 0
        //temporary!!!!!!!!!!!
        return 0.3333333;
    }
}
