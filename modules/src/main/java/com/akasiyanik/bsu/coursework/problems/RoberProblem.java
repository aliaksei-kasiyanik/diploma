package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

/**
 * @author akasiyanik
 */
public class RoberProblem extends SteadyingEquation {

    public RoberProblem(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public RoberProblem(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    @Override
    public double[][] getG() {
        return new double[0][];
    }

    @Override
    public double[] f(double t, double[] arg) {
        double f0 = -0.04 * arg[0] + Math.pow(10, 4) * arg[1] * arg[2];
        double f1 = 0.04 * arg[0] - Math.pow(10, 4) * arg[1] * arg[2] - 3 * Math.pow(10, 7) * Math.pow(arg[1], 2);
        double f2 = 3 * Math.pow(10, 7) * Math.pow(arg[1], 2);


        double[] res = {f0, f1, f2};
        return res;
    }

    @Override
    public double[][] JacobiMatrix(double[] arg) {
        return new double[0][];
    }

    @Override
    public double getJacobiMatrMaxEigenvalue(double t, double[] arg) {
        double[] m = new double[3];
        m[0] = 0.04 + Math.pow(10, 4)* arg[2] + Math.pow(10, 4)* arg[1];
        m[1] = 0.04 + Math.pow(10, 4) * arg[2] + 3 * 2 * Math.pow(10, 7) * arg[1] + Math.pow(10, 4) * arg[1];
        m[2] = 2 * 3 * Math.pow(10, 7) * arg[1];

        double max_j = 0.0;
        for (int i = 0; i < 3; i++)
            if (m[i] > max_j)
                max_j = m[i];

        return max_j;
    }

    @Override
    public double getAMatrMaxEigenvalue() {
        return 0.3333333;//Radau - 3
//        return 0.2748888295956;//Radau - 5
    }
}
