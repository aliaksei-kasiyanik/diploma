package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

/**
 * author: Aliaksei Kasiyanik
 * date: 20.05.14
 */
public class SimpleProblemWithClarifying extends SteadyingEquation {

    public SimpleProblemWithClarifying(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public SimpleProblemWithClarifying(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    @Override
    public double[] f(double t, double[] arg) {
        double f0 = -0.0745 * arg[0];
        double f1 = -17829.0 * arg[1] - 100.0 * arg[2];
        double f2 = 100.0 * arg[1] - 17829.0 * arg[2];

        double[] res = {f0, f1, f2};
        return res;
    }

    @Override
    public double[][] JacobiMatrix(double[] arg) {
        return new double[0][];
    }

    @Override
    public double getJacobiMatrMaxEigenvalue(double t, double[] arg) {
        double[] m = new double[8];
        m[0] = 1.17 + 0.43 + 8.32;
        m[1] = 1.71 + 8.75;
        m[2] = 10.03 + 0.43 + 0.035;
        m[3] = 8.32 + 1.17 + 1.12;
        m[4] = 1.745 + 0.43 + 0.43;
        m[5] = 0.69 + 1.71 + 0.43 + 280. * Math.abs(arg[7]) + 0.69 + 280. * Math.abs(arg[5]);
        m[6] = 280. * arg[7] + 1.81 + 280. * arg[5];
        m[7] = 280. * arg[7] + 1.81 + 280. * arg[5];

        double max_j = 0.0;
        for (int i = 0; i < 8; i++)
            if (m[i] > max_j)
                max_j = m[i];

        return max_j;

//             return -17.1078;
    }

    @Override
    public double getAMatrMaxEigenvalue() {
        return 0.3333333;
    }

    @Override
    public double[][] getG() {
        return new double[0][];
    }
}
