package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

/**
 * @author: akasiyanik
 */
public class PuassonProblem extends SteadyingEquation {

    public PuassonProblem(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public PuassonProblem(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    @Override
        public double[] f(double t, double[] arg) {
            double f0 = -1.71 * arg[0] + 0.43 * arg[1] + 8.32 * arg[2] + 0.0007;
            double f1 = 1.71 * arg[0] - 8.75 * arg[1];
            double f2 = -10.03 * arg[2] + 0.43 * arg[3] + 0.035 * arg[4];
            double f3 = 8.32 * arg[1] + 1.71 * arg[2] - 1.12 * arg[3];
            double f4 = -1.745 * arg[4] + 0.43 * arg[5] + 0.43 * arg[6];
            double f5 = -280 * arg[5] * arg[7] + 0.69 * arg[3] + 1.71 * arg[4] - 0.43 * arg[5] + 0.69 * arg[6];
            double f6 = 280 * arg[5] * arg[7] - 1.81 * arg[6];
            double f7 = -280 * arg[5] * arg[7] + 1.81 * arg[6];

            double[] res = {f0, f1, f2, f3, f4, f5, f6, f7};
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
}
