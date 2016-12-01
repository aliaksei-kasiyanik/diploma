package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.equations.LinearSteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author: akasiyanik
 */
public class PoissonProblem extends LinearSteadyingEquation {

    public PoissonProblem(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public PoissonProblem(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    //todo
//    @Override
//    public double getW() {
//        return 1.0 / 79.0826;
//    }

    @Override
    public double[] f(double t, double[] arg) {
//        int n = 30;
//        double[] res = new double[30];
//        double[] res = new double[n];
//        double tmp = (n + 1) * (n + 1);
//        for (int i = 0; i < n; i++) {
//            res[i] += arg[i] * tmp * (- 2);
//            if (i < n - 1) {
//                res[i] += arg[i + 1] * tmp;
//            }
//            if (i > 0) {
//                res[i] += arg[i - 1] * tmp;
//            }
//        }
//        return res;
        double[][] J = getJMatrix();
        return MatrixUtils.multiply(J, arg);
    }

    @Override
    public double[][] JacobiMatrix(double[] arg) {
        return new double[0][];
    }

    @Override
    public double getJacobiMatrMaxEigenvalue(double t, double[] arg) {
//        int n = 30;
        List<Double> evs = new ArrayList<Double>();
        for (int i = 0; i < n; i++) {
            double e = -4 * Math.pow(n - 1.0, 2.0) * Math.pow(Math.sin(Math.PI * i / (2 * n + 2)), 2);
            evs.add(Math.abs(e));
        }
        return Collections.max(evs);

//             return 17.1078;
    }

    @Override
    public double getAMatrMaxEigenvalue() {
        return 0.408248;
    }

    @Override
    public double[][] getJMatrix() {
        double[][] J = new double[n][n];
        double tmp = (n + 1) * (n + 1) ;
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

//    @Override
//    public double[][] getJMatrix() {
//        double[][] J = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            J[i][i] = -1922.0;
//            if (i < n - 1) {
//                J[i][i + 1] = 961.0;
//            }
//            if (i > 0) {
//                J[i][i - 1] = 961.0;
//            }
//        }
//        return J;
//    }
}
