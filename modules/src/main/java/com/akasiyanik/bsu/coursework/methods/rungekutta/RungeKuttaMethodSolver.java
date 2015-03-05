package com.akasiyanik.bsu.coursework.methods.rungekutta;

/**
 * @author: akasiyanik
 */
public class RungeKuttaMethodSolver {

    private RungeKuttaMethod method;
    double x0;
    double y0;
    double tau;

    public double solve() {
        return y1();
    }

    private double y1() {
        int s = method.getS();
        double[] b = method.getB();
        double[] k = calculateK();

        double sum = 0.0;
        for (int i = 0; i < s; i++) {
            sum += b[i] * k[i];
        }
        return y0 + tau * sum;
    }

    private double[] calculateK() {
        int s = method.getS();
        double[] c = method.getC();
        double[][] a = method.getA();

        double[] k = new double[s];

        for (int i = 0; i < s; i++) {
            double sum = 0.0;
            for (int j = 0; j < i - 1; j++) {
                sum += a[i][j]*k[j];
            }
            double arg1 = x0 + c[i] * tau;
            double arg2 = y0 + tau * sum;
            k[i] = f(arg1, arg2);

        }

        return k;

    }

    private double f(double arg1, double arg2) {
        // type your function
        return arg1 + arg2;
    }


}
