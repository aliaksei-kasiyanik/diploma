package com.akasiyanik.bsu.coursework.methods.power;

import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import static com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithOpression.opress;
/**
 * @author Aliaksei Kasiyanik
 */
public  class PowerMethod {

    private double[][] A;

    private double[] y0;

    private double eps;

    private double[] coeffs;

    public PowerMethod(double[][] A, double[] y0, double[] coeffs,  double eps) {
        this.A = A;
        this.y0 = y0;
        this.coeffs = coeffs;
        this.eps = eps;
    }

    public double solve(double[][] G) {
        double[] y_k = y0;
        double[] y_k_1 = opress(A, MatrixUtils.multiply(G, y_k), coeffs);
        double[] y_k_2 = opress(A, MatrixUtils.multiply(G, y_k_1), coeffs);
        double[] y_k_3 = opress(A, MatrixUtils.multiply(G, y_k_2), coeffs);
        double r = r(y_k, y_k_1, y_k_2, y_k_3);
        boolean not_first_iteration = false;
        double new_r = 0.0;
        do {
            if (not_first_iteration) {
                r = new_r;
            }
            y_k = y_k_1;
            y_k_1 = opress(A, MatrixUtils.multiply(G, y_k), coeffs);
            y_k_2 = opress(A, MatrixUtils.multiply(G, y_k_1), coeffs);
            y_k_3 = opress(A, MatrixUtils.multiply(G, y_k_2), coeffs);
            new_r = r(y_k, y_k_1, y_k_2, y_k_3);
            not_first_iteration = true;
        } while (Math.abs(new_r - r) > eps);
        return r;
    }

    private double r(double[] y_k__1, double[] y_k, double[] y_k_1, double[] y_k_2) {
        double y__1 = y_k__1[0];
        double y_0 = y_k[0];
        double y_1 = y_k_1[0];
        double y_2 = y_k_2[0];
        return (y_0 * y_2 - Math.pow(y_1, 2)) / (y__1 * y_1 - Math.pow(y_0, 2));
    }
}

