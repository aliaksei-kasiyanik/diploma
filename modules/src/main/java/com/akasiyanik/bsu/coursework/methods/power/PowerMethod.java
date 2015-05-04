package com.akasiyanik.bsu.coursework.methods.power;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;

/**
 * @author Aliaksei Kasiyanik
 */
public class PowerMethod {

    private double[][] A;

    private double[] y0;

    private double eps;

    public PowerMethod(double[][] A, double[] y0, double eps) {
        this.A = A;
        this.y0 = y0;
        this.eps = eps;
    }

    public double solve() {
        double[] y_k = y0;
        double[] y_k_1 = multiply(A, y_k);
        double[] y_k_2 = multiply(A, y_k_1);
        double[] y_k_3 = multiply(A, y_k_2);
        double r = r(y_k, y_k_1, y_k_2, y_k_3);
        double new_r;
        do {
            y_k = y_k_1;
            y_k_1 = multiply(A, y_k);
            y_k_2 = multiply(A, y_k_1);
            y_k_3 = multiply(A, y_k_2);
            new_r = r(y_k, y_k_1, y_k_2, y_k_3);
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
