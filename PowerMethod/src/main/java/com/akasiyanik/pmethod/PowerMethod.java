package com.akasiyanik.pmethod;

import static com.akasiyanik.utils.MatrixUtils.*;

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
         int i = 1;
         double[] y_k = y0;
         double[] y_k_1 = norm(multiply(A, y_k));
         double[] y_k_2 = norm(multiply(A, y_k_1));
         double[] y_k_3 = norm(multiply(A, y_k_2));
         double[] r = r(y_k, y_k_1, y_k_2, y_k_3);
         boolean not_first_iteration = false;
         double[] new_r = new double[r.length];
         do {
             if (not_first_iteration) {
                 r = new_r;
             }
             y_k = y_k_1;
             y_k_1 = norm(multiply(A, y_k));
             y_k_2 = norm(multiply(A, y_k_1));
             y_k_3 = norm(multiply(A, y_k_2));
             new_r = r(y_k, y_k_1, y_k_2, y_k_3);
             System.out.println("iter = " + i + "; new_r = "); output(new_r);
             not_first_iteration = true;
             i++;
         } while (calculateMaxError(new_r, r) > eps);
         System.out.println("power method iter count : " + i);
         System.out.println("power method error : " + calculateMaxError(new_r, r));
         System.out.println("power method r^2 = "); output(r);
         return maxComponent(r);
     }

    private double[] r(double[] y_k__1, double[] y_k, double[] y_k_1, double[] y_k_2) {
        double r[] = new double[y_k__1.length];
        for (int i = 0; i < y_k__1.length; i++) {
            double y__1 = y_k__1[i];
            double y_0 = y_k[i];
            double y_1 = y_k_1[i];
            double y_2 = y_k_2[i];
            if ( (y__1 * y_1 - Math.pow(y_0, 2)) == 0.0) {
                r[i] = 0.0;
            } else {
                r[i] = (y_0 * y_2 - Math.pow(y_1, 2)) / (y__1 * y_1 - Math.pow(y_0, 2));
            }
        }
        return r;
    }

    private double[] norm(double[] y) {
        double max = maxComponent(y);
        return scalarMultipy(1.0 / max, y);
    }
}
