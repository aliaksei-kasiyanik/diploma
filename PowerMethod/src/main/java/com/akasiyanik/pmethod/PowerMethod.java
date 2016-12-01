package com.akasiyanik.pmethod;

import java.util.LinkedList;
import java.util.List;

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

//     public double solve() {
//         int i = 1;
//         double[] y_k = y0;
//         double[] y_k_1 = multiply(A, y_k);
//         double[] y_k_2 = multiply(A, y_k_1);
//         double[] y_k_3 = multiply(A, y_k_2);
//         double[] r = r(y_k, y_k_1, y_k_2, y_k_3);
//         boolean not_first_iteration = false;
//         double[] new_r = new double[r.length];
//         do {
//             if (not_first_iteration) {
//                 r = new_r;
//             }
//             y_k = divide(y_k_3, norm(y_k_3));
//             y_k_1 = multiply(A, y_k);
//             y_k_2 = multiply(A, y_k_1);
//             y_k_3 = multiply(A, y_k_2);
//             new_r = r(y_k, y_k_1, y_k_2, y_k_3);
//             System.out.println("iter = " + i + "; new_r = "); output(new_r);
//             not_first_iteration = true;
//             i++;
//         } while (calculateMaxError(new_r, r) > eps);
//         System.out.println("power method iter count : " + i);
//         System.out.println("power method error : " + calculateMaxError(new_r, r));
//         System.out.println("power method r^2 = "); output(r);
//         return maxComponent(r);
//     }

//    public double solve() {
//        double[] u = y0;
//        double[] v;
//        double lambda = 100.0;
//        double new_lambda = -100.0;
//        while (Math.abs(new_lambda - lambda) > eps) {
//            v = multiply(A, u);
//            lambda = new_lambda;
//            new_lambda = Math.sqrt(maxComponent(v));
//            u = divide(v, new_lambda);
//            System.out.println(new_lambda);
//        }
//        return new_lambda;
//
//    }

//    public double solve() {
//        double[] u = y0;
//        double[] v;
//        double lambda = 100.0;
//        double new_lambda = -100.0;
//        while (Math.abs(lambda - new_lambda) > eps) {
//            lambda = new_lambda;
//            v = multiply(A, u);
//            u = multiply(A, v);
//            new_lambda = Math.sqrt(maxComponent(u));
//            System.out.println(new_lambda);
//            u = add(scalarMultipy(new_lambda, v), u);
//            u = divide(u, maxComponent(u));
////            v = divide(v, maxComponent(v));
//        }
//        return new_lambda;
//
//    }

    public double solve() {
         int i = 1;
         double[] y_k = y0;
         double[] y_k_1 = multiply(A, y_k);
         double[] y_k_2 = multiply(A, y_k_1);
         double[] y_k_3 = multiply(A, y_k_2);
         double[] r = r(y_k, y_k_1, y_k_2, y_k_3);
         double lambda = calculateRealEigenvalue(y_k_1, y_k_3);
        //---- Norm --
//        y_k_1 = divide(y_k_1, getNorm(y_k_1));
//        y_k_2 = divide(y_k_2, getNorm(y_k_2));
//        y_k_3 = divide(y_k_3, getNorm(y_k_3));
        //------------

         boolean not_first_iteration = false;
         double[] new_r = new double[r.length];
         double new_lambda = 0.0;
         do {
             if (not_first_iteration) {
                 r = new_r;
                 lambda = new_lambda;
             }
             y_k = divide(y_k_3, getNorm(y_k_3));
             y_k_1 = multiply(A, y_k);
             y_k_2 = multiply(A, y_k_1);
             y_k_3 = multiply(A, y_k_2);
             new_r = r(y_k, y_k_1, y_k_2, y_k_3);
             new_lambda = calculateRealEigenvalue(y_k_1, y_k_3);
             //---- Norm --
//             y_k_3 = divide(y_k_3, getNorm(y_k_3));
             //------------
             System.out.println("iter = " + i + "; new_r = "); output(new_r);
             System.out.println("iter = " + i + "; new_lambda = " + new_lambda);
             not_first_iteration = true;
             i++;
         } while (checkComplexConvergence(new_r, r) || checkRealConvergence(lambda, new_lambda));
         System.out.println("power method iter count : " + i);
         System.out.println("power method error : " + calculateMaxError(new_r, r));
         System.out.println("power method r^2 = "); output(r);

        double complexEigenvalue = Math.sqrt(maxComponent(r));
        double realEigenvalue = new_lambda;
        System.out.println("Complex: " + complexEigenvalue);
        System.out.println("Real: " + realEigenvalue);
        return chooseEigenvalue(y_k_3, complexEigenvalue, realEigenvalue);
     }

    private boolean checkRealConvergence(double lambda, double new_lambda) {
        double tmp = Math.abs(Math.abs(lambda) - Math.abs(new_lambda));
        if (Double.isNaN(tmp) || Double.isInfinite(tmp)) {
            return false;
        } else {
            return Math.abs((Math.abs(lambda) - Math.abs(new_lambda))) > eps;
        }
    }

    private boolean checkComplexConvergence(double[] new_r, double[] r) {
        double tmp = calculateMaxError(new_r, r);
        if (Double.isNaN(tmp) || Double.isInfinite(tmp)) {
            return false;
        } else {
            return tmp > eps;
        }
    }

    private double chooseEigenvalue(double[] y_k_3, double complexEigenvalue, double realEigenvalue) {
//        evector = zkkk/Norm@zkkk;
        if (Double.isNaN(realEigenvalue)) return complexEigenvalue;

        double[] evector = divide(y_k_3, getNorm(y_k_3));
        double[] tmp1 = subtract(multiply(A, evector), scalarMultipy(realEigenvalue, evector));
        double[] tmp2 = add(multiply(A, evector), scalarMultipy(realEigenvalue, evector));
        if (maxComponent(tmp1) < 1.0 || maxComponent(tmp2) < 1.0) {
            return realEigenvalue;
        } else {
            return complexEigenvalue;
        }
    }

    private double calculateRealEigenvalue(double[] y_k_1, double[] y_k_3) {
//        real = (Abs@(zkkk/zk))^0.5;
        return Math.sqrt(getAverage(divide(y_k_3, y_k_1)));
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

    private double[] normalize(double[] y) {
        double max = maxComponent(y);
        return divide(y, max);
    }

    private double getNorm(double[] y) {
        return maxComponent(y);
    }
}
