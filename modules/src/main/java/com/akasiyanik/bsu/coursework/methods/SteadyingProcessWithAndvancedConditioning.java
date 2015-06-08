package com.akasiyanik.bsu.coursework.methods;

import com.akasiyanik.bsu.coursework.equations.LinearSteadyingEquation;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.add;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.multiply;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.scalarMultipy;

/**
 * @author akasiyanik
 */
public class SteadyingProcessWithAndvancedConditioning extends SteadyingProcess {

    private static int matrixVectorMultimplicationCount;

    protected double[] opressorCoefficients;
    protected double j_eighen;

    public SteadyingProcessWithAndvancedConditioning(AuxRungeKuttaMethod auxMethod, double EPS, double w, SteadyingEquation steadyingEquation, double[] opressorCoefficients, double j_eig) {
        super(auxMethod, EPS, w, steadyingEquation);
        this.opressorCoefficients = opressorCoefficients;
        this.j_eighen = j_eig;
    }

    protected double[] F(double[] Y0) {
        int n = auxMethod.getN();
        double[][] betta = auxMethod.getBetta();
        double[] Y = Y0;
        double gamma = auxMethod.getGamma();
//        double newW = - 1.0 / 176.0 ; //Poisson
//        double newW = - 1.0 / 1000.0;

        if (auxMethod.getSigma() % 2 == 0) { //четное сигма
            for (int m = 0; m < n; m++) {
                double[] r_after_cond = r_with_cond(Y);
                double[] sum1 = scalarMultipy(betta[m][0], r_after_cond);
                double[] wr = scalarMultipy(w, r_after_cond);
                double[] sum2 = scalarMultipy(betta[m][1], r_with_cond(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
        } else { // нечетное сигма
            for (int m = 0; m < n - 1; m++) {  // n-1 двухстадийные методы
                double[] r_after_cond = r_with_cond(Y);
                double[] sum1 = scalarMultipy(betta[m][0], r_after_cond);
                double[] wr = scalarMultipy(w, r_after_cond);
                double[] sum2 = scalarMultipy(betta[m][1], r_with_cond(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
            Y = add(Y, scalarMultipy(w * gamma, r_with_cond(Y)));   //n-ый метод - метод Эйлера
        }
        return Y;
    }

//    public double[] r_with_opressor(double[] Y) {   //dim Y: s*n x 1
//        double[] r = steadyingEquation.r(Y);
//        double[][] G = scalarMultipy(steadyingEquation.getW(), steadyingEquation.getG());
//        return opress(G, r, opressorCoefficients);
//    }

    public double[] r_with_cond(double[] Y) {   //dim Y: s*n x 1
        double[] r = steadyingEquation.r(Y);
//        double[][] J = scalarMultipy(steadyingEquation.getW(), ((LinearSteadyingEquation)steadyingEquation).getJMatrix());
        double[][] J = scalarMultipy(1.0 / j_eighen, ((LinearSteadyingEquation) steadyingEquation).getJMatrix());
        return advancedOpress(steadyingEquation.getS(), J, r, opressorCoefficients);
    }

//    public static double[] opress(double[][] A, double[] r0, double[] coeffs) {
//        int coeffNumber = coeffs.length;
//        double[] res = new double[r0.length];
//        double[] tmp = new double[r0.length];
//        for (int i = coeffNumber - 1; i >= 0; i--) {
//            tmp = add(res, scalarMultipy(coeffs[i], r0));
//            res = multiply(A, tmp);
//            matrixVectorMultimplicationCount++;
//        }
//        return tmp;
//    }

    public static double[] advancedOpress(int s, double[][] A, double[] r0, double[] coeffs) {
        int coeffNumber = coeffs.length;
        double[][] commonRes = new double[s][];
        double[][] splittedR0 = split(r0, s);
        for (int j = 0; j < s; j++) {
            double[] tmp = new double[r0.length / s];
            double[] res = new double[r0.length / s];
            for (int i = coeffNumber - 1; i >= 0; i--) {
                tmp = add(res, scalarMultipy(coeffs[i], splittedR0[j]));
                res = multiply(A, tmp);
                matrixVectorMultimplicationCount++;
            }
            commonRes[j] = tmp;
        }
        return flatten(commonRes);
    }

    public static int getMatrixVectorMultimplicationCount() {
        return matrixVectorMultimplicationCount;
    }

    public static void resetMatrixVectorMultimplicationCount() {
        matrixVectorMultimplicationCount = 0;
    }

    public static class PowerMethod {

        private double[][] A;

        private double[] y0;

        private double eps;

        private double[] coeffs;

        private int s;

        public PowerMethod(int s, double[][] A, double[] y0, double[] coeffs,  double eps) {
            this.A = A;
            this.y0 = y0;
            this.coeffs = coeffs;
            this.eps = eps;
            this.s = s;
        }

        public double solve() {
            double[] u = y0;
            double[] v;
            double lambda = 100.0;
            double new_lambda = -100.0;
            while (Math.abs(new_lambda - lambda) > eps) {
//                v = cond(A, multiply(A, u));
                v = multiply(A, u);
                lambda = new_lambda;
                new_lambda = Math.sqrt(maxComponent(v));
                u = divide(v, new_lambda);
//                System.out.println(new_lambda);
            }
            return new_lambda;
        }

        public double solve(double[][] J) {
            double[] u = y0;
            double[] v;
            double lambda = 100.0;
            double new_lambda = -100.0;
            while (Math.abs(new_lambda - lambda) > eps) {
                v = cond(J, multiply(A, u));
//                v = multiply(A, u);
                lambda = new_lambda;
                new_lambda = Math.sqrt(maxComponent(v));
                u = divide(v, new_lambda);
//                System.out.println(new_lambda);
            }
            return new_lambda;

        }

//        public double solve(double[][] G) {
//            int i = 1;
//            double[] y_k = y0;
//            //A = G_with_w
//            double[] y_k_1 = cond(A, multiply(G, y_k));
//            double[] y_k_2 = cond(A, multiply(G, y_k_1));
//            double[] y_k_3 = cond(A, multiply(G, y_k_2));
//            double r = r(y_k, y_k_1, y_k_2, y_k_3);
//            boolean not_first_iteration = false;
//            double new_r = 0.0;
//            do {
//                if (not_first_iteration) {
//                    r = new_r;
//                }
//                y_k = y_k_1;
//                y_k_1 = cond(A, multiply(G, y_k));
//                y_k_2 = cond(A, multiply(G, y_k_1));
//                y_k_3 = cond(A, multiply(G, y_k_2));
//                new_r = r(y_k, y_k_1, y_k_2, y_k_3);
////            System.out.println("iter = " + i + "; new_r = " + new_r);
//                not_first_iteration = true;
//                i++;
//            } while (Math.abs(new_r - r) > eps);
////        System.out.println("power method iter count : " + i);
////        System.out.println("power method error : " + Math.abs(new_r - r));
////        System.out.println("power method r^2 = " + r);
//            return r;
//        }

        private double r(double[] y_k__1, double[] y_k, double[] y_k_1, double[] y_k_2) {
            double y__1 = y_k__1[0];
            double y_0 = y_k[0];
            double y_1 = y_k_1[0];
            double y_2 = y_k_2[0];
            return (y_0 * y_2 - Math.pow(y_1, 2)) / (y__1 * y_1 - Math.pow(y_0, 2));
        }

        protected double[] cond(double[][] A, double[] r0) {
            int coeffNumber = coeffs.length;
            double[][] commonRes = new double[s][];
            double[][] splittedR0 = split(r0, s);
//            double[][][] A_blocks = getBlocksFromDiagonalMatrix(A, s);
            for (int j = 0; j < s; j++) {
                double[] tmp = new double[r0.length / s];
                double[] res = new double[r0.length / s];
                for (int i = coeffNumber - 1; i >= 0; i--) {
                    tmp = add(res, scalarMultipy(coeffs[i], splittedR0[j]));
                    res = multiply(A, tmp);
                    matrixVectorMultimplicationCount++;
                }
                commonRes[j] = tmp;
            }
            return flatten(commonRes);
        }


    }
}