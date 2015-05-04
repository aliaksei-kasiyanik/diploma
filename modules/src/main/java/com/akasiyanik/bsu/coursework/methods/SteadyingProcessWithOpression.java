package com.akasiyanik.bsu.coursework.methods;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.power.PowerMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.add;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.multiply;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.scalarMultipy;

/**
 * @author akasiyanik
 */
public class SteadyingProcessWithOpression extends SteadyingProcess {

    protected double[] opressorCoefficients;

    public SteadyingProcessWithOpression(AuxRungeKuttaMethod auxMethod, double EPS, double w, SteadyingEquation steadyingEquation, double[] opressorCoefficients) {
        super(auxMethod, EPS, w, steadyingEquation);
        this.opressorCoefficients = opressorCoefficients;
    }

    protected double[] F(double[] Y0) {
        int n = auxMethod.getN();
        double[][] betta = auxMethod.getBetta();
        double[] Y = Y0;
        double gamma = auxMethod.getGamma();
        double newW = - 1.0 / 176.0 ;

//        double newW = - 1.0 / getNewW(Y0);

        if (auxMethod.getSigma() % 2 == 0) { //четное сигма
            for (int m = 0; m < n; m++) {
                double[] sum1 = scalarMultipy(betta[m][0], r_with_opressor(Y));
                double[] wr = scalarMultipy(newW, r_with_opressor(Y));
                double[] sum2 = scalarMultipy(betta[m][1], r_with_opressor(add(Y, wr)));
                Y = add(Y, scalarMultipy(newW, add(sum1, sum2)));
            }
        } else { // нечетное сигма
            for (int m = 0; m < n - 1; m++) {  // n-1 двухстадийные методы
                double[] sum1 = scalarMultipy(betta[m][0], r_with_opressor(Y));
                double[] wr = scalarMultipy(newW, r_with_opressor(Y));
                double[] sum2 = scalarMultipy(betta[m][1], r_with_opressor(add(Y, wr)));
                Y = add(Y, scalarMultipy(newW, add(sum1, sum2)));
            }
            Y = add(Y, scalarMultipy(newW * gamma, r_with_opressor(Y)));   //n-ый метод - метод Эйлера
        }
        return Y;
    }

    public double[] r_with_opressor(double[] Y) {   //dim Y: s*n x 1
        double[] r = steadyingEquation.r(Y);
        double[][] G = scalarMultipy(steadyingEquation.getW(), steadyingEquation.getG());
        return opress(G, r, opressorCoefficients);
    }

    protected static double[] opress(double[][] A, double[] r0, double[] coeffs) {
        int coeffNumber = coeffs.length;
        double[] res = new double[r0.length];
        double[] tmp = new double[r0.length];
        for (int i = coeffNumber - 1; i >= 0; i--) {
            tmp = add(res, scalarMultipy(coeffs[i], r0));
            res = multiply(A, tmp);
        }
        return tmp;
    }

    private double getNewW(double[] y0) {
        double[][] G = steadyingEquation.getG();
        double eps = Math.pow(10, -6);
        return Math.sqrt(new PowerMethod(G, y0, opressorCoefficients, eps).solve());
    }

    class PowerMethod {

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

        public double solve() {
            double[] y_k = y0;
            double[] y_k_1 = opress(A, y_k, coeffs);
            double[] y_k_2 = opress(A, y_k_1, coeffs);
            double[] y_k_3 = opress(A, y_k_2, coeffs);
            double r = r(y_k, y_k_1, y_k_2, y_k_3);
            double new_r;
            do {
                y_k = y_k_1;
                y_k_1 = opress(A, y_k, coeffs);
                y_k_2 = opress(A, y_k_1, coeffs);
                y_k_3 = opress(A, y_k_2, coeffs);
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
//    protected double[] opress(double[][] A, double[] r0, double[] coeffs) {
//        int coeffNumber = coeffs.length;
//        double[] res = scalarMultipy(coeffs[0], r0);
//        double[] r = r0;
//        for (int i = 1; i < coeffNumber; i++) {
//            r = multiply(A, r);
//            res = add(res, scalarMultipy(coeffs[i], r));
//        }
//        return res;
//    }


}
