package com.akasiyanik.bsu.coursework.methods;

import com.akasiyanik.bsu.coursework.equations.LinearSteadyingEquation;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.kroneckerProduct;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.scalarMultipy;

/**
 * @author akasiyanik
 */
public class SteadyingProcessWithAndvancedConditioningA extends SteadyingProcess {

    private static int matrixVectorMultimplicationCount;

    protected double[] opressorCoefficients;
    protected double j_eighen;

    protected double[][] invA_I;
    protected double[][] invA;

    public SteadyingProcessWithAndvancedConditioningA(AuxRungeKuttaMethod auxMethod, double EPS, double w, SteadyingEquation steadyingEquation, double[] opressorCoefficients, double j_eig) {
        super(auxMethod, EPS, w, steadyingEquation);
        this.opressorCoefficients = opressorCoefficients;
        this.j_eighen = j_eig;
        this.invA = steadyingEquation.getImplicitRungeKuttaMethod().getInverseA();
        this.invA_I = kroneckerProduct(invA, getIdentityMatrix(steadyingEquation.getN()));
    }

    protected double[] F(double[] Y0) {
        int n = auxMethod.getN();
        double[][] betta = auxMethod.getBetta();
        double[] Y = Y0;
        double gamma = auxMethod.getGamma();

        if (auxMethod.getSigma() % 2 == 0) { //четное сигма
            for (int m = 0; m < n; m++) {
                double[] sum1 = scalarMultipy(betta[m][0], r_cond(Y));
                double[] wr = scalarMultipy(w, r_cond(Y));
                double[] sum2 = scalarMultipy(betta[m][1], r_cond(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
        } else { // нечетное сигма
            for (int m = 0; m < n - 1; m++) {  // n-1 двухстадийные методы
                double[] sum1 = scalarMultipy(betta[m][0], r_cond(Y));
                double[] wr = scalarMultipy(w, r_cond(Y));
                double[] sum2 = scalarMultipy(betta[m][1], r_cond(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
            Y = add(Y, scalarMultipy(w * gamma, r_cond(Y)));   //n-ый метод - метод Эйлера
        }
        return Y;
    }

//    public double[] r_with_opressor(double[] Y) {   //dim Y: s*n x 1
//        double[] r = steadyingEquation.r(Y);
//        double[][] G = scalarMultipy(steadyingEquation.getW(), steadyingEquation.getG());
//        return opress(G, r, opressorCoefficients);
//    }

    public double[] r_cond(double[] Y) {   //dim Y: s*n x 1
        double[] e = getIdentityVector(steadyingEquation.getS()); //dim: s x 1
        double[] tmp1 = multiply(invA_I, kroneckerProduct(e, steadyingEquation.getY0()));
        double[] tmp2 = scalarMultipy(steadyingEquation.getTau(), steadyingEquation.funct_F(steadyingEquation.getT0(), Y));
        double[] tmp3 = multiply(invA_I, Y);
        return subtract(add(tmp1, tmp2), tmp3);
    }



//    public double[] r(double[] Y) {   //dim Y: s*n x 1
//        double[][] tmp1 = kroneckerProduct(getIdentityMatrix(steadyingEquation.getS()), ((LinearSteadyingEquation) steadyingEquation).getJMatrix());
//        double[][] G = subtract(scalarMultipy(steadyingEquation.getTau(), tmp1), invA_I);
//        double[] tmp2 = multiply(((LinearSteadyingEquation) steadyingEquation).getJMatrix(), steadyingEquation.getY0());
//        double[] g = multiply(invA_I, kroneckerProduct(getIdentityMatrix(steadyingEquation.getS()), tmp2));
//        return add(multiply(G, Y), g);
//    }

//    public double[] r_with_cond(double[] Y) {   //dim Y: s*n x 1
//        double[] r = r(Y);
//        double[][] invA = steadyingEquation.getImplicitRungeKuttaMethod().getInverseA();
////        double[][] J = scalarMultipy(steadyingEquation.getW(), ((LinearSteadyingEquation)steadyingEquation).getJMatrix());
////        invA = scalarMultipy(w, invA);
//        invA = scalarMultipy(1.0 / j_eighen, invA);
//        return advancedOpress(steadyingEquation.getN(), invA, r);
//    }

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

//    public static double[] advancedOpress(int n, double[][] A, double[] r0) {
//        double[][] commonRes = new double[n][];
//        double[][] splittedR0 = split(r0, n);
//        for (int j = 0; j < n; j++) {
//            double[] res = multiply(A, splittedR0[j]);
//            commonRes[j] = res;
//        }
//        return flatten(commonRes);
//    }

    public static int getMatrixVectorMultimplicationCount() {
        return matrixVectorMultimplicationCount;
    }

    public static void resetMatrixVectorMultimplicationCount() {
        matrixVectorMultimplicationCount = 0;
    }

}