package com.akasiyanik.bsu.coursework.methods;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import java.util.Locale;


/**
 * @author akasiyanik
 */
public class SteadyingProcessWithClarification extends SteadyingProcess {

    public final static int CLARIFICATION_ITERATION = 1;
    public SteadyingProcessWithClarification(AuxRungeKuttaMethod auxMethod, double EPS, double w, SteadyingEquation steadyingEquation) {
        super(auxMethod, EPS, w, steadyingEquation);
    }

    @Override
    public double[] getY() {
        int s = steadyingEquation.getS();
        int n = steadyingEquation.getN();
        double[] y0 = steadyingEquation.getY0();
        int dim = s * n;

        //составляем начальный вектор Y, состоящий из s векторов y0, размерности n
        double[] Y = new double[dim];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < n; j++) {
                Y[i * n + j] = y0[j];
            }
        }
        double err;
        iterationCount = 0;
        iterationCountWithClarifying = 0;
        do {
            double[] Y1;
            if (iterationCount % CLARIFICATION_ITERATION == 0) {
//                System.out.println("CLARIFICATION");
                Y1 = doSteadyingIterationWithClarifying(Y);
            } else {
                Y1 = doSimpleSteadyingIteration(Y);
            }
            //SOUT: ||r||
//            System.out.println(MatrixUtils.maxComponent(steadyingEquation.r(Y)));
//            System.out.printf(Locale.ENGLISH, "%.9f,", MatrixUtils.maxComponent(steadyingEquation.r(Y1)));
//            System.out.println();
            err = calculateErr(Y, Y1);
//            System.out.println("POST In method Err = " + err);
            Y = Y1;
            iterationCount++;
        } while (err >= EPS);
//        System.out.println("iteration count into st. process - " + iterationCount);

        return Y;
    }



    private double[] doSteadyingIterationWithClarifying(double[] Y) {
        double[] Y1 = F(Y);
//        System.out.println("---------------- Y1 ----------------");
//        printArray(Y1);
//        double err = calculateErr(Y, Y1);
//        System.out.println("\nPRE In method Err = " + err);
//        if (checkingConvergence()) {
            return clarify(Y, Y1);
//        } else {
//            return Y1;
//        }
    }

    private double[] doSimpleSteadyingIteration(double[] Y) {
       return F(Y);
    }

    private boolean checkingConvergence() { // проверка сходимости
        // todo akasiyanik
        return true;
    }

    private double[] clarify(double[] Y, double[] Y1) {
        double[] r_k = steadyingEquation.r(Y);
        double[] r_k1 = steadyingEquation.r(Y1);
        double lambda_m = getEigenvalueApproximation(r_k, r_k1);
        if (lambda_m > 0 || lambda_m < -1) {
            iterationCountWithClarifying++;
            Y1 = MatrixUtils.subtract(Y1, MatrixUtils.scalarMultipy(1.0 / lambda_m, r_k1));
        }
        return Y1;
    }

    //среднее значение центрального квантиля

    private double getEigenvalueApproximation(double[] r_k, double[] r_k1) { //getting lambda_m
        double derRin0 = auxMethod.getDerivativeRin0();
        double[] v = MatrixUtils.perComponentDivisionWithZerosInFraction(r_k1, r_k);
        //todo akasiyanik instead of average it is required to use average of central quantile
        double t = MatrixUtils.getMedianWithoutZeros(v);
//        System.out.println("!!!median: " + t);
        double eigenvalueApproximation = (t - 1.0) / (w * derRin0);
//        System.out.println("!!!eigenvalueApproximation: " + eigenvalueApproximation);

        return eigenvalueApproximation;
    }

    private int getMostHeavyComponentNumber(double[] Y, double[] Y1) {
        double maxDiff = 0.0;
        int mostHeavyComponent = 0;

        for (int i = 0; i < Y.length; i++) {
            if (Math.abs(Y1[i] - Y[i]) >= maxDiff) {
                maxDiff = Math.abs(Y1[i] - Y[i]);
                mostHeavyComponent = i;
            }
        }
        return mostHeavyComponent;
    }

}
