package com.akasiyanik.bsu.coursework.methods;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

/**
 * @author: akasiyanik
 */
public class SteadyingProcess {

    protected static double EPS;
    protected SteadyingEquation steadyingEquation;
    protected double w;
    protected AuxRungeKuttaMethod auxMethod;

    protected int iterationCount = 0;
    protected int iterationCountWithClarifying = 0;

    public SteadyingProcess(AuxRungeKuttaMethod auxMethod, double EPS, double w, SteadyingEquation steadyingEquation) {
        this.auxMethod = auxMethod;
        this.EPS = EPS;
        this.w = w;
        this.steadyingEquation = steadyingEquation;
    }

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
            double[] Y1 = doSimpleSteadyingIteration(Y);
            err = calculateErr(Y, Y1);
//            System.out.println("POST In method Err = " + err);
            Y = Y1;
            iterationCount++;
        } while (err >= EPS);
//        System.out.println("iteration count into st. process - " + iterationCount);

        return Y;
    }

    private double[] doSimpleSteadyingIteration(double[] Y) {
        return F(Y);
    }


    protected double calculateErr(double[] arr1, double[] arr2) {
        double max_err = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            double err = Math.abs(arr1[i] - arr2[i]);
            if (err > max_err) {
                max_err = err;
            }
        }
        return max_err;
    }

    protected double[] F(double[] Y0) {
        int n = auxMethod.getN();
        double[][] betta = auxMethod.getBetta();
        double[] Y = Y0;
        double gamma = auxMethod.getGamma();

        if (auxMethod.getSigma() % 2 == 0) { //четное сигма
            for (int m = 0; m < n; m++) {
                double[] sum1 = scalarMultipy(betta[m][0], steadyingEquation.r(Y));
                double[] wr = scalarMultipy(w, steadyingEquation.r(Y));
                double[] sum2 = scalarMultipy(betta[m][1], steadyingEquation.r(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
        } else { // нечетное сигма
            for (int m = 0; m < n - 1; m++) {  // n-1 двухстадийные методы
                double[] sum1 = scalarMultipy(betta[m][0], steadyingEquation.r(Y));
                double[] wr = scalarMultipy(w, steadyingEquation.r(Y));
                double[] sum2 = scalarMultipy(betta[m][1], steadyingEquation.r(add(Y, wr)));
                Y = add(Y, scalarMultipy(w, add(sum1, sum2)));
            }
            Y = add(Y, scalarMultipy(w * gamma, steadyingEquation.r(Y)));   //n-ый метод - метод Эйлера
        }
        return Y;
    }

    public int getIterationCount() {
        return iterationCount;
    }

    public void setIterationCount(int iterationCount) {
        this.iterationCount = iterationCount;
    }

    public int getIterationCountWithClarifying() {
        return iterationCountWithClarifying;
    }

    public void setIterationCountWithClarifying(int iterationCountWithClarifying) {
        this.iterationCountWithClarifying = iterationCountWithClarifying;
    }

//    private double[][] calculateK(double[] Y) {
//        int sigma = auxMethod.getS();
//        int n = steadyingEquation.getN();
//        int s = steadyingEquation.getS();
//        double[][] a = auxMethod.getA();
//        double[][] K = new double[sigma][s * n];
//        for (int p = 0; p < sigma; p++) {
//            double[] sum = new double[s * n];
//            for (int q = 0; q <= p - 1; q++) {
//                sum = add(sum, scalarMultipy(a[p][q], K[q]));
//            }
//            K[p] = steadyingEquation.r(add(Y, scalarMultipy(w, sum)));
//        }
//        return K;
//    }


}
