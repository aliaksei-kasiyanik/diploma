package com.akasiyanik.bsu.coursework;

import com.akasiyanik.bsu.coursework.equations.LinearSteadyingEquation;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithAndvancedConditioning;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithOpression;
import com.akasiyanik.bsu.coursework.methods.power.AdvPowerMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;
import com.akasiyanik.bsu.coursework.problems.PoissonProblem;
import com.akasiyanik.bsu.coursework.utils.FileUtils;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import java.io.*;

/**
 * @author akasiyanik
 */
public class PoissonSolverWithAdvancedOpression implements Solver {

    private int iterationCount = 0;
    private int iterationCountWithClarifying = 0;
    private static double eps = Math.pow(10, -6);

    public static void main(String[] args) throws FileNotFoundException {
        PoissonSolverWithAdvancedOpression solver = new PoissonSolverWithAdvancedOpression();
        double tau = 0.05;
        int n = 30;

        double t0 = 0.0;

        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0;
        }
        solver.solve(t0, y0, tau, eps);
    }

    public double[] solve(double t0, double[] y0, double tau, double eps) {
        RungeKuttaMethod baseRungeKuttaMethod = null;
        AuxRungeKuttaMethod auxRungeKuttaMethod = null;
        double[] opressionCoeffs = null;

        InputStream base = null;
        InputStream aux = null;
        InputStream coeff = null;
        try {
            base = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\RadauIIA-5-Order-Method.txt"));
            aux = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\9-Order-Generated-AuxMethod.txt"));
            coeff = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\opression_coefficients.txt"));

            opressionCoeffs = FileUtils.readOpressionCoefficients(coeff);
            baseRungeKuttaMethod = FileUtils.readBaseRungeKuttaMethod(base);
            auxRungeKuttaMethod = FileUtils.readAuxRungeKuttaMethod(aux);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } finally {
            try {
                if (base != null) {
                    base.close();
                }
                if (aux != null) {
                    aux.close();
                }
                if (coeff != null) {
                    coeff.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

        LinearSteadyingEquation steadyingEquation = new PoissonProblem(tau, t0, y0, baseRungeKuttaMethod);
        double JEigenvalue =  steadyingEquation.getJacobiMatrMaxEigenvalue(0.0, null);
//        double JEigenvalue = new AdvPowerMethod(steadyingEquation.getS(), steadyingEquation.getJMatrix(), steadyingEquation.getY0(), opressionCoeffs, eps).solve();
        System.out.println("JEigenvalue " + JEigenvalue);


        double w = - 0.7 / getNewW(steadyingEquation, opressionCoeffs, JEigenvalue);
        System.out.println("newW " + w);
        setIterationCount(0);

        SteadyingProcessWithAndvancedConditioning steadyingProcess = new SteadyingProcessWithAndvancedConditioning(auxRungeKuttaMethod, eps, w, steadyingEquation, opressionCoeffs, JEigenvalue);
        double[] Y = steadyingProcess.getY();

        setIterationCount(steadyingProcess.getIterationCount());
        setIterationCountWithClarifying(steadyingProcess.getIterationCountWithClarifying());

        double[] solution = steadyingEquation.getSolution(Y);

//        System.out.println("With clarifying: " + getIterationCountWithClarifying());
//
//        System.out.println("Y:");
//        for (int i = 0; i < Y.length; i++) {
//            System.out.println(Y[i]);
//        }

        System.out.println("MATRIX x VECTOR count - " + SteadyingProcessWithOpression.getMatrixVectorMultimplicationCount());

        System.out.println("Solution:");
        for (int i = 0; i < solution.length; i++) {
            System.out.println(solution[i]);
        }

        return solution;
    }

    private double getNewW(LinearSteadyingEquation steadEquation, double[] opressCoeffs, double j_eighen) {
        double[][] J = steadEquation.getJMatrix();
        double eps = Math.pow(10, -3);
        double[][] J_norm = MatrixUtils.scalarMultipy(1.0 / j_eighen, J);
//        double[][] G = MatrixUtils.scalarMultipy(steadEquation.getW(), steadEquation.getG());
//        return new AdvPowerMethod(steadEquation.getS(), G, getY0(steadEquation), opressCoeffs, eps).solve(J_norm);
        return new AdvPowerMethod(steadEquation.getS(), steadEquation.getG(), getY0(steadEquation), opressCoeffs, eps).solve(J_norm);

    }

    public double[] getY0(SteadyingEquation steadyingEquation) {
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
}

