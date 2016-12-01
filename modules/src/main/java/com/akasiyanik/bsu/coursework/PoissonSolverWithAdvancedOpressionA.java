package com.akasiyanik.bsu.coursework;

import com.akasiyanik.bsu.coursework.equations.LinearSteadyingEquation;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithAndvancedConditioning;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithAndvancedConditioningA;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithOpression;
import com.akasiyanik.bsu.coursework.methods.power.AdvPowerMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;
import com.akasiyanik.bsu.coursework.problems.PoissonProblem;
import com.akasiyanik.bsu.coursework.utils.FileUtils;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

/**
 * @author akasiyanik
 */
public class PoissonSolverWithAdvancedOpressionA implements Solver {

    private int iterationCount = 0;
    private int iterationCountWithClarifying = 0;
    private static double eps = Math.pow(10, -6);

    public static void main(String[] args) throws FileNotFoundException {
        PoissonSolverWithAdvancedOpressionA solver = new PoissonSolverWithAdvancedOpressionA();
        double tau = 0.05;
        int n = 20;

        double t0 = 0.0;

        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0 + i;
        }
        solver.solve(t0, y0, tau, eps);
    }

    public double[] solve(double t0, double[] y0, double tau, double eps) {
        RungeKuttaMethod baseRungeKuttaMethod = null;
        AuxRungeKuttaMethod auxRungeKuttaMethod = null;
        double[] opressionCoeffs = null;

        InputStream base = null;
        InputStream aux = null;
        try {
            base = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\RadauIIA-5-Order-Method.txt"));
            aux = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\9-Order-Generated-AuxMethod.txt"));

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
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

        PoissonProblem steadyingEquation = new PoissonProblem(tau, t0, y0, baseRungeKuttaMethod);
        double invAEigen = 4.0611980714;

        double w =  0.7 / getNewW(steadyingEquation, invAEigen);
        setIterationCount(0);

        System.out.println("w " + w);

        SteadyingProcessWithAndvancedConditioningA steadyingProcess = new SteadyingProcessWithAndvancedConditioningA(auxRungeKuttaMethod, eps, w, steadyingEquation, opressionCoeffs, invAEigen);
        double[] Y = steadyingProcess.getY();

        setIterationCountWithClarifying(steadyingProcess.getIterationCountWithClarifying());
        setIterationCount(steadyingProcess.getIterationCount());

        double[] solution = steadyingEquation.getSolution(Y);

//        System.out.println("With clarifying: " + getIterationCountWithClarifying());
//
        System.out.println("Y:");
        for (int i = 0; i < Y.length; i++) {
            System.out.println(Y[i]);
        }

        System.out.println("MATRIX x VECTOR count - " + SteadyingProcessWithOpression.getMatrixVectorMultimplicationCount());

        System.out.println("Solution:");
        for (int i = 0; i < solution.length; i++) {
            System.out.println(solution[i]);
        }

        return solution;
    }

    private double getNewW(LinearSteadyingEquation steadEquation, double invA_eighen) {
        double lambda = steadEquation.getJacobiMatrMaxEigenvalue(0.0, null);
        return steadEquation.getTau() * lambda + 1.0 / invA_eighen;
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

