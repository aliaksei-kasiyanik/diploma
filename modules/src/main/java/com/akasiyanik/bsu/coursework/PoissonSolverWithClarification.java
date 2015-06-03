package com.akasiyanik.bsu.coursework;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcess;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithClarification;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;
import com.akasiyanik.bsu.coursework.problems.PoissonProblem;
import com.akasiyanik.bsu.coursework.utils.FileUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

/**
 * @author: akasiyanik
 */
public class PoissonSolverWithClarification implements Solver {

    private int iterationCount = 0;
    private int iterationCountWithClarifying = 0;
    private static double eps = Math.pow(10, -6);

    public static void main(String[] args) throws FileNotFoundException {
        PoissonSolverWithClarification solver = new PoissonSolverWithClarification();
        double tau = 0.05;
        double t0 = 0.0;

        double[] y0 = new double[20];
        for (int i = 0; i < 20; i++) {
            y0[i] = 1.0 + i;
        }
        solver.solve(t0, y0, tau, eps);
    }

    public double[] solve(double t0, double[] y0, double tau, double eps) {
        RungeKuttaMethod baseRungeKuttaMethod = null;
        AuxRungeKuttaMethod auxRungeKuttaMethod = null;
        try {
            InputStream base = new FileInputStream(new File("E:\\university\\Coursework\\modules\\src\\main\\resources\\RadauIIA-3-Order-Method.txt"));
            InputStream aux = new FileInputStream(new File("E:\\university\\Coursework\\modules\\src\\main\\resources\\9-Order-Generated-AuxMethod.txt"));

            baseRungeKuttaMethod = FileUtils.readBaseRungeKuttaMethod(base);
            auxRungeKuttaMethod = FileUtils.readAuxRungeKuttaMethod(aux);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        SteadyingEquation steadyingEquation = new PoissonProblem(tau, t0, y0, baseRungeKuttaMethod);

        double w = steadyingEquation.getW();
        setIterationCount(0);
        setIterationCountWithClarifying(0);

        SteadyingProcessWithClarification steadyingProcess = new SteadyingProcessWithClarification(auxRungeKuttaMethod, eps, w, steadyingEquation);
        double[] Y = steadyingProcess.getY();

        setIterationCount(steadyingProcess.getIterationCount());
        setIterationCountWithClarifying(steadyingProcess.getIterationCountWithClarifying());

        double[] solution = steadyingEquation.getSolution(Y);

        System.out.println("TOTAL: " + getIterationCount());
        System.out.println("With clarifying: " + getIterationCountWithClarifying());

//        System.out.println("Y:");
//        for (int i = 0; i < Y.length; i++) {
//            System.out.println(Y[i]);
//        }

//        System.out.println("Solution:");
//        for (int i = 0; i < solution.length; i++) {
//            System.out.println(solution[i]);
//        }

        return solution;
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

