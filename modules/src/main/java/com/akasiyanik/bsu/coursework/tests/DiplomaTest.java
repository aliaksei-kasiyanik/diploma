package com.akasiyanik.bsu.coursework.tests;

import com.akasiyanik.bsu.coursework.PoissonSolver;
import com.akasiyanik.bsu.coursework.PoissonSolverWithClarification;
import com.akasiyanik.bsu.coursework.PoissonSolverWithOpression;
import com.akasiyanik.bsu.coursework.PoissonSolverWithOpressionAndClarification;
import org.apache.commons.lang3.time.StopWatch;

/**
 * @author Aliaksei Kasiyanik
 */
public class DiplomaTest {

    private static double eps = Math.pow(10, -8);

    public static void main(String[] args) {
        double tau = 0.05;
        int[] testSizes = {20, 30, 50, 100, 250};

        double t0 = 0.0;

        for (int n : testSizes) {
            double[] y0 = getY0(n);
            System.out.println("-------- SIZE N = " + n + "-----------");

            System.out.println("--- SP ---");

            StopWatch w = new StopWatch();
            w.start();
            PoissonSolver solver = new PoissonSolver();
            solver.solve(t0, y0, tau, eps);
            w.split();
            System.out.println("Time : " + w.toSplitString());
            w.stop();
            w.reset();

            System.out.println("--- SP with CLARIFICATION ---");

            w.start();
            PoissonSolverWithClarification solverWithClarif = new PoissonSolverWithClarification();
            solverWithClarif.solve(t0, y0, tau, eps);
            w.split();
            System.out.println("Time : " + w.toSplitString());
            w.stop();
            w.reset();

            System.out.println("--- SP with PRECONDITIONING ---");

            w.start();
            PoissonSolverWithOpression solverWithOpr = new PoissonSolverWithOpression();
            solverWithOpr.solve(t0, y0, tau, eps);
            w.split();
            System.out.println("Time : " + w.toSplitString());
            w.stop();
            w.reset();

//        PoissonSolverWithOpressionAndClarification superSolver = new PoissonSolverWithOpressionAndClarification();
//        superSolver.solve(t0, y0, tau, eps);
        }

    }

    private static double[] getY0(int n) {
        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0;
        }
        return y0;
    }
}
