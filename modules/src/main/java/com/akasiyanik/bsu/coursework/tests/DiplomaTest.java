package com.akasiyanik.bsu.coursework.tests;

import com.akasiyanik.bsu.coursework.PoissonSolver;
import com.akasiyanik.bsu.coursework.PoissonSolverWithClarification;
import com.akasiyanik.bsu.coursework.PoissonSolverWithOpression;
import com.akasiyanik.bsu.coursework.PoissonSolverWithOpressionAndClarification;

/**
 * @author Aliaksei Kasiyanik
 */
public class DiplomaTest {

    private static double eps = Math.pow(10, -6);

    public static void main(String[] args) {
        double tau = 0.05;
        int n = 20;

        double t0 = 0.0;

        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0;
        }


        System.out.println("--- SP ---");

        PoissonSolver solver = new PoissonSolver();
        solver.solve(t0, y0, tau, eps);


        System.out.println("--- SP with CLARIFICATION ---");

        PoissonSolverWithClarification solverWithClarif = new PoissonSolverWithClarification();
        solverWithClarif.solve(t0, y0, tau, eps);


        System.out.println("--- SP with PRECONDITIONING ---");

        PoissonSolverWithOpression solverWithOpr = new PoissonSolverWithOpression();
        solverWithOpr.solve(t0, y0, tau, eps);

//        PoissonSolverWithOpressionAndClarification superSolver = new PoissonSolverWithOpressionAndClarification();
//        superSolver.solve(t0, y0, tau, eps);

    }
}
