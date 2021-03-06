package com.akasiyanik.bsu.coursework.tests.diploma;

import com.akasiyanik.bsu.coursework.PoissonSolver;
import com.akasiyanik.bsu.coursework.PoissonSolverWithAdvancedOpression;
import com.akasiyanik.bsu.coursework.PoissonSolverWithClarification;
import com.akasiyanik.bsu.coursework.Solver;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;
import org.apache.commons.lang3.time.StopWatch;

/**
 * @author Aliaksei Kasiyanik
 */
public class DiplomaTest3 {

    private static double eps = Math.pow(10, -8);
    private static double tau = 0.05;
    private static double t0 = 0.0;
    private static double[] y0;

    public static void main(String[] args) {

//        int[] testSizes = {300};
        int[] testSizes = {20, 30, 50, 100, 150, 200, 300};

        for (int n : testSizes) {
            y0 = getY0(n);
            System.out.println("-------- SIZE N = " + n + "-----------");

            System.out.println("--- SP ---");
            solve(new PoissonSolver());
////
            System.out.println("--- SP with CLARIFICATION ---");
            solve(new PoissonSolverWithClarification());
////
            System.out.println("--- SP with PRECONDITIONING ---");
            solve(new PoissonSolverWithAdvancedOpression());

//        PoissonSolverWithOpressionAndClarification superSolver = new PoissonSolverWithOpressionAndClarification();
//        superSolver.solve(t0, y0, tau, eps);
        }

    }

    private static void solve(Solver solver) {
        MatrixUtils.resetMatrixVectorCount();
        StopWatch w = new StopWatch();
        w.start();
        solver.solve(t0, y0, tau, eps);
        System.out.println("Iter COUNT: " + solver.getIterationCount());
        System.out.println("MATRIX x VECTOR count = " + MatrixUtils.getMatrixVectorCount());
        w.split();
        System.out.println("Time : " + w.toSplitString());
        w.stop();
        w.reset();

    }

    private static double[] getY0(int n) {
        double[] y0 = new double[n];
        for (int i = 0; i < n; i++) {
            y0[i] = 1.0 + i;
        }
        return y0;
    }
}
