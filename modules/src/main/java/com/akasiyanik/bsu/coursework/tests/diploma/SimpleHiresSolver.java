package com.akasiyanik.bsu.coursework.tests.diploma;

import com.akasiyanik.bsu.coursework.HiresSolver;
import com.akasiyanik.bsu.coursework.Solver;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;
import org.apache.commons.lang3.time.StopWatch;

/**
 * @author Aliaksei Kasiyanik
 */
public class SimpleHiresSolver {
    private static double tau = 0.25;
    private static double t0 = 0.0;
    private static double[] y0;
    private static double EPS = Math.pow(10, -4);

    public static void main(String[] args) {
        y0 = getY0();
        solve(new HiresSolver());
    }

    private static double[] getY0() {
        double[] y0 = new double[8];
        y0[0] = 1.0;
        for (int i = 1; i < 7; i++){
            y0[i] = 0.0;
        }
        y0[7] = 0.0057;
        return y0;
    }

    private static void solve(Solver solver) {

        MatrixUtils.resetMatrixVectorCount();
        StopWatch w = new StopWatch();
        w.start();
        solver.solve(t0, y0, tau, EPS);
        w.split();
        System.out.println("Time : " + w.toSplitString());
        w.stop();
        w.reset();

    }

}
