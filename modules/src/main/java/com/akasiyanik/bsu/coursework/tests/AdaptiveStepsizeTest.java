package com.akasiyanik.bsu.coursework.tests;

import com.akasiyanik.bsu.coursework.HiresSolver;
import com.akasiyanik.bsu.coursework.SimpleSolver;
import com.akasiyanik.bsu.coursework.Solver;
import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;

/**
 * @author: akasiyanik
 */
public class AdaptiveStepsizeTest {

    public static void main(String[] args) {
        AdaptiveStepsizeTest adaptiveStepsizeTest = new AdaptiveStepsizeTest();
        double[] apprSolution = adaptiveStepsizeTest.test();
        double[] exactSolution = {
                0.00073713125733260241156,
                0.00014424857263162552644,
                5.8887297409677703125e-05,
                0.0011756513432832147847,
                0.0023863561988375892563,
                0.0062389682529970192856,
                0.0028499983948810357105,
                0.0028500016051188838302,
        };

        double[] relError = new double[8];
        for (int i = 0; i < 8; i++) {
            relError[i] = Math.abs(exactSolution[i] - apprSolution[i]) / Math.abs(exactSolution[i]);
        }

        double max = 0.0;
        try {
            max = MatrixUtils.calculateMaxRelativeError(exactSolution, apprSolution);
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("RELATIVE ERROR = " + max);

    }

    public double[] test() {
        HiresSolver solver = new HiresSolver();

//        double tau = 0.05;
//        double t0 = 0.0;

        //-------- The initial vector --------
        double[] y0 = new double[8];
        y0[0] = 1.0;
        for (int i = 1; i < 7; i++){
            y0[i] = 0.0;
        }
        y0[7] = 0.0057;

        //-------- The end point of the integration interval --------
        double H = 321.8122;

        double tol = Math.pow(10, -10);
        int p = 3;
        double alpha = 0.8; //[0.7, 0.9]

        double h0 = 0.00001;
        double x0 = 0.0;

        double x = x0;
        double[] y = y0;
        double h = h0;
        double X = x0 + H;
        double delta;
        double[] y2;
        double h_new;
        double[] y2_with_wave = y0;

        int j = 0;

        int total_iteration_count = 0;
        int total_iteration_with_clarifying_count = 0;

        while (x != X) {
            do {
                y2 = solver.solve(x + h, solver.solve(x, y, h, tol), h, tol);
                total_iteration_count += solver.getIterationCount();
                total_iteration_with_clarifying_count += solver.getIterationCountWithClarifying();

                y2_with_wave = solver.solve(x, y, 2 * h, tol);
                total_iteration_count += solver.getIterationCount();
                total_iteration_with_clarifying_count += solver.getIterationCountWithClarifying();

                double err = norm(subtract(y2, y2_with_wave)) / (Math.pow(2.0, 3.0) - 1.0);
                delta = Math.pow(tol / Math.abs(err), 1.0 / (p + 1.0));
                System.out.println("err = " + err + "\tdelta = " + delta);
                h_new = alpha * delta * h;
                j++;
                if (delta < 1) {
                    h = h_new;
                    System.out.println("i =" + j + "Ejected!\n--------------");
                }
            } while (delta < 1);

            System.out.println("i =" + j + "\nAccepted.\n x + 2h: " + (x + 2 * h) + "\ny2:");
            for (int i = 0; i < y2.length; i++) {
                System.out.println(y2[i]);
            }
            System.out.println("--------------");
            x = x + 2 * h;
            y = y2;
            h = Math.min(h_new, (X - x) / 2.0);
        }
        System.out.println("y2_with_wave:");
        for (int i = 0; i < y2_with_wave.length; i++) {
            System.out.println(y2_with_wave[i]);
        }

        System.out.println("total_iteration_count: " + total_iteration_count);
        System.out.println("total_iteration_with_clarifying_count: " + total_iteration_with_clarifying_count);


        return y2_with_wave;
    }
}
