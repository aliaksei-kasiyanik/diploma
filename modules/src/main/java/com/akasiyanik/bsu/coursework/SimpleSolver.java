package com.akasiyanik.bsu.coursework;

import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.problems.SimpleProblem;
import com.akasiyanik.bsu.coursework.utils.FileUtils;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcess;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

/**
 * @author: akasiyanik
 */
public class SimpleSolver implements Solver {

    public static void main(String[] args) throws FileNotFoundException {
        SimpleSolver solver = new SimpleSolver();
        double tau = 0.05;
        double t0 = 0.0;
        double[] y0 = new double[] {1.0, 1.0};
        solver.solve(t0, y0, tau, Math.pow(10, -2));
    }

    public double[] solve(double t0, double[] y0, double tau, double eps) {
        RungeKuttaMethod baseRungeKuttaMethod = null;
        AuxRungeKuttaMethod auxRungeKuttaMethod = null;
        try {
            InputStream base = new FileInputStream(new File("C:\\coursework\\modules\\src\\main\\resources\\RadauIIA-3-Order-Method.txt"));
            InputStream aux = new FileInputStream(new File("C:\\coursework\\modules\\src\\main\\resources\\9-Order-Generated-AuxMethod.txt"));

            baseRungeKuttaMethod = FileUtils.readBaseRungeKuttaMethod(base);
            auxRungeKuttaMethod = null;

            auxRungeKuttaMethod = FileUtils.readAuxRungeKuttaMethod(aux);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        SteadyingEquation steadyingEquation = new SimpleProblem(tau, t0, y0, baseRungeKuttaMethod);

//        double eps = 0.01;
        double w = steadyingEquation.getW();
//        System.out.println("w = " + w);

        SteadyingProcess steadyingProcess = new SteadyingProcess(auxRungeKuttaMethod, eps, w, steadyingEquation);
        double[] Y = steadyingProcess.getY();

//        System.out.println("Y:");
//        for (int i = 0; i < Y.length; i++) {
//            System.out.println(Y[i]);
//        }

        double[] solution = steadyingEquation.getSolution(Y);
//        System.out.println("Solution:");
//        for (int i = 0; i < solution.length; i++) {
//            System.out.println(solution[i]);
//        }

        return solution;
    }

    @Override
    public int getIterationCount() {
        return 0;
    }

    @Override
    public int getIterationCountWithClarifying() {
        return 0;
    }
}
