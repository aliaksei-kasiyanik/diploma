package com.akasiyanik.bsu.coursework.solvers;

import com.akasiyanik.bsu.coursework.Solver;
import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcess;
import com.akasiyanik.bsu.coursework.methods.SteadyingProcessWithClarification;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;
import com.akasiyanik.bsu.coursework.problems.HiresProblem;
import com.akasiyanik.bsu.coursework.utils.FileUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

/**
 * @author Aliaksei Kasiyanik
 */
public class HiresSolverWithClarification implements Solver {

    private int iterationCount = 0;
    private int iterationCountWithClarifying = 0;

    public static void main(String[] args) throws FileNotFoundException {
        HiresSolverWithClarification solver = new HiresSolverWithClarification();
        double tau = 0.5;
        double t0 = 0.0;

        double[] y0 = new double[8];
        y0[0] = 1.0;
        for (int i = 1; i < 7; i++){
            y0[i] = 0.0;
        }
        y0[7] = 0.0057;
        solver.solve(t0, y0, tau, Math.pow(10, -10));
    }

    public double[] solve(double t0, double[] y0, double tau, double eps) {
        RungeKuttaMethod baseRungeKuttaMethod = null;
        AuxRungeKuttaMethod auxRungeKuttaMethod = null;
        try {
            InputStream base = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\RadauIIA-5-Order-Method.txt"));
            InputStream aux = new FileInputStream(new File("E:\\university\\diploma\\modules\\src\\main\\resources\\9-Order-Generated-AuxMethod.txt"));

            baseRungeKuttaMethod = FileUtils.readBaseRungeKuttaMethod(base);
            auxRungeKuttaMethod = FileUtils.readAuxRungeKuttaMethod(aux);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        SteadyingEquation steadyingEquation = new HiresProblem(tau, t0, y0, baseRungeKuttaMethod);

//        double eps = Math.pow(10, -5);
        double w = steadyingEquation.getW();
        setIterationCount(0);
//        setIterationCountWithClarifying(0);
//        System.out.println("w = " + w);

        SteadyingProcess steadyingProcess = new SteadyingProcessWithClarification(auxRungeKuttaMethod, eps, w, steadyingEquation);
        double[] Y = steadyingProcess.getY();

        setIterationCount(steadyingProcess.getIterationCount());
        setIterationCountWithClarifying(steadyingProcess.getIterationCountWithClarifying());
        System.out.println("ITER COUNT" + getIterationCount());

        double[] solution = steadyingEquation.getSolution(Y);

//        System.out.println("Y:");
//        for (int i = 0; i < Y.length; i++) {
//            System.out.println(Y[i]);
//        }
//
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
