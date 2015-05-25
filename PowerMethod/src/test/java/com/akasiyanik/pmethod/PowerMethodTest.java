package com.akasiyanik.pmethod;

import com.akasiyanik.problems.PoissonProblem;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple App.
 */
public class PowerMethodTest extends TestCase {

    private int n = 30;

    private static final double EPS = Math.pow(10, -6);

    public static Test suite() {
        return new TestSuite(PowerMethodTest.class);
    }

    public void testPowerMethod() {
        PoissonProblem problem = new PoissonProblem(n);
        double[][] J = problem.getJ();
        double[] y0 = problem.getY0();

        PowerMethod powerMethod = new PowerMethod(J, y0, EPS);
        System.out.println(powerMethod.solve());
    }

}
