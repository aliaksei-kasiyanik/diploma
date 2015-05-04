package com.akasiyanik.bsu.coursework.equations;

import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.multiply;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.subtract;

/**
 * @author Aliaksei Kasiyanik
 */
public abstract class LinearSteadyingEquation extends SteadyingEquation {

    public LinearSteadyingEquation(RungeKuttaMethod implMethod) {
        super(implMethod);
    }

    public LinearSteadyingEquation(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        super(tau, t0, y0, implMethod);
    }

    public abstract double[][] getJMatrix();

    public double[] r(double[] Y) {   //dim Y: s*n x 1
        double[] sum1 = multiply(getJMatrix(), getY0());
        double[] sum2 = new double[n * s];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < n; j++) {
                sum2[i * n + j] = sum1[j];
            }
        }
        return add(multiply(getG(), Y), sum2);
    }

    public double[][] getG() {
        double[][] A = implicitRungeKuttaMethod.getA(); //dim: s x s
        double[][] J = getJMatrix(); //dim: n x n

        double[][] res = kroneckerProduct(A, J);  //A(x)J
        res = scalarMultipy(tau, res);  // tau A(x)J
        double[][] I = getIdentityMatrix(res.length);
        return subtract(res, I);// tau A(x)J - I
    }


}
