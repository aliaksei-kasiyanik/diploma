package com.akasiyanik.bsu.coursework.equations;

import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.*;
import static com.akasiyanik.bsu.coursework.utils.MatrixUtils.subtract;

import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

/**
 * @author akasiyanik
 */
public abstract class SteadyingEquation {

    protected RungeKuttaMethod implicitRungeKuttaMethod;
    protected double tau;
    protected double t0;
    protected double[] y0;
    protected int n;
    protected int s;

    public abstract double[] f(double t, double[] arg);

    public abstract double[][] JacobiMatrix (double[] arg);

    public abstract double getJacobiMatrMaxEigenvalue(double t, double[] arg);

    public abstract double getAMatrMaxEigenvalue();

    public abstract double[][] getG();

    public RungeKuttaMethod getImplicitRungeKuttaMethod() {
        return implicitRungeKuttaMethod;
    }

    public double getW() {
        double nu = tau * Math.abs(getJacobiMatrMaxEigenvalue(t0, y0) * getAMatrMaxEigenvalue() - 1.0);
        return 1.0 / nu * 0.7;
    }

    public double[] getY0() {
        return y0;
    }

    public SteadyingEquation(RungeKuttaMethod implMethod) {
        tau = 0.05;
        t0 = 0.0;
        y0 = new double[]{1.0, 1.0};
        implicitRungeKuttaMethod = implMethod;
        s = implicitRungeKuttaMethod.getS();
        n = y0.length;
    }

    public SteadyingEquation(double tau, double t0, double[] y0, RungeKuttaMethod implMethod) {
        this.tau = tau;
        this.t0 = t0;
        this.y0 = y0;
        n = y0.length;
        s = implMethod.getS();
        implicitRungeKuttaMethod = implMethod;
    }

    public double[] r(double[] Y) {   //dim Y: s*n x 1

        double[][] A = implicitRungeKuttaMethod.getA(); //dim: s x s
        double[] e = getIdentityVector(s); //dim: s x 1
        double[][] I = getIdentityMatrix(n); //dim: n x n

        double[] F = funct_F(t0, Y);  // s*n x 1

        double[][] sum1 = kroneckerProduct(A, I);  //A(x)I
        sum1 = scalarMultipy(tau, sum1);// tau*(A(x)I)

        double[] sum2 = multiply(sum1, F); //tau*(A(x)I)*F
        sum2 = subtract(sum2, Y); //tau*(A(x)I)*F - Y

        double[] sum3 = kroneckerProduct(e, y0);//e(x)y0

        sum2 = add(sum2, sum3); //tau*(A(x)I)*F - Y + e(x)y0

        return sum2;
    }


    public double[] funct_F(double t, double[] Y) {
        double[] c = implicitRungeKuttaMethod.getC();
        double[][] tmp = new double[s][n];

        for (int i = 0; i < s; i++) {
            double[] Yi = new double[n];
            for (int j = 0; j < n; j++) {
                Yi[j] = Y[i * n + j];
            }
            tmp[i] = f(t + c[i] * tau, Yi);
        }
        return flatten(tmp);
    }

    public double[] getSolution(double[] Y) {
        double[] b = implicitRungeKuttaMethod.getB(); // s x 1
        double[][] I = getIdentityMatrix(n); //dim: n x n

        double[] F = funct_F(t0, Y);  // s*n x 1

        double[] res = multiply(kroneckerProduct(transpose(b), I), F);
        res = scalarMultipy(tau, res);// tau*(A(x)I)

        res = add(y0, res); //tau*(A(x)I)*F - Y + e(x)y0

        return res;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public int getS() {
        return s;
    }

    public void setS(int s) {
        this.s = s;
    }

    public double getT0() {
        return t0;
    }

    public void setT0(double t0) {
        this.t0 = t0;
    }

    public double getTau() {
        return tau;
    }

    public void setTau(double tau) {
        this.tau = tau;
    }
}
