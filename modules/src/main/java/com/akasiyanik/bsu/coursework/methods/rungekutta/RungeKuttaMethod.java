package com.akasiyanik.bsu.coursework.methods.rungekutta;

import org.apache.commons.math3.linear.*;

import java.io.OutputStream;

/**
 *
 *
 * @author Aliaksei Kasiyanik
 */
public class RungeKuttaMethod {

    /** The number of stages */
    private int s;

    /** the Runge–Kutta matrix */
    private double[][] a;

    /** The weights */
    private double[] b;

    /** The nodes */
    private double[] c;

    public RungeKuttaMethod(int dimension, double[][] a,  double[] b, double[] c) {
        s = dimension;

        this.a = a;
        this.b = b;
        this.c = c;
    }

    public void showButcherTableau(){
        for (int i = 0; i < s; i++) {
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.append(b[i]).append("|");
            for (int j = 0; j < s; j++) {
                stringBuilder.append(a[i][j]).append(" ");
            }
            System.out.println(stringBuilder.toString());
        }

        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("_________________").append("\n");
        stringBuilder.append("  ").append("|");
        for (int i = 0; i < s; i++) {
            stringBuilder.append(c[i]).append(" ");
        }
        System.out.println(stringBuilder.toString());
    }

    // Getters and Setters


    public int getS() {
        return s;
    }

    public void setS(int s) {
        this.s = s;
    }

    public double[][] getA() {
        return a;
    }

    public void setA(double[][] a) {
        this.a = a;
    }

    public double[] getB() {
        return b;
    }

    public void setB(double[] b) {
        this.b = b;
    }

    public double[] getC() {
        return c;
    }

    public void setC(double[] c) {
        this.c = c;
    }
}