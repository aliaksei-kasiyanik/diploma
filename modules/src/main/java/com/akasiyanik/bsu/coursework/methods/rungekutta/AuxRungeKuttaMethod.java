package com.akasiyanik.bsu.coursework.methods.rungekutta;

/**
 * @author: akasiyanik
 */
public class AuxRungeKuttaMethod {
    private int sigma; // количество стадий
    private int n; // количество уравнений
    private double[] alpha;
    private double[][] betta;
    private double gamma;

    private double derivativeRin0; // производная многочлена перехода в нуле

    public double getGamma() {
        return gamma;
    }

    public void setGamma(double gamma) {
        this.gamma = gamma;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public double[][] getBetta() {
        return betta;
    }

    public void setBetta(double[][] betta) {
        this.betta = betta;
    }

    public int getSigma() {
        return sigma;
    }

    public void setSigma(int sigma) {
        this.sigma = sigma;
    }

    public double[] getAlpha() {
        return alpha;
    }

    public void setAlpha(double[] alpha) {
        this.alpha = alpha;
    }

    public double getDerivativeRin0() {
        return derivativeRin0;
    }

    public void setDerivativeRin0(double derivativeRin0) {
        this.derivativeRin0 = derivativeRin0;
    }
}