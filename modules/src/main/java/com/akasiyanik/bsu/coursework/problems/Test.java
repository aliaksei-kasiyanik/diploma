package com.akasiyanik.bsu.coursework.problems;

import com.akasiyanik.bsu.coursework.utils.MatrixUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Aliaksei Kasiyanik
 */
public class Test {

    public static void main(String[] args) {
        double t = 1;
        double[] arg = new double[30];
//        for (int i = 0; i < 30; i++) {
//            arg[i] = 1.0;
//            if (i == 0) arg[i] = 0;
//        }
//
//        double f[] = f(t, arg);
//        MatrixUtils.printArray(f);

        System.out.println(getJacobiMatrMaxEigenvalue());

    }

    public static double getJacobiMatrMaxEigenvalue() {
        int n = 30;
        List<Double> evs = new ArrayList<Double>();
        for (int i = 0; i < n; i++) {
            double e = -4 * Math.pow(n - 1.0, 2.0) * Math.pow(Math.sin(Math.PI * i / (2 * n + 2)), 2);
            evs.add(Math.abs(e));
        }
        return Collections.max(evs);

//             return -17.1078;
    }

    public static double[] f(double t, double[] arg) {
        int n = 30;
        double[] res = new double[30];
        for (int i = 0; i < n; i++) {
            res[i] += arg[i] * -1922.0;
            if (i < n - 1) {
                res[i] += arg[i + 1] * 961;
            }
            if (i > 0) {
                res[i] += arg[i - 1] * 961;
            }
        }
        return res;


//        double f0 = -1.71 * arg[0] + 0.43 * arg[1] + 8.32 * arg[2] + 0.0007;
//        double f1 = 1.71 * arg[0] - 8.75 * arg[1];
//        double f2 = -10.03 * arg[2] + 0.43 * arg[3] + 0.035 * arg[4];
//        double f3 = 8.32 * arg[1] + 1.71 * arg[2] - 1.12 * arg[3];
//        double f4 = -1.745 * arg[4] + 0.43 * arg[5] + 0.43 * arg[6];
//        double f5 = -280 * arg[5] * arg[7] + 0.69 * arg[3] + 1.71 * arg[4] - 0.43 * arg[5] + 0.69 * arg[6];
//        double f6 = 280 * arg[5] * arg[7] - 1.81 * arg[6];
//        double f7 = -280 * arg[5] * arg[7] + 1.81 * arg[6];
//
//        double[] res = {f0, f1, f2, f3, f4, f5, f6, f7};
//        return res;
    }
}
