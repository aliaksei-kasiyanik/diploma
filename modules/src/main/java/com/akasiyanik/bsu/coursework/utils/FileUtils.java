package com.akasiyanik.bsu.coursework.utils;

import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;
import com.akasiyanik.bsu.coursework.methods.rungekutta.RungeKuttaMethod;

import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Locale;
import java.util.Scanner;

/**
 * @author: akasiyanik
 */
public final class FileUtils {

    public FileUtils() {
    }

    public static RungeKuttaMethod readBaseRungeKuttaMethod(InputStream inputStream) throws FileNotFoundException {

        Scanner scanner = new Scanner(inputStream);

        scanner.useLocale(Locale.US);
        int s = scanner.nextInt();
//        scanner.nextLine();

        double[][] a = new double[s][s];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < s; j++) {
                a[i][j] = scanner.nextDouble();
            }
            scanner.nextLine();
        }

        scanner.nextLine();
        double[] b = new double[s];
        for (int i = 0; i < s; i++) {
            b[i] = scanner.nextDouble();
        }

        scanner.nextLine();
        double[] c = new double[s];
        for (int i = 0; i < s; i++) {
            c[i] = scanner.nextDouble();
        }

        return new RungeKuttaMethod(s, a, b, c);
    }

    public static AuxRungeKuttaMethod readAuxRungeKuttaMethod(InputStream inputStream) throws FileNotFoundException {
        Scanner scanner = new Scanner(inputStream);
        scanner.useLocale(Locale.US);
        int sigma = scanner.nextInt();
        scanner.nextLine();

        int n;
        double[] alpha;
        double[][] betta;
        AuxRungeKuttaMethod auxRungeKuttaMethod = new AuxRungeKuttaMethod();

        if (sigma % 2 == 0) {
            n = sigma / 2;
            alpha = new double[n];
            betta = new double[n][2];
            for (int i = 0; i < n; i++) {
                alpha[i] = scanner.nextDouble();
                betta[i][0] = scanner.nextDouble();
                betta[i][1] = scanner.nextDouble();
                scanner.nextLine();
            }
        }
        else {
            n = sigma / 2 + 1;
            alpha = new double[n - 1];
            betta = new double[n - 1][2];
            for (int i = 0; i < n - 1; i++) {
                alpha[i] = scanner.nextDouble();
                betta[i][0] = scanner.nextDouble();
                betta[i][1] = scanner.nextDouble();
                scanner.nextLine();
            }
            double gamma = scanner.nextDouble();
            auxRungeKuttaMethod.setGamma(gamma);
        }
        scanner.nextLine();

        double derRin0 = scanner.nextDouble();

        auxRungeKuttaMethod.setSigma(sigma);
        auxRungeKuttaMethod.setN(n);
        auxRungeKuttaMethod.setAlpha(alpha);
        auxRungeKuttaMethod.setBetta(betta);
        auxRungeKuttaMethod.setDerivativeRin0(derRin0);

        return auxRungeKuttaMethod;
    }

    public static RungeKuttaMethod readExplicitRungeKuttaInitalData(InputStream inputStream) throws FileNotFoundException {

        Scanner scanner = new Scanner(inputStream);
        int s = scanner.nextInt();
        scanner.nextLine();

        double[][] a = new double[s][s];
        for (int i = 1; i < s; i++) {
            for (int j = 0; j < i; j++) {
                a[i][j] = scanner.nextDouble();
            }
            scanner.nextLine();
        }

//        scanner.nextLine();
        double[] b = new double[s];
        for (int i = 0; i < s; i++) {
            b[i] = scanner.nextDouble();
        }

        double[] c = new double[s];

        return new RungeKuttaMethod(s, a, b, c);
    }

}
