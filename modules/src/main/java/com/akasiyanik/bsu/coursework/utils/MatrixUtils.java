package com.akasiyanik.bsu.coursework.utils;

import com.sun.corba.se.spi.activation._InitialNameServiceImplBase;
import org.apache.commons.math3.linear.*;

import java.util.Arrays;

/**
 * @author akasiyanik
 */
public final class MatrixUtils {

    private static int matrixVectorCount;

    public static final double ZERO = Math.pow(10, -10);

    public MatrixUtils() {
    }

    public static double[][] kroneckerProduct(double[][] arg1, double[][] arg2) {
        RealMatrix a = new Array2DRowRealMatrix(arg1);
        RealMatrix b = new Array2DRowRealMatrix(arg2);
        return kroneckerProduct(a, b).getData();
    }

    public static double[][] kroneckerProduct(double[] arg1, double[][] arg2) {
        RealMatrix a = new Array2DRowRealMatrix(arg1);
        RealMatrix b = new Array2DRowRealMatrix(arg2);
        return kroneckerProduct(a, b).getData();
    }

    public static double[] kroneckerProduct(double[][] arg1, double[] arg2) {
        RealMatrix a = new Array2DRowRealMatrix(arg1);
        RealMatrix b = new Array2DRowRealMatrix(arg2);
        return kroneckerProduct(a, b).getColumn(0);
    }

    public static double[] kroneckerProduct(double[] arg1, double[] arg2) {
        RealMatrix a = new Array2DRowRealMatrix(arg1);
        RealMatrix b = new Array2DRowRealMatrix(arg2);
        return kroneckerProduct(a, b).getColumn(0);
    }

    private static RealMatrix kroneckerProduct(RealMatrix a, RealMatrix b){
        int resRowDimension = a.getRowDimension() * b.getRowDimension();
        int resColumnDimension = a.getColumnDimension() * b.getColumnDimension();
        RealMatrix resMatrix = new Array2DRowRealMatrix(resRowDimension, resColumnDimension);

        for (int a_i = 0; a_i < a.getRowDimension(); a_i++) {
            for (int a_j = 0; a_j < a.getColumnDimension(); a_j++) {
                double t = a.getEntry(a_i, a_j);
                RealMatrix tempMatrix = b.scalarMultiply(t);
                resMatrix.setSubMatrix(tempMatrix.getData(), a_i * b.getRowDimension(), a_j * b.getColumnDimension());
            }
        }

        return resMatrix;

    }

    public static double[][] getIdentityMatrix(int dimension) {
        double[][] identityMatrix = new double[dimension][dimension];
        for (int i = 0; i < dimension; i++) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;

    }

    public static double[][] transpose (double[] v) {
        RealMatrix m = new Array2DRowRealMatrix(v);
        m = m.transpose();
        return m.getData();
    }

    public static double[] getIdentityVector(int dimension) {
        double[] identityVector = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            identityVector[i] = 1;
        }
        return identityVector;

    }

    public static double[][] scalarMultipy(double scalar, double[][] m) {
        RealMatrix matrix = new Array2DRowRealMatrix(m);
        return matrix.scalarMultiply(scalar).getData();
    }

    public static double[] scalarMultipy(double scalar, double[] v) {
        RealVector vector = new ArrayRealVector(v);
        return vector.mapMultiply(scalar).toArray();
    }

    public static double[][] multiply(double[][] m1, double[][] m2) {
        RealMatrix a = new Array2DRowRealMatrix(m1);
        RealMatrix b = new Array2DRowRealMatrix(m2);

        return a.multiply(b).getData();
    }

    public static double[] multiply(double[][] m, double[] v) {
        matrixVectorCount++;
        RealMatrix a = new Array2DRowRealMatrix(m);
        RealMatrix b = new Array2DRowRealMatrix(v);

        return a.multiply(b).getColumn(0);
    }

    public static double[] multiply(double[] m, double[] v) {
        RealMatrix a = new Array2DRowRealMatrix(m);
        RealMatrix b = new Array2DRowRealMatrix(v);

        return a.multiply(b).getColumn(0);
    }

    public static double[][] subtract(double[][] m1, double[][] m2) {
        RealMatrix a = new Array2DRowRealMatrix(m1);
        RealMatrix b = new Array2DRowRealMatrix(m2);

        return a.subtract(b).getData();
    }

    public static double[] subtract(double[] v1, double[] v2) {
        RealVector a = new ArrayRealVector(v1);
        RealVector b = new ArrayRealVector(v2);

        return a.subtract(b).toArray();
    }

    public static double[][] add(double[][] m1, double[][] m2) {
        RealMatrix a = new Array2DRowRealMatrix(m1);
        RealMatrix b = new Array2DRowRealMatrix(m2);

        return a.add(b).getData();
    }

    public static double[] add(double[] v1, double[] v2) {
        RealVector a = new ArrayRealVector(v1);
        RealVector b = new ArrayRealVector(v2);

        return a.add(b).toArray();
    }

    public static double[] flatten (double[][] matr) {
        int n = matr.length;
        int m = matr[0].length;

        double[] v = new double [n * m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                v[i * m + j] = matr[i][j];
            }
        }
        return v;
    }

    public static double[][] toMatrix(double[] v, int s, int n) {
        double[][] matr = new double [s][n];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < n; j++) {
                matr[i][j] = v[i * n + j];
            }
        }
        return matr;
    }

    public static double norm(double[] v) {
        RealVector vector = new ArrayRealVector(v);
        return vector.getNorm();
    }

    public static double calculateMaxRelativeError(double[] apprSolution, double[] exactSolution) throws Exception {
        if (apprSolution.length != exactSolution.length) {
            throw new Exception("Solution dimensions are not matching. ApprSolution dimension = " + apprSolution.length +
                    "; ExactSolution dimension = " + exactSolution.length + ".");
        }

        int dim = apprSolution.length;

        double[] relError = new double[dim];
        for (int i = 0; i < dim; i++) {
            relError[i] = Math.abs(exactSolution[i] - apprSolution[i]) / Math.abs(exactSolution[i]);
        }

        double max = 0;
        for (int i = 0; i < dim; i++) {
            if (relError[i] > max) {
                max = relError[i];
            }
        }
        return max;
    }


    public static double calculateMaxError(double[] apprSolution, double[] exactSolution) {

        int dim = apprSolution.length;

        double[] relError = new double[dim];
        for (int i = 0; i < dim; i++) {
            relError[i] = Math.abs(exactSolution[i] - apprSolution[i]);
        }

        double max = 0.0;
        for (int i = 0; i < dim; i++) {
            if (relError[i] > max) {
                max = relError[i];
            }
        }
        return max;
    }

    public static double getAverage(double[] v) {
        double t = v[0];
        for (int i = 1; i < v.length; i++) {
            t += v[i];
        }
        return t / v.length;
    }

    public static double[] perComponentDivision(double[] v1, double[] v2) {
        double[] v = new double[v1.length];
        for (int i = 0; i < v1.length; i++) {
            v[i] = v1[i] / v2[i];
        }
        return v;
    }

    public static double[] perComponentDivisionWithZerosInFraction(double[] v1, double[] v2) {
        double[] v = new double[v1.length];
//        System.out.println("----------------------------[");
        for (int i = 0; i < v1.length; i++) {
            if (equalToZero(v2[i])) {
                v[i] = 0.0;
            } else {
                v[i] = v1[i] / v2[i];
            }
//            System.out.println("v: " + v[i] + ", v1: " + v1[i] + ", v2: " + v2[i]) ;
        }
//        System.out.println("]");


        return v;
    }

    public static double getMedian(double[] l) {
        if (l.length == 1) {
            return l[0];
        }
        Arrays.sort(l);
        int middle = l.length / 2;
        if (l.length % 2 == 0) {
            double left = l[middle - 1];
            double right = l[middle];
            return (left + right) / 2;
        } else {
            return l[middle];
        }
    }

    public static double getMedianWithoutZeros(double[] l) {
        int lengthWithoutZeros = 0;
        for (double t : l) {
            if (!equalToZero(t)) {
                lengthWithoutZeros++;
            }
        }
        if (lengthWithoutZeros == 0) {
            return 0.0;
        }

        double[] array = new double[lengthWithoutZeros];
        int i = 0;
        for (double t : l) {
            if (!equalToZero(t)) {
                array[i] = t;
                i++;
            }
        }
        if (array.length == 1) {
            return array[0];
        }
        Arrays.sort(array);
        int middle = array.length / 2;
        if (array.length % 2 == 0) {
            double left = array[middle - 1];
            double right = array[middle];
            return (left + right) / 2;
        } else {
            return array[middle];
        }
    }

    public static boolean equalToZero(double value) {
        return Math.abs(value) <= ZERO;
    }


    public static void printArray(double[] array) {
        System.out.print("----------------[");
        for (double v : array) {
            System.out.print(v + " ");
        }
        System.out.print("]----------------");
    }

    public static double maxComponent(double[] array) {
        double max = -1.0;
        for (double c : array) {
            if (Math.abs(c) > max) {
                max = Math.abs(c);
            }
        }
        return max;
    }

    public static double[][] split(double[] v, int s) {
        double[][] res = new double[s][];
        int partSize = v.length / s;
        for (int i = 0; i < s; i++) {
            double[] tmp = new double[partSize];
            for (int j = 0; j < partSize; j++) {
                tmp[j] = v[i * partSize + j];
            }
            res[i] = tmp;
        }
        return res;
    }

    public static void output(double[] v) {
        System.out.print("[");
        for (double c : v) {
            System.out.print(c + "  ");
        }
        System.out.println("]");
    }

    public static double[] divide(double[] y1, double[] y2) {
        double[] res = new double[y1.length];
        for (int i = 0; i < y1.length; i++) {
            res[i] = y1[i] / y2[i];
        }
        return res;
    }

    public static double[] divide(double[] y1, double l) {
        double[] res = new double[y1.length];
        for (int i = 0; i < y1.length; i++) {
            res[i] = y1[i] / l;
        }
        return res;
    }

    public static double[][][] getBlocksFromDiagonalMatrix(double[][] A, int s) {
        int blockHeight = A.length / s;
        int blockWeidth = A[0].length / s;
        double[][][] blocks = new double[s][][];
        for (int i = 0; i < s; i++) {
            blocks[i] = getBlock(i * blockHeight, i * blockWeidth, blockHeight, blockWeidth, A);
        }
        return blocks;

    }

    private static double[][] getBlock(int ii, int jj, int blockHeight, int blockWeidth, double[][] A) {
        double[][] block = new double[blockHeight][blockWeidth];
        for (int i = ii, i_new = 0; i < ii + blockHeight; i++, i_new++) {
            for (int j = jj, j_new = 0; j < jj + blockWeidth; j++, j_new++) {
                block[i_new][j_new] = A[i][j];
            }
        }
        return block;
    }

    public static double[] abs(double[] a) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            res[i] = Math.abs(a[i]);
        }
        return res;
    }


    public static int getMatrixVectorCount() {
        return matrixVectorCount;
    }

    public static void resetMatrixVectorCount() {
        MatrixUtils.matrixVectorCount = 0;
    }
}
