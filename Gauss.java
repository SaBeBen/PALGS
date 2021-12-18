import java.util.Arrays;
import java.util.Stack;

public class Gauss {

    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        int n = R.length - 1;

        if (R.length == 0) {
            return null;
        }
        double[] result = new double[R.length];
        for (int i = n; i >= 0; i--) {
            for (int j = i + 1; j <= n; j++) {
                b[i] -= R[i][j] * result[j];
            }
            result[i] = b[i] / R[i][i];
        }
        return result;
    }

    /**
     * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
     * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden.
     * PARAMETER: A:
     * Eine regulaere Matrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] solve(double[][] A, double[] b) {
        if (A.length == 0)
            return null;

        double[][] helperA = new double[A.length][A.length + 1];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                helperA[i][j] = A[i][j];
            }
        }
        for (int i = 0; i < helperA.length; i++) {
            helperA[i][helperA[0].length - 1] = b[i];
        }

        for (int i = 0; i < helperA.length; i++) {
            if (helperA[i][i] == 0) {       // Falls bereits ein Diagonalfeld 0 ist...
                for (int j = 0; j < helperA.length; j++) {
                    if (helperA[j][i] != 0) {
                        swap(helperA, i, j); // muss dieses natürlich != 0 sein
                    }
                }
            }
            swap(helperA, i, findPivot(helperA, i, 0));       // Reihe mit pivot finden und mit derzeitiger swapen
            for (int j = i + 1; j < helperA.length; j++) {
                double factor = (helperA[j][i] / helperA[i][i]);    // Faktor bestimmen, sodass ...
                for (int k = i; k < helperA[0].length; k++) {
                    helperA[j][k] = helperA[j][k] - factor * helperA[i][k];     // ...helper[j][i] == 0
                }
            }
        }
        double[] helperB = Arrays.copyOf(b, b.length);

        for (int i = 0; i < helperA.length; i++) {
            helperB[i] = helperA[i][helperA[0].length - 1];
        }

        double[][] elimMatrix = new double[A.length][A[0].length];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                elimMatrix[i][j] = helperA[i][j];
            }
        }
        helperB = backSubst(elimMatrix, helperB);
        return helperB;
    }

    private static void swap(double[][] helperA, int i, int pivotRow) {
        if (i == pivotRow)
            return;
        double[] temp = Arrays.copyOf(helperA[i], helperA[i].length);
        helperA[i] = Arrays.copyOf(helperA[pivotRow], helperA[pivotRow].length);
        helperA[pivotRow] = temp;
    }

//    private static void swap(double[][] helperA, int i, int pivotRow, ) {
//        if (i == pivotRow)
//            return;
//        double[] temp = Arrays.copyOf(helperA[i], helperA[i].length);
//        helperA[i] = Arrays.copyOf(helperA[pivotRow], helperA[pivotRow].length);
//        helperA[pivotRow] = temp;
//    }

    // wenn die Matrix augmentiert ist, d h Ax = b -> Ab, dann muss die letzte Spalte nicht untersucht werden
    private static int findPivot(double[][] helperA, int startingRow, int isAugmented) {
        int currentRow = startingRow;
        double current = 0;
        for (int i = startingRow; i < helperA.length; i++) {
            for (int j = 0; j < helperA[0].length - 1 + isAugmented; j++) {
                if (Math.abs(helperA[i][j]) > current) {
                    currentRow = i;
                    current = Math.abs(helperA[i][j]);
                }
            }
        }
        return currentRow;
    }

    /**
     * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
     * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
     * <p>
     * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt):
     * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung
     * solange durch, bis in einem Schritt alle moeglichen Pivotelemente
     * numerisch gleich 0 sind (d.h. <1E-10)
     * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
     * loesen Sie Tx = -v durch Rueckwaertssubstitution
     * -Geben Sie den Vektor (x,1,0,...,0) zurueck
     * <p>
     * Sollte A doch intvertierbar sein, kann immer ein Pivot-Element gefunden werden(>=1E-10).
     * In diesem Fall soll der 0-Vektor zurueckgegeben werden.
     * PARAMETER:
     * A: Eine singulaere Matrix der Groesse n x n
     */
    public static double[] solveSing(double[][] A) {
        // Following case does not seem to work
        //  |1.74 3.01 -17.0 |			|0.0|
        //  |3.12 9.0 1.11 |	x	=	|0.0|
        //  |4.86 12.01 -15.89 |		|0.0
        if (A.length == 0)
            return null;

        double[][] helperA = Arrays.copyOf(A, A.length);

        double pivot = 1;
        int pivotRow = 0;
        int lastRow = 0;

        while (pivot > Math.pow(1, -10)) {
            for (int i = 0; i < helperA.length; i++) {
                pivotRow = findPivot(helperA, i, 1);
                swap(helperA, i, pivotRow);       // Reihe mit pivot finden und mit derzeitiger swapen

                for (int j = i + 1; j < helperA.length; j++) {
                    double factor = (helperA[j][i] / helperA[i][i]);    // Faktor bestimmen, sodass ...
                    for (int k = i; k < helperA[0].length; k++) {
                        helperA[j][k] = helperA[j][k] - factor * helperA[i][k];     // ...helper[j][i] == 0
                    }
                }

                for (int j = findPivot(helperA, i, 1); j < A[0].length; j++) {       // find pivot element
                    if (pivot > helperA[i][j])
                        pivot = helperA[i][j];
                }
                lastRow = i;
            }
        }

        // im Folgenden werden v und T ermittelt, um die Rückwärtssubst durchzuführen
        double[][] T = new double[lastRow + 1][lastRow + 1];
        double[] v = new double[lastRow + 1];
        for (int i = 0; i <= lastRow; i++) {
            for (int j = 0; j <= lastRow; j++) {
                T[i][j] = helperA[i][j];
            }
        }
        for (int i = 0; i <= lastRow; i++) {
            v[i] = -helperA[i][lastRow + 1];
        }
        double[] result = new double[A.length];
        v = backSubst(T, v);
        if (v == null)
            throw new IllegalStateException("backSubst war nicht erfolgreich!");
        for (int i = 0; i < v.length; i++) {
            result[i] = v[i];
        }
        //anfügen der 1 und der 0en
        if (result.length == v.length)
            return result;
        result[v.length] = 1;
        for (int i = v.length + 1; i < result.length; i++) {
            result[i] = 0;
        }
        return result;
    }

    /**
     * Diese Methode berechnet das Matrix-Vektor-Produkt A*x mit A einer nxm
     * Matrix und x einem Vektor der Laenge m. Sie eignet sich zum Testen der
     * Gauss-Loesung
     */
    public static double[] matrixVectorMult(double[][] A, double[] x) {
        int n = A.length;
        int m = x.length;

        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            y[i] = 0;
            for (int j = 0; j < m; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }
}
