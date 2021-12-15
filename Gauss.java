import javax.crypto.AEADBadTagException;
import java.lang.reflect.Array;
import java.util.Arrays;

public class Gauss {

    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        int n = R.length-1;

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
        //TODO: Diese Methode ist zu implementieren: edge cases?
        if(A.length == 0)
            return null;
        double[] helperB = Arrays.copyOf(b,b.length);

        double[][] helperA = Arrays.copyOf(A,A.length);

        for(int i = 0; i < A.length; i++) {
            swap(helperA, i, findPivot(A,i));       // Reihe mit pivot finden und mit derzeitiger swapen
            for (int j = i+1; j < A[0].length; j++) {
                helperA[i][j] = helperA[i][j] / helperA[i][i];
                helperB[j] -= helperA[i][j] * helperB[i];
                for(int k = 0; k < A.length; k++){
                    helperA[j][k] = helperA[j][k] - helperA[i][k] * helperA[j][i];
                }
            }
        }
        helperB = backSubst(helperA, helperB);
        return helperB;
    }

    private static void swap(double[][] helperA, int i, int pivotRow) {
        if(i == pivotRow)
            return;
        double[] temp = Arrays.copyOf(helperA[i], helperA.length);
        helperA[i] = Arrays.copyOf(helperA[pivotRow], helperA[pivotRow].length);
        helperA[pivotRow] = temp;
    }

    private static int findPivot(double[][] A, int startingRow){
        int currentRow = startingRow;
        double current = 0;
        for(int i = startingRow; i < A.length; i++){
            for(int j = 0; j < A[0].length; j++){
                if(Math.abs(A[i][j]) > current){
                    currentRow = i;
                    current = A[i][j];
                }
            }
        }
        return currentRow;
    }

    /**
     * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
     * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
     *
     * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt):
     * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung
     *  solange durch, bis in einem Schritt alle moeglichen Pivotelemente
     *  numerisch gleich 0 sind (d.h. <1E-10)
     * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
     *  loesen Sie Tx = -v durch Rueckwaertssubstitution
     * -Geben Sie den Vektor (x,1,0,...,0) zurueck
     *
     * Sollte A doch intvertierbar sein, kann immer ein Pivot-Element gefunden werden(>=1E-10).
     * In diesem Fall soll der 0-Vektor zurueckgegeben werden.
     * PARAMETER:
     * A: Eine singulaere Matrix der Groesse n x n
     */
    public static double[] solveSing(double[][] A) {
        //TODO: Diese Methode ist zu implementieren
        if(A.length == 0)
            return null;
        double[] result = new double[A.length];

        double[][] helperA = Arrays.copyOf(A,A.length);
        double[] helperB = new double[A.length];
//        for(int i = 0; i < helperB.length; i++){
//            helperB[i] = 0;
//        }

        double pivot = 1;

        while (pivot > Math.pow(1,-10)) {
            for (int i = 0; i < A.length; i++) {
                int pivotRow = findPivot(A, i);
                swap(helperA, i, pivotRow);       // Reihe mit pivot finden und mit derzeitiger swapen
                for (int j = i + 1; j < A[0].length; j++) {
                    helperA[i][j] = helperA[i][j] / helperA[i][i];
                    helperB[j] -= helperA[i][j] * helperB[i];
                    for (int k = 0; k < A.length; k++) {
                        helperA[j][k] = helperA[j][k] - helperA[i][k] * helperA[j][i];
                    }
                }
                for(int j = i; j < A[0].length; j++){       // find pivot element
                    if(pivot < helperA[i][j])
                        pivot = helperA[i][j];
                }
            }
        }
        for(int i = 0; i < helperB.length; i++){
            helperB[i] = -helperB[i];
        }

        //TODO: Vektor (x,1,0,...,0) zurückgeben
        result = backSubst(helperA, helperB);
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
