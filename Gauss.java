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

        // Derzeitige Reihe ermitteln
        int currentRow = 0;
//        double current = 0;
//        int currentRow = 0;
//        for(int i = 0; i < A.length; i++){
//            for(int j = 0; j < A[0].length; j++){
//                if(Math.abs(A[i][j]) > current){
//                    currentRow = i;
//                    current = A[i][j];
//                }
//            }
//        }

        for(int i = 0; i < A.length; i++) {
            swap(helperA, i, findPivot(A,i));
            for (int j = i+1; j < A[0].length; j++) {
                helperA[i][j] = helperA[i][j] / helperA[i][i];
                for(int k = 0; k < A.length; k++){
                    helperA[j][k] = helperA[j][k] - helperA[i][k] * helperA[j][i];
                }
            }
        }
        backSubst(helperA, helperB);
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
        return null;
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
