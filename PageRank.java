import java.util.Arrays;
import java.util.Comparator;

public class PageRank {

    /**
     * Diese Methode erstellt die Matrix A~ fuer das PageRank-Verfahren
     * PARAMETER:
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     * zufaellig irgendeine Seite zu besuchen
     */
    public static double[][] buildProbabilityMatrix(int[][] L, double rho) {
        //TODO: Diese Methode ist zu implementieren
        if (L.length == 0)
            return null;
        double[][] probMatr = new double[L.length][L[0].length];
        for (int i = 0; i < L.length; i++) {
            for (int j = 0; j < L[0].length; j++) {
                double temp = 0;
                for (int k = 0; k < L.length; k++)
                    temp += L[k][j];
                if (temp != 0)
                    probMatr[i][j] = 1 / temp;
                else
                    probMatr[i][j] = 0;
                probMatr[i][j] = (1 - rho) * probMatr[i][j] + rho / L.length;
            }

//            for (int j = 0; j < L[0].length; j++) {
//
//            }
        }
        return probMatr;
    }

    /**
     * Diese Methode berechnet die PageRanks der einzelnen Seiten,
     * also das Gleichgewicht der Aufenthaltswahrscheinlichkeiten.
     * (Entspricht dem p-Strich aus der Angabe)
     * Die Ausgabe muss dazu noch normiert sein.
     * PARAMETER:
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, zufaellig irgendeine Seite zu besuchen
     * ,anstatt einem Link zu folgen.
     */
    public static double[] rank(int[][] L, double rho) {
        double[][] probM = buildProbabilityMatrix(L, rho);
        if(probM == null)
            return null;
        for(int i = 0; i < probM.length; i++){      // (A^~ - I) berechnen
            probM[i][i] -= 1;
        }
        double[] p = Gauss.solveSing(probM);
        if(p == null)
            return null;
        double lambda = 0;
        for (double v : p) {
            lambda += v;
        }
        for(int i = 0; i < p.length; i++){
            p[i] /= lambda;
        }
        return p;
    }

    /**
     * Diese Methode erstellt eine Rangliste der uebergebenen URLs nach
     * absteigendem PageRank.
     * PARAMETER:
     * urls: Die URLs der betrachteten Seiten
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     * zufaellig irgendeine Seite zu besuchen
     */
    public static String[] getSortedURLs(String[] urls, int[][] L, double rho) {
        int n = L.length;

        double[] p = rank(L, rho);

        RankPair[] sortedPairs = new RankPair[n];
        for (int i = 0; i < n; i++) {
            sortedPairs[i] = new RankPair(urls[i], p[i]);
        }

        Arrays.sort(sortedPairs, new Comparator<RankPair>() {

            @Override
            public int compare(RankPair o1, RankPair o2) {
                return -Double.compare(o1.pr, o2.pr);
            }
        });

        String[] sortedUrls = new String[n];
        for (int i = 0; i < n; i++) {
            sortedUrls[i] = sortedPairs[i].url;
        }

        return sortedUrls;
    }

    /**
     * Ein RankPair besteht aus einer URL und dem zugehoerigen Rang, und dient
     * als Hilfsklasse zum Sortieren der Urls
     */
    private static class RankPair {
        public String url;
        public double pr;

        public RankPair(String u, double p) {
            url = u;
            pr = p;
        }
    }
}
