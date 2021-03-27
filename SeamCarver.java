/* *****************************************************************************
 *  Name:
 *  Date:
 *  Description:
 **************************************************************************** */

import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Stack;
import edu.princeton.cs.algs4.StdOut;
import java.util.LinkedList;

public class SeamCarver {
    private static final double MAX = Double.POSITIVE_INFINITY;
    private Picture picture;
    private double[][] energy;
    private double[][] energyTransposed;

    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException();
        this.picture = picture;
        energy = new double[height()][width()];
        energyTransposed = new double[width()][height()]; 
        for (int row = 0; row < this.height(); row++) {
            for (int col = 0; col < this.width(); col++) {
               double en = energy(col, row);
               energy[row][col] = en;
               energyTransposed[col][row] = en;
            }
        }
    }

    public Picture picture() {
        return this.picture;
    }

    public int width() {
        return this.picture.width();
    }

    public int height() {
        return this.picture.height();
    }

    private void validateColumnIndex(int col) {
        if (col < 0 || col >= width()) throw new IllegalArgumentException("Column index out of range.");
    }

    private void validateRowIndex(int row) {
        if (row < 0 || row >= height()) throw new IllegalArgumentException("Row index is out of range.");
    }

    private double calculate(int rgbX, int rgbY) {
        int xr = (rgbX >> 16) & 0xFF;
        int xg = (rgbX >>  8) & 0xFF;
        int xb = rgbX & 0xFF;
        int yr = (rgbY >> 16) & 0xFF;
        int yg = (rgbY >>  8) & 0xFF;
        int yb = rgbY & 0xFF;
        return Math.pow((xr - yr), 2) + Math.pow((xg - yg), 2) + Math.pow((xb - yb), 2);
    }
    public double energy(int x, int y) {
        try {
            validateColumnIndex(x - 1);
            validateColumnIndex(x + 1);
            validateRowIndex(y - 1);
            validateRowIndex(y + 1);
            return Math.sqrt(calculate(picture.getRGB(x-1, y), picture.getRGB(x+1, y)) +
                            calculate(picture.getRGB(x, y-1), picture.getRGB(x, y+1)));
        } catch (IllegalArgumentException e) {
            return 1000.00;
        }
    }

    public int[] findVerticalSeam() {
        double minWeight = MAX;
        Stack<Integer> seam = new Stack<>();
        for (int col = 0; col < width(); col++) {
            double[][] distTo = new double[height()][width()];
            int[][] edgeTo = new int[height()][width()];
            LinkedList<int[]> orderedPoints = new LinkedList<>();
            for (int r = 0; r < height(); r++) {
                for (int c = 0; c < width(); c++) {
                    distTo[r][c] = MAX;
                    edgeTo[r][c] = -1;
                }
                for (int t = col - r - 1; t < col + r + 1; t++) {
                    try {
                        validateColumnIndex(t);
                        orderedPoints.add(new int[]{t, r});
                    } catch (IllegalArgumentException e) {
                       // nothing to do
                    }
                }
            }
            distTo[0][col] = energy[0][col];
            for (int[] p : orderedPoints) {
                LinkedList<int[]> tpp = new LinkedList<>();
                tpp.add(new int[]{p[0]-1, p[1]+1});
                tpp.add(new int[]{p[0], p[1]+1});
                tpp.add(new int[]{p[0]+1, p[1]+1});
                for (int[] tp : tpp) {
                    try {
                        validateColumnIndex(tp[0]);
                        validateRowIndex(tp[1]);
                        if (distTo[tp[1]][tp[0]] > distTo[p[1]][p[0]] + energy[p[1]][p[0]]) 
                        {
                            distTo[tp[1]][tp[0]] = distTo[p[1]][p[0]] + energy[p[1]][p[0]];
                            edgeTo[tp[1]][tp[0]] = p[0];
                        }
                    } catch (IllegalArgumentException e) {
                        // nothing to do
                    }
                }
            }

            for (int i = 0; i < width(); i++) {
                int lastRow = height() - 1;
                if (distTo[lastRow][i] < MAX && distTo[lastRow][i] < minWeight) {
                    minWeight = distTo[lastRow][i];
                    Stack<Integer> path = new Stack<>();
                    for (int e = i; e != -1; e = edgeTo[lastRow--][e])
                        path.push(e);
                    seam = path;
                }
            }
        }
        int[] seamArr = new int[seam.size()];
        int idx = 0;
        while (!seam.isEmpty()) {
           seamArr[idx] = seam.pop();
           idx++;
        }
        return seamArr;
    }

    public int[] findHorizontalSeam() {
        double minWeight = MAX;
        int height = energyTransposed.length;
        int width = energyTransposed[0].length;
        Stack<Integer> seam = new Stack<>();
        for (int col = 0; col < width; col++) {
            double[][] distTo = new double[height][width];
            int[][] edgeTo = new int[height][width];
            LinkedList<int[]> orderedPoints = new LinkedList<>();
            for (int r = 0; r < height; r++) {
                for (int c = 0; c < width; c++) {
                    distTo[r][c] = MAX;
                    edgeTo[r][c] = -1;
                }
                for (int t = col - r - 1; t < col + r + 1; t++) 
                    if (t >= 0 && t < width) orderedPoints.add(new int[]{t, r});
            }
            distTo[0][col] = energyTransposed[0][col];
            for (int[] p : orderedPoints) {
                LinkedList<int[]> tpp = new LinkedList<>();
                tpp.add(new int[]{p[0]-1, p[1]+1});
                tpp.add(new int[]{p[0], p[1]+1});
                tpp.add(new int[]{p[0]+1, p[1]+1});
                for (int[] tp : tpp) {
                    if (tp[0] >= 0 && tp[0] < width && tp[1] >= 0 && tp[1] < height) {
                        if (distTo[tp[1]][tp[0]] > distTo[p[1]][p[0]] + energyTransposed[p[1]][p[0]]) 
                        {
                            distTo[tp[1]][tp[0]] = distTo[p[1]][p[0]] + energyTransposed[p[1]][p[0]];
                            edgeTo[tp[1]][tp[0]] = p[0];
                        }
                    }
                }
            }

            for (int i = 0; i < width; i++) {
                int lastRow = height - 1;
                if (distTo[lastRow][i] < MAX && distTo[lastRow][i] < minWeight) {
                    minWeight = distTo[lastRow][i];
                    Stack<Integer> path = new Stack<>();
                    for (int e = i; e != -1; e = edgeTo[lastRow--][e])
                        path.push(e);
                    seam = path;
                }
            }
        }
        int[] seamArr = new int[seam.size()];
        int idx = 0;
        while (!seam.isEmpty()) {
           seamArr[idx] = seam.pop();
           idx++;
        }
        return seamArr;
    }

    public static void main(String[] args) {
        Picture pic = new Picture(args[0]);
        SeamCarver sc = new SeamCarver(pic);
        StdOut.print("Vertical Seam:\n");
        int[] s = sc.findVerticalSeam();
        for (int i : s) StdOut.printf("%4d", i);
        StdOut.println();
        StdOut.print("Horizontal Seam:\n");
        int[] h = sc.findHorizontalSeam();
        for (int i : h) StdOut.printf("%4d", i);
        StdOut.println();
    }
}
