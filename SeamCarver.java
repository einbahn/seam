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

    public SeamCarver(Picture picture) {
        if (picture == null)
            throw new IllegalArgumentException();
        this.picture = picture;
    }

    public Picture picture() {
        return picture;
    }

    public int width() {
        return picture.width();
    }

    public int height() {
        return picture.height();
    }

    private double calculate(int rgbX, int rgbY) {
        int xr = (rgbX >> 16) & 0xFF;
        int xg = (rgbX >> 8) & 0xFF;
        int xb = rgbX & 0xFF;
        int yr = (rgbY >> 16) & 0xFF;
        int yg = (rgbY >> 8) & 0xFF;
        int yb = rgbY & 0xFF;
        return Math.pow((xr - yr), 2) + Math.pow((xg - yg), 2) + Math.pow((xb - yb), 2);
    }

    public double energy(int x, int y) {
        if (x-1 >= 0 && x + 1 < width() && y - 1 >= 0 && y + 1 < height()) {
            return Math.sqrt(
                calculate(picture.getRGB(x-1, y), picture.getRGB(x+1, y)) + 
                calculate(picture.getRGB(x, y-1), picture.getRGB(x, y+1))
            );
        } else {
            return 1000;
        }
    }

    public int[] findVerticalSeam() {
        double minWeight = MAX;
        int height = height();
        int width = width();
        Stack<Integer> seam = new Stack<>();
        double[][] energy = new double[height][width];
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
                    if (t >= 0 && t < width) orderedPoints.add(new int[] { t, r });
            }
            if (energy[0][col] == 0) energy[0][col] = energy(col, 0);
            distTo[0][col] = energy[0][col];
            for (int[] p : orderedPoints) {
                LinkedList<int[]> tpp = new LinkedList<>();
                tpp.add(new int[] { p[0] - 1, p[1] + 1 });
                tpp.add(new int[] { p[0], p[1] + 1 });
                tpp.add(new int[] { p[0] + 1, p[1] + 1 });
                for (int[] tp : tpp) {
                    if (tp[0] >= 0 && tp[0] < width && tp[1] >= 0 && tp[1] < height) {
                        if (energy[p[1]][p[0]] == 0) energy[p[1]][p[0]] = energy(p[0], p[1]);
                        double en = energy[p[1]][p[0]];
                        if (distTo[tp[1]][tp[0]] > distTo[p[1]][p[0]] + en) {
                            distTo[tp[1]][tp[0]] = distTo[p[1]][p[0]] + en;
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

    public int[] findHorizontalSeam() {
        double minWeight = MAX;
        int height = width(); 
        int width = height(); 
        Stack<Integer> seam = new Stack<>();
        double[][] energy = new double[height][width]; 
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
                    if (t >= 0 && t < width) orderedPoints.add(new int[] { t, r });
            } 
            if (energy[0][col] == 0) energy[0][col] = energy(0, col);
            distTo[0][col] = energy[0][col]; 
            for (int[] p : orderedPoints) { 
                LinkedList<int[]> tpp = new LinkedList<>();
                tpp.add(new int[] { p[0] - 1, p[1] + 1 }); 
                tpp.add(new int[] { p[0], p[1] + 1 }); 
                tpp.add(new int[] { p[0] + 1, p[1] + 1 }); 
                for (int[] tp : tpp) { // 1, 1
                    if (tp[0] >= 0 && tp[0] < width && tp[1] >= 0 && tp[1] < height) {
                        if (energy[p[1]][p[0]] == 0) energy[p[1]][p[0]] = energy(p[1], p[0]);
                        double en = energy[p[1]][p[0]];
                        if (distTo[tp[1]][tp[0]] > distTo[p[1]][p[0]] + en) {
                            distTo[tp[1]][tp[0]] = distTo[p[1]][p[0]] + en;
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
        for (int i : s)
            StdOut.printf("%4d", i);
        StdOut.println();
        StdOut.print("Horizontal Seam:\n");
        int[] h = sc.findHorizontalSeam();
        for (int i : h)
            StdOut.printf("%4d", i);
        StdOut.println();
    }
}
