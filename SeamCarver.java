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
    private static final boolean HORIZONTAL   = true;
    private static final boolean VERTICAL     = false;
    private static final double MAX = Double.POSITIVE_INFINITY;
    private int[][] pixels;

    public SeamCarver(Picture picture) {
        if (picture == null)
            throw new IllegalArgumentException();
        pixels = new int[picture.height()][picture.width()];
        for (int i = 0; i < picture.height(); i++) {
            for (int j = 0; j < picture.width(); j++) {
                pixels[i][j] = picture.getRGB(j, i);
            }
        }
    }

    public Picture picture() {
        Picture rtn = new Picture(pixels[0].length, pixels.length);
        for (int i = 0; i < pixels.length; i++) 
            for (int j = 0; j < pixels[0].length; j++) 
               rtn.setRGB(j, i, pixels[i][j]);
        return rtn;
    }

    public int width() {
        // return picture.width();
        return pixels[0].length;
    }

    public int height() {
        // return picture.height();
        return pixels.length;
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
        if (x < 0 || x >= width() || y < 0 || y >= height()) throw new IllegalArgumentException();
        if (x-1 >= 0 && x + 1 < width() && y - 1 >= 0 && y + 1 < height()) {
            return Math.sqrt(
                calculate(pixels[y][x-1], pixels[y][x+1]) + 
                calculate(pixels[y-1][x], pixels[y+1][x])
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

    private void validateSeam(int[] seam, boolean direction) {
        int length = 0;
        int range = 0;
        if (direction == VERTICAL) {
            length = height();
            range = width();
        } else {
            length = width();
            range = height();
        }
        if (seam == null) throw new IllegalArgumentException("Input is null");
        if (seam.length != length) throw new IllegalArgumentException("Seam is not of right length.");
        int entry = seam[0];
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= range) throw new IllegalArgumentException("Seam is not in range:" + seam[i]);
            if (Math.abs(seam[i] - entry) > 1) throw new IllegalArgumentException("successive entries in seam must differ by -1, 0, or +1");
            entry = seam[i];
        }
    }

    public void removeVerticalSeam(int[] seam) {
        validateSeam(seam, VERTICAL);
        int[][] np = new int[height()][width()-1];
        for (int i = 0; i < pixels.length; i++) {
            int index = seam[i];
            int[] src = pixels[i];
            int[] dest = np[i];
        
            if (index == 0) {
                System.arraycopy(src, 1, dest, 0, dest.length);
            } else if (index == width()-1) {
                System.arraycopy(src, 0, dest, 0, dest.length);
            } else {
                System.arraycopy(src, 0, dest, 0, index);
                System.arraycopy(src, index+1, dest, index, dest.length-index);
            }
        }
        pixels = np;
    }

    public void removeHorizontalSeam(int[] seam) {
        validateSeam(seam, HORIZONTAL);
        int[][] np = new int[height()-1][width()];
        for (int col = 0; col < width(); col++) {
            int nprow = 0;
            for (int row = 0; row < height(); row++) {
                if (row != seam[col]) {
                    np[nprow++][col] = pixels[row][col];
                }
            }
        }
        pixels = np;
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
        // sc.removeVerticalSeam(s);
        sc.removeHorizontalSeam(h);
    }
}
