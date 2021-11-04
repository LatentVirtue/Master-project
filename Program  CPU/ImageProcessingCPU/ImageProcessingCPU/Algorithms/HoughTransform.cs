using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Numerics;
using System.Linq;

namespace ImageProcessingCPU.Algorithms
{
    struct Pixel
    {
        public int x_index;
        public int y_index;

        public double x;
        public double y;

        public Pixel(int x_i, int y_i, double xx, double yy)
        {
            x = xx;
            y = yy;
            x_index = x_i;
            y_index = y_i;
        }
    }

    /*
     struct Pixel
    {
        public int x;
        public int y;
        public bool val;

        public Pixel(int xx, int yy, bool vall)
        {
            x = xx;
            y = yy;
            val = vall;
        }
        public Vector2 ToVector2()
        {
            return new Vector2(x, y);
        }
    }
     */

    struct Kernel
    {
        public double rho;
        public double theta;
        public double[] lambda;
        public int rho_index;
        public int theta_index;
        public LinkedList<Pixel> cluster;

        public double height;
        public Kernel(double rho, double theta, double[] lambda, int rho_index, int theta_index, double height, LinkedList<Pixel> cluster)
        {
            this.rho = rho;
            this.theta = theta;
            this.lambda = lambda;
            this.rho_index = rho_index;
            this.theta_index = theta_index;
            this.height = height;
            this.cluster = cluster;
        }
    }
    class Accumulator
    {
        public int m_image_width = 0;
        public int m_image_height = 0;
        public double m_delta = 0;
        public int m_bins = 0;
        public int m_width = 0;
        public int m_height = 0;
        public double[] m_ro;
        public double roboundsLow = 0;
        public double roboundsHigh = 0;
        public double[] m_theta;
        public double thetaboundsLow = 0;
        public double thetaboundsHigh= 0;
        public double[,] bins;
        public Accumulator(int image_width = 0, int image_height = 0, double delta = 0)
        {
            init(image_width, image_height, delta);
        }
        void init(int image_width, int image_height, double delta)
        {
            if (m_delta != delta || m_image_width != image_width || m_image_height != image_height)
            {
                m_delta = delta;
                m_image_height = image_height;
                m_image_width = image_width;

                double r = Math.Sqrt(image_width * image_width + image_height * image_height);
                m_width = (int)((r + 1) / delta);
                m_ro = new double[m_width + 2];
                m_ro[1] = -0.5 * r;
                for (int i = 2; i <= m_width; ++i)
                {
                    m_ro[i] = m_ro[i - 1] + delta;
                }
                m_ro[0] = m_ro[m_width];
                m_ro[m_width + 1] = m_ro[1];
                roboundsLow = -0.5 * r;
                roboundsHigh = 0.5 * r;

                m_height = (int)(180 / delta);
                m_theta = new double[m_height + 2];
                m_theta[1] = 0.0;
                for (int i = 2; i <= m_height; ++i)
                {
                    m_theta[i] = m_theta[i - 1] + delta;
                }
                m_theta[0] = m_theta[m_height];
                m_theta[m_height + 1] = m_theta[1];
                thetaboundsLow = 0.0;
                thetaboundsHigh = 180 - delta;
                bins = new double[m_height + 2, m_width + 2];
            }
        }
    }
    struct Bin
    {
        public int rho_index;
        public int theta_index;
        public int votes;
        public Bin(int _rho_index = 0, int _theta_index = 0, int _votes = 0)
        {
            rho_index = _rho_index;
            theta_index = _theta_index;
            votes = _votes;
        }
    }
    struct Line
    {
        public double rho;
        public double theta;
        public Line(double _rho, double _theta)
        {
            rho = _rho;
            theta = _theta;
        }
    }
    class VisitedMap
    {
        // Initializes the map.

        // Sets a given accumulator bin as visited.
        public void set_visited(int rho_index, int theta_index)
        {
            m_map[theta_index, rho_index] = true;
        }

        // Class constructor.
        public VisitedMap(int accumulator_width, int accumulator_height)
        {
            if (m_rho_capacity < (accumulator_width + 2) || m_theta_capacity < (accumulator_height + 2))
            {
                m_rho_capacity = accumulator_width + 2;
                m_theta_capacity = accumulator_height + 2;
            }
            m_map = new bool[m_theta_capacity, m_rho_capacity];
        }
        // Returns whether a neighbour bin was visited already.
        public bool Visited_neighbour(int rho_index, int theta_index)
        {
            for (int i = -1; i < 2; i++)
            {
                for (int j = -1; j < 2; j++)
                {
                    if (i == j && j == 0)
                    {
                        continue;
                    }
                    if (m_map[theta_index + i, rho_index + j])
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        // The map of flags ([1,theta_size][1,rho_size] range).
        bool[,] m_map;

        // Specifies the size of allocated storage for the map (rho dimention).
        int m_rho_capacity;

        // Specifies the size of allocated storage for the map (theta dimention).
        int m_theta_capacity;
    }
    class HoughTransform : IAlgoInterface, IDisposable
    {
        Graphics g;
        Canny cannyEdge = new Canny(0, 0.1, 0.3);
        Bitmap actual;
        bool[,] m;
        int cluster_min_size = 5;
        readonly double cluster_min_deviation;
        readonly double delta;
        readonly double kernel_min_height;
        readonly double n_sigmas;
        int[,] gaussKernel = { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } };
        List<LinkedList<Pixel>> chains = new List<LinkedList<Pixel>>();
        List<LinkedList<Pixel>> clusters = new List<LinkedList<Pixel>>();
        Accumulator accumulator;
        List<Line> lines = new List<Line>();
        List<Kernel> used_kernels;
        public HoughTransform(ref Image a, double cluster_min_deviation = 2.0, double delta = 0.5, double kernel_min_height = 0.01, double n_sigmas = 2, Graphics gg = null)
        {
            actual = (Bitmap)a;
            if (gg == null)
            {
                g = Graphics.FromImage(actual);
            }
            else
            {
                g = gg;
            }
            this.cluster_min_deviation = cluster_min_deviation;
            this.delta = delta;
            this.kernel_min_height = kernel_min_height;
            this.n_sigmas = n_sigmas;
            m = cannyEdge.matrixDirect(ref actual);
        }
        public void Dispose()
        {
            g = null;
            cannyEdge = null;
            actual = null;
            m = null;
            lines = null;
            chains = null;
            clusters = null;
            accumulator = null;
        }
        //1. Identify clusters of approx. collinear feature pixels
        //1.a linking
        //1.b subdivision
        /*
        LinkedList<Pixel> LinkOld(Pixel pr)
        {
            LinkedList<Pixel> S = new LinkedList<Pixel>();
            Pixel p = pr;
            while (p.x != -1)
            {
                S.AddLast(p);
                m[p.y, p.x] = false;
                p = Next(p);
            }
            p = Next(pr);
            while (p.x != -1)
            {
                S.AddFirst(p);
                m[p.y, p.x] = false;
                p = Next(p);
            }
            return S;
        }*/
        //Complementing Link()
        /*Pixel NextOld(Pixel ps)
        {
            for (int i = -1; i < 2; i++)
            {
                for (int j = -1; j < 2; j++)
                {
                    if (i == 0 && j == 0)
                    {
                        continue;
                    }
                    if (i + ps.y < 0 || i + ps.y >= m.GetLength(0) || j + ps.x < 0 || j + ps.x >= m.GetLength(1) || (i == 0 && j == 0))
                    {
                        continue;
                    }
                    if (m[i + ps.y, j + ps.x])
                    {
                        return new Pixel(j + ps.x, i + ps.y, m[i + ps.y, j + ps.x]);
                    }
                }
            }
            return new Pixel(-1, -1, false);
        }*/
        int[] X_OFFSET = { 0, 1, 0, -1, 1, -1, -1, 1 };
        int[] Y_OFFSET = { 1, 0, -1, 0, 1, 1, -1, -1 };
        bool Next(ref int x_seed, ref int y_seed)
        {
            for (int i = 0; i != 8; ++i)
            {
                int x = x_seed + X_OFFSET[i];
                if (0 <= x && x < actual.Width)
                {
                    int y = y_seed + Y_OFFSET[i];
                    if (0 <= y && y < actual.Height)
                    {
                        if (m[x, y])
                        {
                            m[x, y] = false;
                            x_seed = x;
                            y_seed = y;
                            return true;
                        }
                    }
                }
            }
            return false;
        }
        void Link(ref LinkedList<Pixel> chain, int x_ref, int y_ref)
        {
            chain.Clear();
            int x = x_ref;
            int y = y_ref;

            do
            {
                //Pixel t = new Pixel(x, y, x - actual.Width / 2.0, y - actual.Height / 2.0);
                Pixel t = new Pixel(y, x, y - actual.Height / 2.0, x - actual.Width / 2.0);
                chain.AddFirst(t);
            } while (Next(ref x, ref y));

            x = x_ref;
            y = y_ref;
            if (Next(ref x, ref y))
            {
                do
                {
                    //Pixel t = new Pixel(x, y, x - actual.Width / 2.0, y - actual.Height / 2.0);
                    Pixel t = new Pixel(y, x, y - actual.Height / 2.0, x - actual.Width / 2.0);
                    chain.AddLast(t);
                } while (Next(ref x, ref y));
            }
        }

        void Find_chains()
        {
            for (int y = 1, y_end = actual.Height - 1; y != y_end; ++y)
            {
                for (int x = 1, x_end = actual.Width - 1; x != x_end; ++x)
                {
                    if (m[x, y])
                    {
                        LinkedList<Pixel> chain = new LinkedList<Pixel>();
                        Link(ref chain, x, y);
                        if (chain.Count >= cluster_min_size)
                        {
                            chains.Add(chain);
                        }
                    }
                }
            }
        }

        /*
        void Find_chainsOld(int image_width, int image_height, int min_size)
        {
            for (int y = 0; y < image_height; y++)
            {
                for (int x = 0; x < image_width; x++)
                {
                    if (m[y, x])
                    {
                        var t = Link(new Pixel(x, y, true));
                        if (t.Count > min_size)
                        {
                            chains.Add(t);
                        }
                    }
                }
            }
        }*/
        double Subdivision_procedure(ref List<LinkedList<Pixel>> clusters, ref LinkedList<Pixel> chain, int first_index, int last_index)
        {
            int clusters_count = clusters.Count;
            if (chain.Count == 0)
            {
                return 0;
            }
            Pixel first = chain.First.Value;
            Pixel last = chain.Last.Value;

            // Compute the length of the straight line segment defined by the endpoints of the cluster.
            int y = first.y_index - last.y_index;
            int x = first.x_index - last.x_index;
            double length = Math.Sqrt(x * x + y * y);

            // Find the pixels with maximum deviation from the line segment in order to subdivide the cluster.
            int max_pixel_index = 0;
            double deviation, max_deviation = -1.0;
            LinkedListNode<Pixel> t = chain.First;
            if (last_index > chain.Count)
            {
                //last_index = chain.Count;
            }
            for (int i = first_index, count = chain.Count; i != last_index; i = (i + 1) % count)
            {
                Pixel current = t.Value;
                t = t.Next;

                deviation = Math.Abs((current.x - first.x) * (first.y - last.y) + (current.y - first.y) * (last.x - first.x));

                if (deviation > max_deviation)
                {
                    max_pixel_index = i;
                    max_deviation = deviation;
                }
            }
            max_deviation /= length;

            // Compute the ratio between the length of the segment and the maximum deviation.
            double ratio = length / Math.Max(max_deviation, cluster_min_deviation);
            if (Math.Abs(first_index - last_index) < 2)
            {
                return ratio;
            }
            // Test the number of pixels of the sub-clusters.
            if ((max_pixel_index - first_index + 1) >= cluster_min_size && (last_index - max_pixel_index + 1) >= cluster_min_size)
            {
                double ratio1 = Subdivision_procedure(ref clusters, ref chain, first_index, max_pixel_index);
                double ratio2 = Subdivision_procedure(ref clusters, ref chain, max_pixel_index, last_index);

                // Test the quality of the sub-clusters against the quality of the current cluster.
                if (ratio1 > ratio || ratio2 > ratio)
                {
                    return Math.Max(ratio1, ratio2);
                }
            }
            // Remove the sub-clusters from the list of clusters.
            clusters.RemoveRange(clusters_count, clusters.Count - clusters_count);
            // Keep current cluster
            LinkedList<Pixel> temp = new LinkedList<Pixel>();
            var c = chain.First;
            for (int i = 1; i < chain.Count; i++)
            {
                if (i >= first_index)
                {
                    if (i > last_index)
                    {
                        break;
                    }
                    temp.AddLast(c.Value);
                }
                c = c.Next;
            }
            clusters.Add(temp);
            return ratio;
        }
        void Find_clusters()
        {
            for (int i = 0; i < chains.Count; i++)
            {
                var chain = chains[i];
                Subdivision_procedure(ref clusters, ref chain, 0, chains[i].Count - 1);
                chains[i] = chain;
            }
        }
        //2. Compute an elliptical kernel for each cluster from its line fitting uncertainty

        Kernel ComputeKernelOld(LinkedList<Pixel> cluster, double rho_max, double one_div_delta)
        {
            double one_div_npixels = 1 / (double)(cluster.Count);

            // Alternative reference system definition.
            double mean_x = 0.0;
            double mean_y = 0.0;
            foreach (Pixel p in cluster)
            {
                mean_x += p.x;
                mean_y += p.y;
            }
            mean_x *= one_div_npixels;
            mean_y *= one_div_npixels;
            double[,] M = new double[2, 2];
            foreach (Pixel p in cluster)
            {
                double x = p.x - mean_x;
                double y = p.y - mean_y;

                M[0, 0] += x * x;
                M[1, 1] += y * y;
                M[0, 1] += x * y;
            }
            M[1, 0] = M[0, 1];

            //eigen(V, S, M);
            double[] temp = eigenValue(M);
            Vector2 t1 = eigenVector(M, temp[0]);
            Vector2 t2 = eigenVector(M, temp[1]);
            Vector2 u = temp[0] > temp[1] ? t1 : t2;
            Vector2 v = temp[0] > temp[1] ? t2 : t1;
            // y_v >= 0 condition verification.
            if (v.Y < 0.0)
            {
                v.X *= -1;
                v.Y *= -1;
            }

            // Normal equation parameters computation 
            double rho = v.X * mean_x + v.Y * mean_y;
            double theta = Math.Acos(v.X) * (180 / Math.PI);
            int rho_index = (int)(Math.Abs((rho + rho_max) * one_div_delta)) + 1;
            int theta_index = (int)(Math.Abs(theta * one_div_delta)) + 1;

            // sigma^2_m' and sigma^2_b' computation
            double aux = Math.Sqrt(1.0 - v.X * v.X);
            double[] nabla = { -(u.X * mean_x + u.Y * mean_y), 1.0, aux != 0.0 ? (u.X / aux) * (180 / Math.PI) : 0.0, 0.0 };

            aux = 0.0;
            foreach (Pixel p in cluster)
            {
                double x = u.X * (p.x - mean_x) + u.Y * (p.y - mean_y);
                aux += x * x;
            }

            double[] lambda = { 1.0 / aux, 0.0, 0.0, one_div_npixels };

            // Uncertainty from sigma^2_m' and sigma^2_b' to sigma^2_rho,  sigma^2_theta and sigma_rho_theta.
            Solve(ref lambda, ref nabla, ref lambda);

            if (lambda[3] == 0.0)
            {
                lambda[3] = 0.1;
            }

            lambda[0] *= n_sigmas * n_sigmas;
            lambda[3] *= n_sigmas * n_sigmas;

            // Compute the height of the kernel.
            double height = Gauss(0.0, 0.0, lambda[0], lambda[3], lambda[1]);

            // Keep kernel.
            return new Kernel(rho, theta, lambda, rho_index, theta_index, height, cluster);
        }

        Kernel ComputeKernel(LinkedList<Pixel> cluster, double rho_max, double one_div_delta)
        {
            double one_div_npixels = 1 / (double)(cluster.Count);

            // Alternative reference system definition.
            double mean_x = 0.0;
            double mean_y = 0.0;
            foreach (Pixel p in cluster)
            {
                mean_x += p.x;
                mean_y += p.y;
            }
            mean_x *= one_div_npixels;
            mean_y *= one_div_npixels;
            double[] M = new double[4];
            double[] V = new double[4];
            double[] S = new double[4];
            foreach (Pixel p in cluster)
            {
                double x = p.x - mean_x;
                double y = p.y - mean_y;

                M[0] += x * x;
                M[3] += y * y;
                M[1] += x * y;
            }
            M[2] = M[1];

            eigen(ref V, ref S, ref M);
            // y_v >= 0 condition verification.
            if (V[3] < 0.0)
            {
                V[1] *= -1.0;
                V[3] *= -1.0;
            }

            // Normal equation parameters computation 
            double rho = V[1] * mean_x + V[3] * mean_y;
            double theta = Math.Acos(V[1]) * (180 / Math.PI);
            int rho_index = (int)(Math.Abs((rho + rho_max) * one_div_delta)) + 1;
            int theta_index = (int)(Math.Abs(theta * one_div_delta)) + 1;

            // sigma^2_m' and sigma^2_b' computation
            double aux = Math.Sqrt(1.0 - V[1] * V[3]);
            double[] nabla = { -(V[0] * mean_x + V[2] * mean_y), 1.0, aux != 0.0 ? (V[0] / aux) * (180 / Math.PI) : 0.0, 0.0 };

            aux = 0.0;
            foreach (Pixel p in cluster)
            {
                double x = V[0] * (p.x - mean_x) + (V[2] * (p.y - mean_y));
                aux += x * x;
            }

            double[] lambda = { 1.0 / aux, 0.0, 0.0, one_div_npixels };

            // Uncertainty from sigma^2_m' and sigma^2_b' to sigma^2_rho,  sigma^2_theta and sigma_rho_theta.
            Solve(ref lambda, ref nabla, ref lambda);

            if (lambda[3] == 0.0)
            {
                lambda[3] = 0.1;
            }

            lambda[0] *= n_sigmas * n_sigmas;
            lambda[3] *= n_sigmas * n_sigmas;

            // Compute the height of the kernel.
            double height = Gauss(0.0, 0.0, lambda[0], lambda[3], lambda[1]);

            // Keep kernel.
            return new Kernel(rho, theta, lambda, rho_index, theta_index, height, cluster);
        }
        //helper function for 2x2 matrix addition / subtraction
        int[,] TbTOperations(int[,] a, int[,] b, char operation)
        {
            int[,] r = new int[2, 2];
            switch (operation)
            {
                case '-':
                    r[0, 0] = a[0, 0] - b[0, 0];
                    r[0, 1] = a[0, 1] - b[0, 1];
                    r[1, 0] = a[1, 0] - b[1, 0];
                    r[1, 1] = a[1, 1] - b[1, 1];
                    return r;
                case '+':
                    r[0, 0] = a[0, 0] + b[0, 0];
                    r[0, 1] = a[0, 1] + b[0, 1];
                    r[1, 0] = a[1, 0] + b[1, 0];
                    r[1, 1] = a[1, 1] + b[1, 1];
                    return r;
                default:
                    return null;
            }
        }
        //helper function for 2x2 matrix composition out of 2 columns
        int[,] multiply(int x1, int y1, int x2, int y2)
        {
            int[,] r = new int[2, 2];
            r[0, 0] = x1 * x2;
            r[0, 1] = x1 * y2;
            r[1, 0] = y1 * x2;
            r[1, 1] = y1 * y2;
            return r;
        }
        //helper function for finding eigenvalues of a 2x2 matrix
        float[] eigenValue(float[,] m)
        {
            float b = -(m[0, 0] + m[1, 1]);
            float c = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
            float[] r = new float[2];
            r[0] = (float)(b - Math.Sqrt(b * b - 4 * c));
            r[1] = (float)(b + Math.Sqrt(b * b - 4 * c));
            return r;
        }
        //alternative eigen decomposition
        void tri_diagonalize(ref double[] Cxd, ref double[] d, ref double[] e, ref double[] A, int L, double tol)
        {
            int i, j, k, l;
            double f, g, h, hh;
            for (i = 0; i < L; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    A[i * L + j] = Cxd[i * L + j];
                }
            }
            for (i = L - 1; i > 0; i--)
            {
                l = i - 2;
                f = A[i * L + i - 1];
                g = 0.0;
                for (k = 0; k <= l; k++)
                {
                    g += A[i * L + k] * A[i * L + k];
                }
                h = g + f * f;
                if (g <= tol)
                {
                    e[i] = f;
                    h = 0.0;
                    d[i] = h;
                    continue;
                }
                l++;
                g = Math.Sqrt(h);
                if (f >= 0.0) g = -g;
                e[i] = g;
                h = h - f * g;
                A[i * L + i - 1] = f - g;
                f = 0.0;
                for (j = 0; j <= l; j++)
                {
                    A[j * L + i] = A[i * L + j] / h;
                    g = 0.0;
                    for (k = 0; k <= j; k++)
                    {
                        g += A[j * L + k] * A[i * L + k];
                    }
                    for (k = j + 1; k <= l; k++)
                    {
                        g += A[k * L + j] * A[i * L + k];
                    }
                    e[j] = g / h;
                    f += g * A[j * L + i];
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++)
                {
                    f = A[i * L + j];
                    g = e[j] - hh * f;
                    e[j] = g;
                    for (k = 0; k <= j; k++)
                    {
                        A[j * L + k] = A[j * L + k] - f * e[k] - g * A[i * L + k];
                    }
                }
                d[i] = h;
            }
            d[0] = e[0] = 0.0;
            for (i = 0; i < L; i++)
            {
                l = i - 1;
                if (d[i] != 0.0)
                {
                    for (j = 0; j <= l; j++)
                    {
                        g = 0.0;
                        for (k = 0; k <= l; k++)
                        {
                            g += A[i * L + k] * A[k * L + j];
                        }
                        for (k = 0; k <= l; k++)
                        {
                            A[k * L + j] = A[k * L + j] - g * A[k * L + i];
                        }
                    }
                }
                d[i] = A[i * L + i];
                A[i * L + i] = 1.0;
                for (j = 0; j <= l; j++)
                {
                    A[i * L + j] = A[j * L + i] = 0.0;
                }
            }
        }
        int calc_eigenstructure(ref double[] d, ref double[] e, ref double[] A, int L, double macheps)
        {
            int i, j, k, l, m;
            double b, c, f, g, h, p, r, s;

            for (i = 1; i < L; i++) e[i - 1] = e[i];
            e[L - 1] = b = f = 0.0;
            for (l = 0; l < L; l++)
            {
                h = macheps * (Math.Abs(d[l]) + Math.Abs(e[l]));
                if (b < h) b = h;
                for (m = l; m < L; m++)
                {
                    if (Math.Abs(e[m]) <= b) break;
                }
                j = 0;
                if (m != l)
                {
                    do
                    {
                        if (j++ == 30) return -1;
                        p = (d[l + 1] - d[l]) / (2.0 * e[l]);
                        r = Math.Sqrt(p * p + 1);
                        h = d[l] - e[l] / (p + (p < 0.0 ? -r : r));
                        for (i = l; i < L; i++) d[i] = d[i] - h;
                        f += h;
                        p = d[m];
                        c = 1.0;
                        s = 0.0;
                        for (i = m - 1; i >= l; i--)
                        {
                            g = c * e[i];
                            h = c * p;
                            if (Math.Abs(p) >= Math.Abs(e[i]))
                            {
                                c = e[i] / p;
                                r = Math.Sqrt(c * c + 1);
                                e[i + 1] = s * p * r;
                                s = c / r;
                                c = 1.0 / r;
                            }
                            else
                            {
                                c = p / e[i];
                                r = Math.Sqrt(c * c + 1);
                                e[i + 1] = s * e[i] * r;
                                s = 1.0 / r;
                                c = c / r;
                            }
                            p = c * d[i] - s * g;
                            d[i + 1] = h + s * (c * g + s * d[i]);
                            for (k = 0; k < L; k++)
                            {
                                h = A[k * L + i + 1];
                                A[k * L + i + 1] = s * A[k * L + i] + c * h;
                                A[k * L + i] = c * A[k * L + i] - s * h;
                            }
                        }
                        e[l] = s * p;
                        d[l] = c * p;
                    }
                    while (Math.Abs(e[l]) > b);
                }
                d[l] = d[l] + f;
            }

            /* order the eigenvectors  */
            for (i = 0; i < L; i++)
            {
                k = i;
                p = d[i];
                for (j = i + 1; j < L; j++)
                {
                    if (d[j] > p)
                    {
                        k = j;
                        p = d[j];
                    }
                }
                if (k != i)
                {
                    d[k] = d[i];
                    d[i] = p;
                    for (j = 0; j < L; j++)
                    {
                        p = A[j * L + i];
                        A[j * L + i] = A[j * L + k];
                        A[j * L + k] = p;
                    }
                }
            }
            return 0;
        }

        // Computes the decomposition of a matrix into matrices composed of its eigenvectors and eigenvalues.
        void eigen(ref double[] vectors, ref double[] values, ref double[] matrix)
        {
            double[] temp = new double[2];

            tri_diagonalize(ref matrix, ref values, ref temp, ref vectors, 2, 1.0e-6);
            calc_eigenstructure(ref values, ref temp, ref vectors, 2, 1.0e-16);

            values[3] = values[1];
            values[1] = values[2] = 0.0;
        }

        double[] eigenValue(double[,] m)
        {
            double b = m[0, 0] + m[1, 1];
            double c = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
            double[] r = new double[2];
            r[0] = -(b + Math.Sqrt(b * b - 4 * c)) / 2;
            r[1] = -(b - Math.Sqrt(b * b - 4 * c)) / 2;
            return r;
        }
        Vector2 eigenVector(float[,] m, float eV)
        {
            return new Vector2((-(m[0, 1]) / (m[0, 0] - eV)), 1);
        }
        Vector2 eigenVector(double[,] m, double eV)
        {
            return new Vector2((float)(-(m[0, 1]) / (m[0, 0] - eV)), 1);
        }
        int convolution(ref double[,] bins, int rho_index, int theta_index)
        {
            double t = 0;
            long offset = -1;
            int n = bins.GetLength(0);
            int m = bins.GetLength(1);
            for (long i = 2; i > -1; i--)
            {
                for (long j = 2; j > -1; j--)
                {
                    long x = theta_index + i + offset;
                    double tt;
                    bool u = true;
                    if (x >= n || x < 0)
                    {
                        u = false;
                    }
                    long y = rho_index + j + offset;
                    if (y >= m || y < 0)
                    {
                        u = false;
                    }
                    if (u)
                    {
                        tt = bins[x, y];
                    }
                    else
                    {
                        tt = bins[theta_index, rho_index];
                    }
                    t += gaussKernel[i, j] * tt;
                }
            }
            return (int)Math.Round(t);
        }
        //3. Vote for kernels with bigger contributions
        //3.a culling
        //3.b voting
        void Solve(ref double[] result, ref double[] nabla, ref double[] lambda)
        {
            double[] temp = new double[4];
            for (int i = 0, i_line = 0; i < 2; ++i, i_line += 2)
            {
                for (int j = 0; j < 2; ++j)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        temp[i_line + j] += nabla[i_line + k] * lambda[k * 2 + j];
                    }
                }
            }

            result = new double[4];
            for (int i = 0, i_line = 0; i < 2; ++i, i_line += 2)
            {
                for (int j = 0; j < 2; ++j)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        result[i_line + j] += temp[i_line + k] * nabla[j * 2 + k];
                    }
                }
            }
        }
        double Gauss(double rho, double theta, double sigma2_rho, double sigma2_theta, double sigma_rho_sigma_theta, double two_r, double a, double b)
        {
            return a * Math.Exp(-(((rho * rho) / sigma2_rho) - ((two_r * rho * theta) / sigma_rho_sigma_theta) + ((theta * theta) / sigma2_theta)) * b);
        }
        double Gauss(double rho, double theta, double sigma2_rho, double sigma2_theta, double sigma_rho_theta)
        {

            //Equation 15

            double sigma_rho_sigma_theta = Math.Sqrt(sigma2_rho * sigma2_theta);
            double r = (sigma_rho_theta / sigma_rho_sigma_theta), two_r = 2.0 * r;
            double a = 1.0 / (2.0 * Math.PI * sigma_rho_sigma_theta * Math.Sqrt(1.0 - r * r));
            double b = 1.0 / (2.0 * (1.0 - r * r));
            return Gauss(rho, theta, sigma2_rho, sigma2_theta, sigma_rho_sigma_theta, two_r, a, b);
        }
        void Vote(int rho_start_index, int theta_start_index, double rho_start, double theta_start, int inc_rho_index, int inc_theta_index, double sigma2_rho, double sigma2_theta, double sigma_rho_theta, double scale)
        {
            //int rho_size = accumulator.m_width;
            //int theta_size = accumulator.m_height;
            //double delta = accumulator.m_delta;
            var bins = accumulator.bins;
            double inc_rho = delta * inc_rho_index, inc_theta = delta * inc_theta_index;

            double sigma_rho_sigma_theta = Math.Sqrt(sigma2_rho * sigma2_theta);
            double r = (sigma_rho_theta / sigma_rho_sigma_theta);
            double a = 1.0 / (2.0 * Math.PI * sigma_rho_sigma_theta * Math.Sqrt(1.0 - r * r));
            double b = 1.0 / (2.0 * (1.0 - r * r));

            bool theta_voted;
            double rho, theta;
            int votes, theta_not_voted = 0;
            long rho_index, theta_index, theta_count = 0;

            // Loop for the theta coordinates of the parameter space.
            theta_index = theta_start_index;
            theta = theta_start;

            do
            {
                // Test if the kernel exceeds the parameter space limits.
                if (theta_index == 0 || theta_index == (accumulator.m_height + 1))
                {
                    rho_start_index = accumulator.m_width - rho_start_index + 1;
                    theta_index = theta_index == 0 ? accumulator.m_height : 1;
                    inc_rho_index = -inc_rho_index;
                }

                // Loop for the rho coordinates of the parameter space.
                theta_voted = false;

                rho_index = rho_start_index;
                rho = rho_start;
                while ((votes = (int)Math.Round((Gauss(rho, theta, sigma2_rho, sigma2_theta, sigma_rho_sigma_theta, r * 2, a, b) * scale) + 0.5)) > 0 && rho_index >= 1 && rho_index <= accumulator.m_width)
                {
                    accumulator.bins[theta_index, rho_index] += votes;
                    theta_voted = true;

                    rho_index += inc_rho_index;
                    rho += inc_rho;
                }
                if (!theta_voted)
                {
                    theta_not_voted++;
                }

                theta_index += inc_theta_index;
                theta += inc_theta;
                theta_count++;
            }
            while (theta_not_voted != 2 && theta_count < accumulator.m_height);
        }
        void Voting(double kernel_min_height)
        {
            List<Kernel> kernels;


            kernels = new List<Kernel>();
            used_kernels = new List<Kernel>();

            foreach (LinkedList<Pixel> cluster in clusters)
            {
                kernels.Add(ComputeKernel(cluster, accumulator.roboundsHigh, 1.0 / accumulator.m_delta));
            }

            // Discard groups with very short kernels.
            double norm = double.MinValue;
            foreach (Kernel kernel in kernels)
            {
                if (norm < kernel.height)
                {
                    norm = kernel.height;
                }
                used_kernels.Add(kernel);
            }
            norm = 1.0 / norm;
            int i = 0;
            for (int k = 0; k < used_kernels.Count; ++k)
            {
                if (kernels[k].height * norm >= kernel_min_height)
                {
                    if (i != k)
                    {
                        Kernel temp = used_kernels[i];
                        used_kernels[i] = used_kernels[k];
                        used_kernels[k] = temp;
                    }
                    i++;
                }
            }
            if (i < used_kernels.Count)
            {
                used_kernels.RemoveRange(i, used_kernels.Count - i);
            }
            // Find the g_min threshold and compute the scale factor for integer votes.
            double kernels_scale = double.MinValue;
            for (int j = 0; j < used_kernels.Count; j++)
            {
                double[] V = new double[4];
                double[] S = new double[4];
                double[] lambda = used_kernels[j].lambda;
                eigen(ref V, ref S, ref lambda);

                double radius = Math.Sqrt(S[3]);
                double scale = Gauss(V[2] * radius, V[3] * radius, lambda[0], lambda[3], lambda[1]);
                scale = scale < 1.0 ? (1.0 / scale) : 1.0;
                used_kernels[j].lambda[0] = lambda[0];
                used_kernels[j].lambda[1] = lambda[1];
                used_kernels[j].lambda[2] = lambda[2];
                used_kernels[j].lambda[3] = lambda[3];
                if (kernels_scale < scale)
                {
                    kernels_scale = scale;
                }
            }

            // Vote for each selected kernel.
            foreach (Kernel kernel in used_kernels)
            {
                Vote(kernel.rho_index, kernel.theta_index, 0.0, 0.0, 1, 1, kernel.lambda[0], kernel.lambda[3], kernel.lambda[1], kernels_scale);
                Vote(kernel.rho_index, kernel.theta_index - 1, 0.0, -delta, 1, -1, kernel.lambda[0], kernel.lambda[3], kernel.lambda[1], kernels_scale);
                Vote(kernel.rho_index - 1, kernel.theta_index, -delta, 0.0, -1, 1, kernel.lambda[0], kernel.lambda[3], kernel.lambda[1], kernels_scale);
                Vote(kernel.rho_index - 1, kernel.theta_index - 1, -delta, -delta, -1, -1, kernel.lambda[0], kernel.lambda[3], kernel.lambda[1], kernels_scale);
            }
        }
        void Peak_detection()
        {
            // Create a list with all cells that receive at least one vote.
            List<Bin> used_bins = new List<Bin>();

            for (int i = 0; i < accumulator.m_height; i++)
            {
                for (int j = 0; j < accumulator.m_width; j++)
                {
                    if (accumulator.bins[i, j] != 0)
                    {
                        //used_bins.Add(new Bin(j, i, convolution(ref accumulator.bins, j, i))); // Convolution of the cells with a 3x3 Gaussian kernel
                        used_bins.Add(new Bin(j, i, (int)Math.Round(accumulator.bins[i, j])));
                    }
                }
            }

            // Sort the list in descending order according to the result of the convolution.
            used_bins.Sort((x, y) => -x.votes.CompareTo(y.votes));
            // Use a sweep plane that visits each cell of the list.
            VisitedMap visited = new VisitedMap(accumulator.m_width, accumulator.m_height);

            foreach (Bin bin in used_bins)
            {
                if (!visited.Visited_neighbour(bin.rho_index, bin.theta_index))
                {
                    lines.Add(new Line(accumulator.m_ro[bin.rho_index], accumulator.m_theta[bin.theta_index]));
                }
                visited.set_visited(bin.rho_index, bin.theta_index);
            }

        }
        public void DrawLines(List<Line> lines, int lineWidth)
        {
            //Point c = new Point(actual.Width / 2, actual.Height / 2);
            int count = 0;
            foreach (Line l in lines)
            {
                if (count >= 25)
                {
                    break;
                }
                count++;
                Pen yPen = new Pen(Color.Yellow, lineWidth);
                double rho = l.rho;
                double theta = l.theta * (Math.PI / 180);
                double cos_theta = Math.Cos(theta), sin_theta = Math.Sin(theta);
                PointF p1 = new Point();
                PointF p2 = new Point();

                if (sin_theta != 0)
                {
                    p1.X = -actual.Width / 2.0f;
                    p1.Y = (float)((rho - p1.X * cos_theta) / sin_theta);
                    p2.X = actual.Width / 2.0f - 1;
                    p2.Y = (float)((rho - p2.X * cos_theta) / sin_theta);
                }
                else
                {
                    p1.X = (float)rho;
                    p1.Y = -actual.Height;
                    p2.X = (float)rho;
                    p2.Y = (float)actual.Height * 0.5f - 1;
                }
                p1.X += actual.Width * 0.5f;
                p1.Y += actual.Height * 0.5f;
                p2.X += actual.Width * 0.5f;
                p2.Y += actual.Height * 0.5f;
                //Point p = new Point((int)(l.rho * Math.Cos(l.theta)), (int)(l.rho * Math.Sin(l.theta)));
                //g.DrawLine(yPen, new Point(p.Y,p.X), p);
                PointF[] niz = { p1, p2 };
                g.DrawLine(yPen, p1, p2);
                //g.DrawCurve(yPen, niz);


            }
        }
        public Image Apply(ref Image x)
        {
            accumulator = new Accumulator(x.Width, x.Height, delta);
            Find_chains();
            Find_clusters();
            //printing clusters for testing purposes
            /*
            Random r = new Random();
            foreach (LinkedList<Pixel> l in clusters)
            {
                Color t = Color.FromArgb((byte)r.Next(0, 255), (byte)r.Next(0, 255), (byte)r.Next(0, 255));
                foreach (Pixel p in l)
                {
                    actual.SetPixel(p.x, p.y, t);
                }
            }
            */

            Voting(kernel_min_height);
            Peak_detection();
            /*
            foreach(Line li in lines)
            {
                double rUp = li.rho * 1.2;
                double rDown = li.rho * 0.8;
                double tUp = li.theta * 1.2;
                double tDown = li.theta * 0.8;
                foreach (Kernel k in used_kernels)
                {
                    if (rDown < k.rho && k.rho < rUp && tDown < k.theta && k.theta < tUp)
                    {
                        Color t = Color.Red;
                        foreach (Pixel p in k.cluster)
                        {
                            actual.SetPixel(p.x, p.y, t);
                        }

                        //var p1 = k.cluster.First.Value;
                        //var p2 = k.cluster.Last.Value;
                        //g.DrawLine(new Pen(Color.Red, 2), p1.x, p1.y, p2.x, p2.y);
                        
                    }
                }
            }
            */
            DrawLines(lines, 1);
            return actual;
        }
    }
}
