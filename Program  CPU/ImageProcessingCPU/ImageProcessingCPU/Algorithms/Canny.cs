using System;
using System.Collections.Generic;
using System.Drawing;
using ImageProcessor.Imaging.Filters.Photo;
using System.Linq;

namespace ImageProcessingCPU.Algorithms
{
    class HelperObject
    {
        //public int label; //in case i need it for conected component labeling or hoshen-kopelman algo
        public bool hasStrong = false;
        public List<Point> body = new List<Point>();
        public bool upIntersect = false;
        public bool downIntersect = false;
        //this will be the magical BFS
        //need to implement levels
        public HelperObject(ref byte[,] x, int i, int j)
        {
            body.Add(new Point(i, j));
            Queue<Point> q = new Queue<Point>();
            q.Enqueue(new Point(i, j));
            x[i, j] += 10;
            while (q.Count > 0)
            {
                Point t = q.Dequeue();
                for (int ii = Math.Max(t.X - 1, 0); ii < Math.Min(t.X + 2, x.GetLength(0)); ii++)
                {
                    for (int jj = Math.Max(t.Y - 1, 0); jj < Math.Min(t.Y + 2, x.GetLength(1)); jj++)
                    {
                        if (x[ii, jj] < 10)
                        {
                            if (x[ii, jj] > 0)
                            {
                                Point n = new Point(ii, jj);
                                q.Enqueue(n);
                                body.Add(n);
                                if (!hasStrong && x[ii, jj] == 2)
                                {
                                    hasStrong = true;
                                }
                            }
                            x[ii, jj] += 10;
                        }
                    }
                }
            }
            //magical LINQ sorter for fixing the fact that I changed body from List<Point[]> to List<Point> because I'm out of time (mora trcim sa skoletom u nemir)
            //dont ask, it works
            //body = body.OrderByDescending(p => p.X).ThenBy(p => p.Y).ToList();
        }
        //this will set all pixels within the object to 0 and visited
        public void MutuallyAssuredDestruction(ref byte[,] x)
        {
            //kill everyone else
            //edit for point lists
            foreach (Point p in body)
            {
                x[p.X, p.Y] = 10;
            }
            //kill self
            body = null;
        }
    }
    class Canny : IAlgoInterface
    {
        //Sobel operator kernels
        readonly int[,] sobelX = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
        readonly int[,] sobelY = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };

        //Prewitt
        readonly int[,] prewittX = { { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 } };
        readonly int[,] prewittY = { { 1, 1, 1 }, { 0, 0, 0 }, { -1, -1, -1 } };

        //Roberts
        readonly int[,] robertsX = { { 1, 0 }, { 0, -1 } };
        readonly int[,] robertsY = { { 0, 1 }, { -1, 0 } };

        //Scharr
        readonly int[,] scharrX = { { 3, 0, -3 }, { 10, 0, -10 }, { 3, 0, -3 } };
        readonly int[,] scharrY = { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } };
        //actual
        int[,] opX;
        int[,] opY;
        //optimize these for memory
        Image actual;
        double[,] angle;
        double[,] gIntensity;
        byte[,] label;
        double max;
        readonly double lowerT = 0.1;
        readonly double upperT = 0.3;
        public Canny(int kOperator = 0)
        {
            switch (kOperator)
            {
                case 0:
                    opX = sobelX;
                    opY = sobelY;
                    break;
                case 1:
                    opX = prewittX;
                    opY = prewittY;
                    break;
                case 2:
                    opX = robertsX;
                    opY = robertsY;
                    break;
                case 3:
                    opX = scharrX;
                    opY = scharrY;
                    break;
            }
        }
        void ToGrayscale()
        {
            ImageHandler.factory.Load(actual);
            ImageHandler.factory.Filter(MatrixFilters.GreyScale);
            actual = ImageHandler.factory.Image;
        }
        //1. Gaussian filter
        void GaussianFilter()
        {
            ImageHandler.factory.Load(actual);
            ImageHandler.factory.GaussianBlur(5);
            actual = ImageHandler.factory.Image;
        }
        //2. Intensity gradient
        void Gradient()
        {
            //applying kernels
            Bitmap temp = new Bitmap(actual);
            int[,] Gx = Convolve(ref temp, ref opX);
            int[,] Gy = Convolve(ref temp, ref opY);
            //computing gradient intensity via G = sqrt(Gx^2+Gy^2)
            ComputeGradient(ref Gx, ref Gy);
            //actual = temp;
            ImageHandler.factory.Load(actual);
        }
        //also implement channels
        void ComputeGradient(ref int[,] Gx, ref int[,] Gy)
        {
            gIntensity = new double[Gx.GetLength(0), Gx.GetLength(1)];
            angle = new double[Gx.GetLength(0), Gx.GetLength(1)];
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    gIntensity[i, j] = Math.Sqrt(Gx[i, j] * Gx[i, j] + Gy[i, j] * Gy[i, j]);
                    angle[i, j] = (Math.Atan2(Gx[i, j], Gy[i, j]) * (180 / Math.PI)) % 180;
                }
            }
            max = gIntensity.Cast<double>().Max();
        }
        //also add an overload for color channels
        int[,] Convolve(ref Bitmap target, ref int[,] filter)
        {
            int[,] ret = new int[target.Height, target.Width];
            for (int i = 0; i < target.Height; i++)
            {
                for (int j = 0; j < target.Width; j++)
                {
                    //here
                    int sumR = 0;
                    int offX = (filter.GetLength(0) - 1) / -2;
                    int offY = (filter.GetLength(1) - 1) / -2;
                    for (int a = filter.GetLength(0) - 1; a >= 0; a--)
                    {
                        for (int b = filter.GetLength(1) - 1; b >= 0; b--)
                        {
                            int ci = i + offY + a;
                            int cj = j + offX + b;
                            Color t;
                            if (ci < 0 || ci >= target.Height || cj < 0 || cj >= target.Width)
                            {
                                t = target.GetPixel(j, i);
                            }
                            else
                            {
                                t = target.GetPixel(cj, ci);
                            }
                            sumR += t.R * filter[a, b];
                        }
                    }
                    ret[i, j] = sumR;
                }
            }
            return ret;
        }
        //3. Apply gradient magnitude tresholding or lower-bound cutoff suppression to get rid of spurious response to edge detection
        //non-maximum suppresion
        //also works only with greyscale. consider overloading for color channels
        //adjust this, it eats too many pixels for no reason
        void NonMaxSuppression()
        {
            //inner image
            double[,] x = new double[gIntensity.GetLength(0), gIntensity.GetLength(1)];
            for (int i = 1; i < angle.GetLength(0) - 1; i++)
            {
                for (int j = 1; j < angle.GetLength(1) - 1; j++)
                {
                    //linear interpolation
                    if (angle[i, j] > -45)
                    {
                        if (angle[i, j] > 0)
                        {
                            if (angle[i, j] > 45)
                            {
                                if (angle[i, j] > 90)
                                {
                                    if (angle[i, j] > 135)
                                    {
                                        //135 to 180 1
                                        double l = (angle[i, j] - 135) / 45;
                                        double r = 1 - l;
                                        double a = Math.Max(gIntensity[i - 1, j - 1] * l + gIntensity[i - 1, j] * r, gIntensity[i + 1, j + 1] * l + gIntensity[i + 1, j] * r);
                                        if (gIntensity[i, j] >= a)
                                        {
                                            x[i, j] = gIntensity[i, j];
                                        }
                                    }
                                    else
                                    {
                                        //90 to 135 1
                                        //main diagonal
                                        double l = (angle[i, j] - 90) / 45;
                                        double r = 1 - l;
                                        double a = Math.Max(gIntensity[i, j + 1] * l + gIntensity[i + 1, j + 1] * r, gIntensity[i, j - 1] * l + gIntensity[i - 1, j - 1] * r);
                                        if (gIntensity[i, j] >= a)
                                        {
                                            x[i, j] = gIntensity[i, j];
                                        }
                                    }
                                }
                                else
                                {
                                    //45 to 90 1
                                    //side diagonal
                                    double l = (angle[i, j] - 45) / 45;
                                    double r = 1 - l;
                                    double a = Math.Max(gIntensity[i - 1, j + 1] * l + gIntensity[i, j + 1] * r, gIntensity[i + 1, j - 1] * l + gIntensity[i, j - 1] * r);
                                    if (gIntensity[i, j] >= a)
                                    {
                                        x[i, j] = gIntensity[i, j];
                                    }
                                }
                            }
                            else
                            {
                                //0 to 45 1
                                //vertical right
                                double l = angle[i, j] / 45;
                                double r = 1 - l;
                                double a = Math.Max(gIntensity[i - 1, j] * l + gIntensity[i - 1, j + 1] * r, gIntensity[i + 1, j] * l + gIntensity[i + 1, j - 1] * r);
                                if (gIntensity[i, j] >= a)
                                {
                                    x[i, j] = gIntensity[i, j];
                                }
                            }
                        }
                        else
                        {
                            //-45 to 0 1
                            //vertical left
                            double l = (-45 - (angle[i, j])) / -45;
                            double r = 1 - l;
                            double a = Math.Max(gIntensity[i - 1, j - 1] * l + gIntensity[i - 1, j] * r, gIntensity[i + 1, j + 1] * l + gIntensity[i + 1, j] * r);
                            if (gIntensity[i, j] >= a)
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                        }
                    }
                    else
                    {
                        if (angle[i, j] < -90)
                        {
                            if (angle[i, j] < -135)
                            {
                                //-180 to -135 1
                                double l = (angle[i, j] + 135) / -45;
                                double r = 1 - l;
                                double a = Math.Max(gIntensity[i - 1, j] * r + gIntensity[i - 1, j + 1] * l, gIntensity[i + 1, j] * r + gIntensity[i + 1, j - 1] * l);
                                if (gIntensity[i, j] >= a)
                                {
                                    x[i, j] = gIntensity[i, j];
                                }
                            }
                            else
                            {
                                //-135 to -90 1
                                double l = (angle[i, j] + 90) / -45;
                                double r = 1 - l;
                                double a = Math.Max(gIntensity[i, j + 1] * l + gIntensity[i - 1, j + 1] * r, gIntensity[i, j - 1] * l + gIntensity[i + 1, j - 1] * r);
                                if (gIntensity[i, j] >= a)
                                {
                                    x[i, j] = gIntensity[i, j];
                                }
                            }
                        }
                        else
                        {
                            //-90 to -45
                            double l = (angle[i, j] + 45) / -45;
                            double r = 1 - l;
                            double a = Math.Max(gIntensity[i - 1, j - 1] * l + gIntensity[i, j - 1] * r, gIntensity[i + 1, j + 1] * l + gIntensity[i, j + 1] * r);
                            if (gIntensity[i, j] >= a)
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                        }
                    }
                }
            }
            gIntensity = x;
        }
        Bitmap ReconstructFromGradient()
        {
            if (label != null)
            {
                Bitmap res = new Bitmap(label.GetLength(1), label.GetLength(0), System.Drawing.Imaging.PixelFormat.Format32bppArgb);
                for (int i = 0; i < label.GetLength(0); i++)
                {
                    for (int j = 0; j < label.GetLength(1); j++)
                    {
                        if (label[i, j] % 10 > 0)
                        {
                            res.SetPixel(j, i, Color.White);
                        }
                        else
                        {
                            res.SetPixel(j, i, Color.Black);
                        }
                    }
                }
                return res;
            }
            else
            {
                Bitmap res = new Bitmap(gIntensity.GetLength(1), gIntensity.GetLength(0), System.Drawing.Imaging.PixelFormat.Format32bppArgb);
                //divide gradient by max value, then multiply by 255
                //max = gIntensity.Cast<double>().Max();
                for (int i = 0; i < gIntensity.GetLength(0); i++)
                {
                    for (int j = 0; j < gIntensity.GetLength(1); j++)
                    {
                        int intensity = (int)((gIntensity[i, j] / max) * 255);
                        res.SetPixel(j, i, Color.FromArgb(intensity, intensity, intensity));
                    }
                }
                ImageHandler.factory.Load(res);
                return res;
            }
        }
        public Image EdgeEffect(Image x)
        {
            int[,] sobelX = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
            int[,] sobelY = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
            actual = x; 
            GaussianFilter();
            Bitmap temp = new Bitmap(actual);
            int[,] Gx = Convolve(ref temp, ref sobelX);
            int[,] Gy = Convolve(ref temp, ref sobelY);
            ComputeGradient(ref Gx, ref Gy);
            Bitmap res = new Bitmap(gIntensity.GetLength(1), gIntensity.GetLength(0), System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            double max = gIntensity.Cast<double>().Max();
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    Color t = temp.GetPixel(j, i);
                    double ins = (gIntensity[i, j] / max) * 3;
                    int red = Math.Clamp((int)(t.R * ins), 0, 255);
                    int green = Math.Clamp((int)(t.G * ins), 0, 255);
                    int blue = Math.Clamp((int)(t.B * ins), 0, 255);
                    res.SetPixel(j, i, Color.FromArgb(red, green, blue));
                }
            }
            return res;
        }
        //4. Apply double treshold to determine potential edges
        void DoubleThreshold()
        {
            label = new byte[gIntensity.GetLength(0), gIntensity.GetLength(1)];
            //max = gIntensity.Cast<double>().Max();
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    if (gIntensity[i, j] < lowerT * max)
                    {
                        //weak edges, which are disgarded
                        gIntensity[i, j] = 0;
                        label[i, j] = 0;
                    }
                    else if (gIntensity[i, j] < upperT * max)
                    {
                        //edges for which we use hysteresis
                        label[i, j] = 1;
                    }
                    else
                    {
                        //strong edges, which are always kept
                        label[i, j] = 2;
                    }
                }
            }
        }
        //5. Track edge by Hysteresis. Finalize detection by suppressing all the other edges that are weak and not connected to strong edges

        //idea: split image into horizontal bars for parallelization
        //for each bar, go through the pixels
        //if they aren't 0, and aren't 'occupied', bfs to find all connected pixels(in the current strip), and note if they touch the upper or lower border
        //then, for all objects that contain strong pixels, preserve them, for those than dont, if they touch the border, and through it an object which contains a strong pixel, preserve them
        //delete all the rest
        //boom, parallelized
        //in CPU implementation, the objects can immediately be preserved or destroyed, as there are no bars
        void Hysteresis() //definitely not it
        {
            for (int i = 0; i < label.GetLength(0); i++)
            {
                for (int j = 0; j < label.GetLength(1); j++)
                {
                    if (label[i, j] < 10)
                    {
                        //if label[i,j] is >=10, it means it was visited
                        if (label[i, j] == 0)
                        {
                            label[i, j] += 10;
                        }
                        else
                        {
                            //constructor of h is doing BFS
                            HelperObject h = new HelperObject(ref label, i, j);
                            if (!h.hasStrong || h.body.Count == 1)
                            {
                                h.MutuallyAssuredDestruction(ref label);
                            }
                        }
                    }
                }
            }
        }

        //final - apply 
        public Image Apply(ref Image x)
        {
            if (x == null)
            {
                return null;
            }
            actual = x;
            ToGrayscale();
            GaussianFilter();
            Gradient();
            NonMaxSuppression();
            DoubleThreshold();
            Hysteresis();
            actual = ReconstructFromGradient();
            return actual;
        }
        public bool[,] matrixDirect(ref Image x)
        {
            if (x == null)
            {
                return null;
            }
            actual = x;
            ToGrayscale();
            GaussianFilter();
            Gradient();
            NonMaxSuppression();
            DoubleThreshold();
            Hysteresis();
            bool[,] ret = new bool[label.GetLength(0), label.GetLength(1)];
            for (int i = 0; i < label.GetLength(0); i++)
            {
                for (int j = 0; j < label.GetLength(1); j++)
                {
                    ret[i, j] = label[i, j] % 10 > 0;
                }
            }
            return ret;
        }
    }
}
