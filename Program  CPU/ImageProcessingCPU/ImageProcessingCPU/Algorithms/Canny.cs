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
    class Canny : IDisposable
    {
        int selectedKernel;
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
        Bitmap actual;
        double[,] angle;
        double[,] gIntensity;
        double[,] grayscale;
        byte[,] label;
        double max;
        double lowerT = 0.1;
        double upperT = 0.3;
        public Canny(int kOperator, double lT, double uT)
        {
            selectedKernel = kOperator;
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
            lowerT = lT;
            upperT = uT;
        }
        public void Dispose()
        {
            actual = null;
            angle = null;
            gIntensity = null;
            grayscale = null;
            label = null;
        }
        //here I picked CPU because of one-pass nature of the function. GetPixel is not very efficient, but that isn't solved with GPU implementation
        //after grayscale has been applied, everything is passed to CannyGPU
        void ToGrayscale()
        {
            grayscale = new double[actual.Height, actual.Width];
            for (int i = 0; i < grayscale.GetLength(0); i++)
            {
                for (int j = 0; j < grayscale.GetLength(1); j++)
                {
                    Color t = actual.GetPixel(j, i);
                    grayscale[i, j] = (t.R + t.G + t.B) / 3.0;
                }
            }
        }
        //in order to not pass data back and forth between GPU and CPU, the matrixes used the first few steps of canny will be kept in the GPU buffer until completion
        //1. Gaussian filter
        //originally used CPU, switched to GPU as its a convolution.
        void GaussianFilter()
        {
            //normalized gaussian kernel for convolution
            double[,] gaussianMatrix = {    {0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067 },
                                            {0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
                                            {0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
                                            {0.00038771, 0.01330373, 0.11098164, 0.22508352, 0.11098164, 0.01330373, 0.00038771 },
                                            {0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
                                            {0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
                                            {0.00000067, 0.00002292, 0.00019117, 0.00038771 ,0.00019117, 0.00002292, 0.00000067 }};
            double[,] notNormalizedGaussMatrix = { { 1, 4, 7, 4, 1 }, { 4, 16, 26, 16, 4 }, { 7, 26, 41, 26, 7 }, { 4, 16, 26, 16, 4 }, { 1, 4, 7, 4, 1 } };
            double[,] smallGaussianMatrix = { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } };
            CannyGPU.ApplyDouble(gaussianMatrix, ref grayscale);
            //Convolution.ApplyDouble(notNormalizedGaussMatrix, ref grayscale);
        }
        //2. Intensity gradient
        //this is also done on GPU because it is a simple convolution
        void GradientGPU()
        {
            double[,] Gx = grayscale;
            double[,] Gy = grayscale;
            //convert kernels to double[,]
            int n = opX.GetLength(0);
            int m = opX.GetLength(1);
            double[,] kX = new double[n, m];
            double[,] kY = new double[n, m];
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < m; j++)
                {
                    kX[i, j] = opX[i, j];
                    kY[i, j] = opY[i, j];
                }
            }
            //convolution with converted kernels
            CannyGPU.ApplyDouble(kX, ref Gx);
            CannyGPU.ApplyDouble(kY, ref Gy);
            //computing gradient intensity via G = sqrt(Gx^2+Gy^2)
            ComputeGradient(ref Gx, ref Gy);
        }
        void Gradient()
        {
            //applying kernels
            Bitmap temp = new Bitmap(actual);
            int[,] Gx = Convolve(ref temp, opX);
            int[,] Gy = Convolve(ref temp, opY);
            //computing gradient intensity via G = sqrt(Gx^2+Gy^2)
            ComputeGradient(ref Gx, ref Gy);
            //actual = temp;
            ImageHandler.factory.Load(actual);
        }
        //also implement channels
        //staying on CPU due to one-pass with it is faster than copying it to GPU, then passing back the computed thing
        void ComputeGradient(ref int[,] Gx, ref int[,] Gy)
        {
            gIntensity = new double[Gx.GetLength(0), Gx.GetLength(1)];
            angle = new double[Gx.GetLength(0), Gx.GetLength(1)];
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    gIntensity[i, j] = Math.Sqrt(Gx[i, j] * Gx[i, j] + Gy[i, j] * Gy[i, j]);
                    if (max < gIntensity[i, j])
                    {
                        max = gIntensity[i, j];
                    }
                    angle[i, j] = (Math.Atan2(Gx[i, j], Gy[i, j]) * (180 / Math.PI)) % 180;
                }
            }
        }
        void ComputeGradient(ref double[,] Gx, ref double[,] Gy)
        {
            gIntensity = new double[Gx.GetLength(0), Gx.GetLength(1)];
            angle = new double[Gx.GetLength(0), Gx.GetLength(1)];
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    gIntensity[i, j] = Math.Sqrt(Gx[i, j] * Gx[i, j] + Gy[i, j] * Gy[i, j]);
                    if (max < gIntensity[i, j])
                    {
                        max = gIntensity[i, j];
                    }
                    angle[i, j] = (Math.Atan2(Gx[i, j], Gy[i, j]) * (180 / Math.PI)) % 180;
                }
            }
        }
        //also add an overload for color channels
        int[,] Convolve(ref Bitmap target, int[,] filter)
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
                            
                            //test coloring
                            /*
                            if(label[i,j] <20 && label[i,j]%10 > 1)
                            {
                                res.SetPixel(j, i, Color.White);
                            }
                            else if(label[i,j] < 20)
                            {
                                res.SetPixel(j, i, Color.Green);
                            }
                            else
                            {
                                res.SetPixel(j, i, Color.Red);
                            }
                            */
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
        public Image EdgeEffect(Image x, bool blur)
        {
            Bitmap temp = (Bitmap)x;
            Bitmap res = new Bitmap(x.Width, x.Height);
            CannyGPU.TroupleDouble[,] target = new CannyGPU.TroupleDouble[x.Height, x.Width];
            for (int i = 0; i < x.Height; i++)
            {
                for (int j = 0; j < x.Width; j++)
                {
                    target[i, j] = new CannyGPU.TroupleDouble(temp.GetPixel(j, i));
                }
            }
            double[,] gradientHelp = CannyGPU.edgefix(selectedKernel, target, blur);
            double max = gradientHelp.Cast<double>().Max();
            for (int i = 0; i < gradientHelp.GetLength(0); i++)
            {
                for (int j = 0; j < gradientHelp.GetLength(1); j++)
                {
                    Color t = temp.GetPixel(j, i);
                    double ins = (gradientHelp[i, j] / max) * 3;
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
                        //weak edges, which are discarded
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

        //this has been left as a CPU-only operation, as many memory dumps to and from the GPU will completely negate the time gains from kernel exectution

        //in CPU implementation, the objects can immediately be preserved or destroyed, as there are no bars
        void Hysteresis()
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
                                /*
                                foreach (Point p in h.body)
                                {
                                    label[p.X, p.Y] += 20;
                                }
                                */
                            }
                        }
                    }
                }
            }
        }

        //final - apply 
        public Image Apply(Image x)
        {
            if (x == null)
            {
                return null;
            }
            actual = (Bitmap)x;
            //ToGrayscale();
            //GaussianFilter();
            //GradientGPU();
            //NonMaxSuppression();
            //DoubleThreshold();
            CannyGPU.TroupleDouble[,] target = new CannyGPU.TroupleDouble[actual.Height, actual.Width];
            for(int i = 0; i < actual.Height; i++)
            {
                for(int j = 0; j < actual.Width; j++)
                {
                    target[i, j] = new CannyGPU.TroupleDouble(actual.GetPixel(j, i));
                }
            }
            label = CannyGPU.ApplyAll(selectedKernel, target, lowerT, upperT);
            Hysteresis();
            actual = ReconstructFromGradient();
            return actual;
        }
        public bool[,] matrixDirect(ref Bitmap x, bool u)
        {
            if (x == null)
            {
                return null;
            }
            if(u)
            {
                bool[,] t = new bool[x.Height, x.Width];
                for (int i = 0; i < x.Height; i++)
                {
                    for (int j = 0; j < x.Width; j++)
                    {
                        Color temp = x.GetPixel(j, i);
                        t[i, j] = temp.R > 0 || temp.G > 0 || temp.B > 0;
                    }
                }
                return t;
            }
            //ToGrayscale();
            //GaussianFilter();
            //GradientGPU();
            //NonMaxSuppression();
            //DoubleThreshold();
            CannyGPU.TroupleDouble[,] target = new CannyGPU.TroupleDouble[x.Height, x.Width];
            for (int i = 0; i < x.Height; i++)
            {
                for (int j = 0; j < x.Width; j++)
                {
                    target[i, j] = new CannyGPU.TroupleDouble(x.GetPixel(j, i));
                }
            }
            label = CannyGPU.ApplyAll(selectedKernel, target, lowerT, upperT);
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
