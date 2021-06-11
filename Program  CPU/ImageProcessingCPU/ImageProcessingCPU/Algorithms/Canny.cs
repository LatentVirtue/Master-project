using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using ImageProcessor;
using ImageProcessor.Imaging.Filters.Photo;
using ImageProcessor.Imaging;
using ImageProcessingCPU;
using System.Linq;

namespace ImageProcessingCPU.Algorithms
{
    class Canny : IAlgoInterface
    {
        Image actual;
        int[,] angle;
        double[,] gIntensity;
        double adj = 1;
        public Canny()
        {

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
            //ImageHandler.factory.Load(list[1]);
            ImageHandler.factory.GaussianBlur(5);
            actual = ImageHandler.factory.Image;
        }
        //2. Intensity gradient
        void Gradient()
        {
            //Sobel operator kernels
            int[,] sobelX = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
            int[,] sobelY = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
            //applying kernels
            Bitmap temp = new Bitmap(actual);
            int[,] Gx = Convolve(ref temp, ref sobelX);
            int[,] Gy = Convolve(ref temp, ref sobelY);
            //computing gradient intensity via G = sqrt(Gx^2+Gy^2)
            ComputeGradient(ref Gx, ref Gy);

            //computing gradient angle - needed for tresholding
            GradientDirection(ref Gx, ref Gy);
            //non-max suppresion
            NonMaxSuppression();
            actual = ReconstructFromGradient();
            //actual = temp;
            ImageHandler.factory.Load(actual);
        }
        //also implement channels
        void ComputeGradient(ref int[,] Gx, ref int[,] Gy)
        {
            gIntensity = new double[Gx.GetLength(0), Gx.GetLength(1)];
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    gIntensity[i, j] = Math.Sqrt(Gx[i, j] * Gx[i, j] + Gy[i, j] * Gy[i, j]);
                }
            }
        }
        //this assumes R=G=B, or can work for a specified channel, should write an overload
        void GradientDirection(ref int[,] Gx, ref int[,] Gy)
        {
            angle = new int[Gx.GetLength(0), Gx.GetLength(1)];
            for (int i = 0; i < angle.GetLength(0); i++)
            {
                for (int j = 0; j < angle.GetLength(1); j++)
                {
                    //calculating angle as arctan(Gx,Gy) and converting to degrees by
                    double rad = (Math.Atan2(Gx[i, j], Gy[i, j]) * (180 / Math.PI)) % 180;
                    int direction;
                    //vertical |
                    //diagonal 1 /
                    //horizontal --
                    //diagonal 2 \
                    if (rad > 0)
                    {
                        if (rad < 22.5)
                        {
                            direction = 0;
                        }
                        else
                        {
                            if (rad < 67.5)
                            {
                                direction = 1;
                            }
                            else
                            {
                                if (rad < 112.5)
                                {
                                    direction = 2;
                                }
                                else
                                {
                                    if (rad < 157.5)
                                    {
                                        direction = 3;
                                    }
                                    else
                                    {
                                        direction = 0;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if (rad > -22.5)
                        {
                            direction = 0;
                        }
                        else
                        {
                            if (rad > -67.5)
                            {
                                direction = 3;
                            }
                            else
                            {
                                if (rad > -112.5)
                                {
                                    direction = 2;
                                }
                                else
                                {
                                    if (rad > -157.5)
                                    {
                                        direction = 1;
                                    }
                                    else
                                    {
                                        direction = 0;
                                    }
                                }
                            }
                        }
                    }
                    angle[i, j] = direction;
                }
            }
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
        //nonMaxSuppression fails on diagonal 2. Needs rework or better checks.
        void NonMaxSuppression()
        {
            //inner image
            double[,] x = new double[gIntensity.GetLength(0), gIntensity.GetLength(1)];
            for (int i = 1; i < angle.GetLength(0) - 1; i++)
            {
                for (int j = 1; j < angle.GetLength(1) - 1; j++)
                {
                    switch (angle[i, j])
                    {
                        case 0:
                            //vertical angle = horizontal edge
                            if (gIntensity[i, j] < Math.Max(gIntensity[i + 1, j], gIntensity[i - 1, j]) * adj)
                            {
                                x[i, j] = 0;
                            }
                            else
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                            break;
                        case 1:
                            //diagonal 1
                            if (gIntensity[i, j] < Math.Max(gIntensity[i - 1, j + 1], gIntensity[i + 1, j - 1]) * adj)
                            {
                                x[i, j] = 0;
                            }
                            else
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                            break;
                        case 2:
                            //horizontal 
                            if (gIntensity[i, j] < Math.Max(gIntensity[i, j + 1], gIntensity[i, j - 1]) * adj)
                            {
                                x[i, j] = 0;
                            }
                            else
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                            break;
                        case 3:
                            //diagonal 2
                            if (gIntensity[i, j] < Math.Max(gIntensity[i + 1, j + 1], gIntensity[i - 1, j - 1]) * adj)
                            {
                                x[i, j] = 0;
                            }
                            else
                            {
                                x[i, j] = gIntensity[i, j];
                            }
                            break;
                    }
                }
            }
            gIntensity = x;
        }
        Bitmap ReconstructFromGradient()
        {
            Bitmap res = new Bitmap(gIntensity.GetLength(1), gIntensity.GetLength(0), System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            //divide gradient by max value, then multiply by 255
            double max = gIntensity.Cast<double>().Max();
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    int intensity = (int)((gIntensity[i, j] / max) * 255);
                    res.SetPixel(j, i, Color.FromArgb(intensity, intensity, intensity));
                }
            }
            return res;
        }
        //4. Apply double treshold to determine potential edges
        //5. Track edge by Hysteresis. Finalize detection by suppressing all the other edges that are weak and not connected to strong edges


        //final - apply 
        public Image Apply(Image x)
        {
            if (x == null)
            {
                return null;
            }
            actual = x;
            ToGrayscale();
            GaussianFilter();
            Gradient();
            return actual;
        }
    }
}
