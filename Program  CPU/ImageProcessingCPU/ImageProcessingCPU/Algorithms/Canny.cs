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
        double[,] angle;
        double[,] gIntensity;
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
            Bitmap Gx = new Bitmap(ImageHandler.factory.Image);
            Bitmap Gy = new Bitmap(ImageHandler.factory.Image);
            Gx = Convolve(Gx, ref sobelX);
            Gy = Convolve(Gy, ref sobelY);
            //computing gradient intensity via G = sqrt(Gx^2+Gy^2)
            ComputeGradient(Gx, Gy);

            //computing gradient angle - needed for tresholding
            angle = GradientDirection(Gx, Gy);
            //non-max suppresion
            NonMaxSuppression();
            actual = ReconstructFromGradient();
            //actual = temp;
            ImageHandler.factory.Load(actual);
        }
        //also implement channels
        void ComputeGradient(Bitmap Gx, Bitmap Gy)
        {
            gIntensity = new double[Gx.Height, Gx.Width];
            for (int i = 0; i < gIntensity.GetLength(0); i++)
            {
                for (int j = 0; j < gIntensity.GetLength(1); j++)
                {
                    Color tX = Gx.GetPixel(j, i);
                    Color tY = Gy.GetPixel(j, i);
                    gIntensity[i, j] = Math.Sqrt(tX.R * tX.R + tY.R * tY.R);
                }
            }
        }
        //this assumes R=G=B, or can work for a specified channel, should write an overload
        double[,] GradientDirection(Bitmap Gx, Bitmap Gy)
        {
            double[,] res = new double[Gx.Height, Gx.Width];
            for (int i = 0; i < res.GetLength(0); i++)
            {
                for (int j = 0; j < res.GetLength(1); j++)
                {
                    //calculating angle as arctan(Gx,Gy) and converting to degrees by
                    double rad = Math.Atan2(Gx.GetPixel(j, i).R, Gy.GetPixel(j, i).R) * (180 / Math.PI);
                    int direction;
                    if ((rad >= -22.5 && rad < 22.5) || rad > 157.5 || rad <= -157.5)
                    {
                        direction = 0; //vertical |
                    }
                    else if ((rad >= 22.5 && rad < 67.5) || (rad > -157.5 && rad <= -112.5))
                    {
                        direction = 1; //diagonal 1 /
                    }
                    else if ((rad >= 67.5 && rad < 112.5) || (rad > -112.5 && rad <= -67.5))
                    {
                        direction = 2; //horizontal --
                    }
                    else
                    {
                        direction = 3; //diagonal 2 \
                    }
                    res[i, j] = direction;
                }
            }
            return res;
        }
        Bitmap Convolve(Bitmap target, ref int[,] filter)
        {
            Bitmap ret = new Bitmap(target);
            for (int i = 0; i < target.Height; i++)
            {
                for (int j = 0; j < target.Width; j++)
                {
                    //here
                    int sumR = 0;
                    int sumG = 0;
                    int sumB = 0;
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
                            sumG += t.G * filter[a, b];
                            sumB += t.B * filter[a, b];
                        }
                    }
                    ret.SetPixel(j, i, Color.FromArgb(byte.MaxValue, (byte)Math.Clamp(sumR, 0, 255), (byte)Math.Clamp(sumG, 0, 255), (byte)Math.Clamp(sumB, 0, 255)));
                }
            }
            return ret;
        }
        //3. Apply gradient magnitude tresholding or lower-bound cutoff suppression to get rid of spurious response to edge detection
        //non-maximum suppresion
        //also works only with greyscale. consider overloading for color channels
        void NonMaxSuppression()
        {
            //edges
            //upper
            for (int i = 1; i < angle.GetLength(1) - 1; i++)
            {
                if (angle[0, i] == 0)
                {
                    double p1 = Math.Max(gIntensity[0, i - 1], gIntensity[0, i + 1]);
                    double p2 = gIntensity[0, i];
                    if (p2 != Math.Max(p1, p2))
                    {
                        gIntensity[0, i] = 0;
                    }
                }
            }
            //lower
            for (int i = 1; i < angle.GetLength(1) - 1; i++)
            {
                if (angle[angle.GetLength(0) - 1, i] == 0)
                {
                    double p1 = Math.Max(gIntensity[angle.GetLength(0) - 1, i - 1], gIntensity[angle.GetLength(0) - 1, i + 1]);
                    double p2 = gIntensity[angle.GetLength(0) - 1, i];
                    if (p2 != Math.Max(p1, p2))
                    {
                        gIntensity[angle.GetLength(0) - 1, i] = 0;
                    }
                }
            }
            //left
            for (int i = 1; i < angle.GetLength(0) - 1; i++)
            {
                if (angle[i, 0] == 2)
                {
                    double p1 = Math.Max(gIntensity[i - 1, 0], gIntensity[i + 1, 0]);
                    double p2 = gIntensity[i, 0];
                    if (p2 < p1)
                    {
                        gIntensity[i, 0] = 0;
                    }
                }
            }
            //right
            for (int i = 1; i < angle.GetLength(0) - 1; i++)
            {
                if (angle[i, angle.GetLength(1) - 1] == 2)
                {
                    double p1 = Math.Max(gIntensity[i - 1, angle.GetLength(1) - 1], gIntensity[i + 1, angle.GetLength(1) - 1]);
                    double p2 = gIntensity[i, angle.GetLength(1) - 1];
                    if (p2 < p1)
                    {
                        gIntensity[i, angle.GetLength(1) - 1] = 0;
                    }
                }
            }
            //inner image
            for (int i = 1; i < angle.GetLength(0) - 1; i++)
            {
                for (int j = 1; j < angle.GetLength(1) - 1; j++)
                {
                    double p1;
                    double p2;
                    switch (angle[i, j])
                    {
                        case 0:
                            //vertical
                            p1 = Math.Max(gIntensity[i - 1, j], gIntensity[i + 1, j]);
                            p2 = gIntensity[i, j];
                            if (p2 < p1)
                            {
                                gIntensity[i, j] = 0;
                            }
                            break;
                        case 1:
                            //diagonal 1
                            p1 = Math.Max(gIntensity[i - 1, j + 1], gIntensity[i + 1, j - 1]);
                            p2 = gIntensity[i, j];
                            if (p2 < p1)
                            {
                                gIntensity[i, j] = 0;
                            }
                            break;
                        case 2:
                            //horizontal 
                            p1 = Math.Max(gIntensity[i, j - 1], gIntensity[i, j + 1]);
                            p2 = gIntensity[i, j];
                            if (p2 < p1)
                            {
                                gIntensity[i, j] = 0;
                            }
                            break;
                        case 4:
                            //diagonal 2
                            p1 = Math.Max(gIntensity[i - 1, j - 1], gIntensity[i + 1, j + 1]);
                            p2 = gIntensity[i, j];
                            if (p2 < p1)
                            {
                                gIntensity[i, j] = 0;
                            }
                            break;
                    }
                }
            }
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
