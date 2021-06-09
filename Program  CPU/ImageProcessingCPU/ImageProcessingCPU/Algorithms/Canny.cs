using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using ImageProcessor;
using ImageProcessor.Imaging.Filters.Photo;
using ImageProcessor.Imaging;
using ImageProcessingCPU;

namespace ImageProcessingCPU.Algorithms
{
    class Canny : IAlgoInterface
    {
        Image actual;
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
            int[,] sobelX = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
            int[,] sobelY = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
            Bitmap temp = new Bitmap(ImageHandler.factory.Image);
            temp = Convolve(temp, ref sobelX);
            temp = Convolve(temp, ref sobelY);
            actual = temp;
            ImageHandler.factory.Load(actual);
            //ImageHandler.factory.Filter();
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
                            if (ci < 0 || ci >= target.Height || cj < 0 || cj >= target.Width)
                            {
                                continue;
                                Color t = target.GetPixel(j, i);
                                sumR += Math.Clamp(t.R * filter[a, b], 0, 255);
                                sumG += Math.Clamp(t.G * filter[a, b], 0, 255);
                                sumB += Math.Clamp(t.B * filter[a, b], 0, 255);
                            }
                            else
                            {
                                Color t = target.GetPixel(cj, ci);
                                sumR += t.R * filter[a, b];
                                sumG += t.G * filter[a, b];
                                sumB += t.B * filter[a, b];
                            }
                        }
                    }
                    ret.SetPixel(j, i, Color.FromArgb(byte.MaxValue, (byte)Math.Clamp(sumR, 0, 255), (byte)Math.Clamp(sumG, 0, 255), (byte)Math.Clamp(sumB, 0, 255)));
                }
            }
            return ret;
        }
        //3. Apply gradient magnitude tresholding or lower-bound cutoff suppression to get rid of spurious response to edge detection

        void Tresholding(ref Bitmap target)
        {
            
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
            //ToGrayscale();
            GaussianFilter();
            Gradient();
            return actual;
        }
    }
}
