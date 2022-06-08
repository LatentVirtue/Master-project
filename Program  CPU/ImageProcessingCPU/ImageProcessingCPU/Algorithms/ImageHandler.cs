using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Windows.Forms;
using ImageProcessor;
using ImageProcessingCPU.Algorithms;

namespace ImageProcessingCPU
{
    class ImageHandler
    {
        public static int HoughNLines = 15;
        static bool canny = false;
        bool hough;
        static Image original;
        static Image temporary;
        public static Image Temporary
        {
            get
            {
                return temporary;
            }
            set
            {
                temporary = value;
            }
        }
        static Image preview;
        public static Image Preview
        {
            get
            {
                return preview;
            }
        }
        public static ImageFactory factory = new ImageFactory();
        public static bool Load()
        {
            using OpenFileDialog dialog = new OpenFileDialog();
            dialog.Filter = "Image Files(*.BMP; *.JPG; *.PNG)| *.BMP; *.JPG; *.PNG | All files(*.*) | *.*";
            dialog.RestoreDirectory = true;
            if (dialog.ShowDialog() == DialogResult.OK)
            {
                original = Image.FromFile(dialog.FileName);
                temporary = original;
                factory.Load(temporary);
                return true;
            }
            return false;
        }
        public static bool Save()
        {
            using SaveFileDialog dialog = new SaveFileDialog();
            dialog.Filter = "Bitmap (*.BMP)| *.BMP | JPG (*.JPG)| *.JPG";
            if (dialog.ShowDialog() == DialogResult.OK)
            {
                factory.Load(temporary);
                factory.Save(dialog.FileName);
                return true;
            }
            return false;
        }
        public static void Refresh()
        {
            if (!Check())
            {
                return;
            }
            factory.Load(temporary);
            factory.Resize(new Size(700, 400));
            preview = factory.Image;
            canny = false;
        }
        public static bool Check()
        {
            if (original == null)
            {
                if (Load())
                {
                    if (temporary == null)
                    {
                        temporary = original;
                    }
                    return true;
                }
                else
                {
                    return false;
                }
            }
            return true;
        }
        public static Image Rotate90()
        {
            if (!Check())
            {
                return null;
            }
            TestAlgo x = new TestAlgo();
            temporary = x.Apply(ref temporary);
            Refresh();
            return preview;
        }
        public static Image Canny_test(int c, double lT, double uT)
        {
            if (!Check())
            {
                return null;
            }
            using Canny x = new Canny(c, lT, uT);
            temporary = x.Apply(original);
            Refresh();
            canny = true;
            return preview;
        }
        public static Image Edge_effect(int c, bool blur)
        {
            if (!Check())
            {
                return null;
            }
            using Canny x = new Canny(c, 0.1, 0.3);
            temporary = x.EdgeEffect(original, blur);
            Refresh();
            return preview;
        }
        public static Image Hough()
        {
            if (!Check())
            {
                return null;
            }

            if (canny)
            {
                factory.Load(temporary);

            }
            else
            {
                factory.Load(original);
            }
            if (original.Width != original.Height)
            {
                int n = factory.Image.Height < factory.Image.Width ? factory.Image.Height : factory.Image.Width;
                bool u = factory.Image.Height < factory.Image.Width;
                int sx = 0;
                int sy = 0;
                if (u)
                {
                    sx = n / 4;
                }
                else
                {
                    sy = n / 4;
                }
                temporary = factory.Crop(new Rectangle(sx, sy, n, n)).Image;
            }
            using HoughTransform x = new HoughTransform(temporary, canny, 3, 0.5, 0.01, 2, HoughNLines);
            temporary = x.Apply(temporary);
            Refresh();
            return preview;
        }
        public static Image Convolve()
        {
            if (!Check())
            {
                return null;
            }
            Bitmap target = new Bitmap(temporary);
            Bitmap res = new Bitmap(target.Width,target.Height);
            //int[,] filter = { { -2, -1, 0 }, { -1, 1, 1 }, { 0, 1, 2 } }; //emboss
            //int[,] filter = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } }; //nabla fy
            int[,] filter = { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } }; //blur
            int fsum = 0;
            for (int i = 0; i < filter.GetLength(0); i++)
            {
                for (int j = 0; j < filter.GetLength(1); j++)
                {
                    fsum += filter[i, j];
                }
            }
            if (fsum == 0)
            {
                fsum = 1;
            }
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
                            sumR += (t.R * filter[a, b]) / fsum;
                            sumG += (t.G * filter[a, b]) / fsum;
                            sumB += (t.B * filter[a, b]) / fsum;
                        }
                    }
                    res.SetPixel(j, i, Color.FromArgb(255, Math.Clamp(sumR, 0, 255), Math.Clamp(sumG, 0, 255), Math.Clamp(sumB, 0, 255)));
                }
            }
            temporary = res;
            Refresh();
            return preview;
        }
    }
}
