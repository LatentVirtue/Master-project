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
            temporary = x.Apply(temporary);
            Refresh();
            return preview;
        }
        public static Image Edge_effect(int c, bool blur)
        {
            if (!Check())
            {
                return null;
            }
            using Canny x = new Canny(c, 0.1, 0.3);
            temporary = x.EdgeEffect(temporary, blur);
            Refresh();
            return preview;
        }
        public static Image Hough()
        {
            if (!Check())
            {
                return null;
            }
            using HoughTransform x = new HoughTransform(ref temporary,2,0.5,0.01,2);
            x.Apply(ref temporary);
            Refresh();
            return preview;
        }
    }
}
