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
            dialog.Filter = "Bitmap (*.BMP)| *.BMP";
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
            temporary = x.Apply(temporary);
            Refresh();
            return preview;
        }

        public static Image Canny_test(int c)
        {
            if (!Check())
            {
                return null;
            }
            Canny x = new Canny(c);
            temporary = x.Apply(temporary);
            Refresh();
            return preview;
        }
        public static Image Edge_effect(int c)
        {
            if (!Check())
            {
                return null;
            }
            Canny x = new Canny(c);
            temporary = x.EdgeEffect(temporary);
            Refresh();
            return preview;
        }
    }
}
