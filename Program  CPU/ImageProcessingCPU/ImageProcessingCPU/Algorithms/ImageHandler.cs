using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Windows.Forms;
using ImageProcessor;

namespace ImageProcessingCPU
{
    class ImageHandler
    {
        static Image original;
        static Image temporary;
        static Image preview;
        static ImageFactory factory = new ImageFactory();
        public static bool Load()
        {
            using (OpenFileDialog dialog = new OpenFileDialog())
            {
                dialog.Filter = "Image Files(*.BMP; *.JPG; *.PNG)| *.BMP; *.JPG; *.PNG | All files(*.*) | *.*";
                dialog.RestoreDirectory = true;
                if(dialog.ShowDialog() == DialogResult.OK)
                {
                    original = Image.FromFile(dialog.FileName);
                    return true;
                }
                return false;
            }
        }
        public static Image Refresh()
        {
            if(original == null)
            {
                Load();
                factory.Load(original);
                factory.Resize(new Size(700, 400));
                preview = factory.Image;
            }
            else if(temporary == null)
            {
                factory.Load(temporary);
                factory.Resize(new Size(700, 400));
                preview = factory.Image;
            }
            return preview;
        }
    }
}
