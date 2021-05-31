using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Windows.Forms;

namespace ImageProcessingCPU.Algorithms
{
    class ImageHandler
    {
        Image original;
        public Image temporary;
        Image preview;
        public bool Load()
        {
            string filepath = "";
            using (OpenFileDialog dialog = new OpenFileDialog())
            {
                dialog.Filter = "Image Files(*.BMP; *.JPG; *.PNG)| *.BMP; *.JPG; *.PNG | All files(*.*) | *.*";
                dialog.RestoreDirectory = true;
                if(dialog.ShowDialog() == DialogResult.OK)
                {
                    filepath = dialog.FileName;
                    original = Bitmap.FromFile(filepath);
                    return true;
                }
                return false;
            }
        }
    }
}
