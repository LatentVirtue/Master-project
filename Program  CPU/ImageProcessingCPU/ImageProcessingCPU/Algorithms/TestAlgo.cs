using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;

namespace ImageProcessingCPU.Algorithms
{
    class TestAlgo : IAlgoInterface
    {
        public Image Apply(Image x)
        {

            Bitmap temp = new Bitmap(ImageHandler.Temporary);
            //Bitmap r = temp.Clone(new Rectangle(0, 0, temp.Width, temp.Height), System.Drawing.Imaging.PixelFormat.Format16bppGrayScale);
            
            return temp;
        }
    }
}
