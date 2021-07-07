using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;

namespace ImageProcessingCPU.Algorithms
{
    class HoughTransform : IAlgoInterface
    {
        static Canny cannyEdge;
        Image actual;
        public HoughTransform()
        {

        }
        public Image Apply(ref Image x)
        {
            return actual;
        }
    }
}
