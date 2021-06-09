using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
namespace ImageProcessingCPU.Algorithms
{
    interface IAlgoInterface
    {
        public Image Apply(Image x);
    }
}
