using System.Drawing;
namespace ImageProcessingCPU.Algorithms
{
    interface IAlgoInterface
    {
        public Image Apply(ref Image x);
    }
}
