using System;
using System.Collections.Generic;
using System.Text;
using ILGPU;
using ILGPU.Runtime;
using ILGPU.Runtime.Cuda;
using ILGPU.Algorithms;
using System.Drawing;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

// Mark all internals to be visible to the ILGPU runtime
[assembly: InternalsVisibleTo(Context.RuntimeAssemblyName)]

namespace ImageProcessingCPU.Algorithms
{

    class CannyGPU
    {
        public struct TroupleDouble
        {
            public double D1 { get; set; }
            public double D2 { get; set; }
            public double D3 { get; set; }
            public TroupleDouble(double first, double second, double third)
            {
                D1 = first;
                D2 = second;
                D3 = third;
            }
            public TroupleDouble(Color c)
            {
                D1 = c.R;
                D2 = c.G;
                D3 = c.B;
            }
        }
        static readonly double[,] gaussianMatrix = {    {0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067 },
                                            {0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
                                            {0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
                                            {0.00038771, 0.01330373, 0.11098164, 0.22508352, 0.11098164, 0.01330373, 0.00038771 },
                                            {0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
                                            {0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
                                            {0.00000067, 0.00002292, 0.00019117, 0.00038771 ,0.00019117, 0.00002292, 0.00000067 }};
        //Sobel operator kernels
        static readonly double[,] sobelX = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
        static readonly double[,] sobelY = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };

        //Prewitt
        static readonly double[,] prewittX = { { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 } };
        static readonly double[,] prewittY = { { 1, 1, 1 }, { 0, 0, 0 }, { -1, -1, -1 } };

        //Roberts
        static readonly double[,] robertsX = { { 1, 0 }, { 0, -1 } };
        static readonly double[,] robertsY = { { 0, 1 }, { -1, 0 } };

        //Scharr
        static readonly double[,] scharrX = { { 3, 0, -3 }, { 10, 0, -10 }, { 3, 0, -3 } };
        static readonly double[,] scharrY = { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } };

        //Selected kernel
        static double[,] Xfilter;
        static double[,] Yfilter;

        //Image dimensions
        static int n; //height
        static int m; //width
        //Gradient max
        static double max = 0;
        //Treshold Cutoff Points
        static double lowerT;
        static double upperT;

        //ILGPU context and accelerator
        static void GrayscaleKernel(Index2 position, ArrayView2D<TroupleDouble> actual)
        {
            actual[position].D1 = (actual[position].D1 + actual[position].D2 + actual[position].D3) / 3.0;
        }
        static void GaussKernel(Index2 position, ArrayView2D<TroupleDouble> actual, ArrayView2D<double> gaussKernel, int n, int m)
        {
            double t = 0;
            long offset = (gaussKernel.Height - 1) / -2;
            for (long i = gaussKernel.Height - 1; i > -1; i--)
            {
                for (long j = gaussKernel.Width - 1; j > -1; j--)
                {
                    long x = position.X + i + offset;
                    double tt;
                    bool u = true;
                    if (x >= n || x < 0)
                    {
                        u = false;
                    }
                    long y = position.Y + j + offset;
                    if (y >= m || y < 0)
                    {
                        u = false;
                    }
                    if (u)
                    {
                        tt = actual[x, y].D1;
                    }
                    else
                    {
                        tt = actual[position].D1;
                    }
                    t += gaussKernel[i, j] * tt;
                }
            }
            actual[position].D2 = t;
        }
        static void GradientConvolutionKernel(Index2 position, ArrayView2D<TroupleDouble> actual, ArrayView2D<double> kernel, bool g, int n, int m)
        {
            double t = 0;
            long offset = (kernel.Height - 1) / -2;
            for (long i = kernel.Height - 1; i > -1; i--)
            {
                for (long j = kernel.Width - 1; j > -1; j--)
                {
                    long x = position.X + i + offset;
                    double tt;
                    bool u = true;
                    if (x >= n || x < 0)
                    {
                        u = false;
                    }
                    long y = position.Y + j + offset;
                    if (y >= m || y < 0)
                    {
                        u = false;
                    }
                    if (u)
                    {
                        tt = actual[x, y].D2;
                    }
                    else
                    {
                        tt = actual[position].D2;
                    }
                    t += kernel[i, j] * tt;
                }
            }
            if (g)
            {
                actual[position].D1 = t;
            }
            else
            {
                actual[position].D3 = t;
            }
        }
        static void ComputeGradientKernel(Index2 position, ArrayView2D<TroupleDouble> actual)
        {
            actual[position].D2 = (XMath.Atan2(actual[position].D1, actual[position].D3) * (180 / XMath.PI)) % 180;
            actual[position].D1 = XMath.Sqrt(actual[position].D1 * actual[position].D1 + actual[position].D3 * actual[position].D3);
        }
        static void NMSKernel(Index2 position, ArrayView2D<TroupleDouble> actual)
        {
            actual[position.X+1,position.Y+1].D3 = 0;
            double angle = actual[position.X + 1, position.Y + 1].D2;
            //linear interpolation
            if (angle > -45)
            {
                if (angle > 0)
                {
                    if (angle > 45)
                    {
                        if (angle > 90)
                        {
                            if (angle > 135)
                            {
                                //135 to 180 1
                                double l = (angle - 135) / 45;
                                double r = 1 - l;
                                double a = XMath.Max(actual[position.X, position.Y].D1 * l + actual[position.X, position.Y + 1].D1 * r, actual[position.X + 2, position.Y + 2].D1 * l + actual[position.X + 2, position.Y + 1].D1 * r);
                                if (actual[position.X+1,position.Y+1].D1 >= a)
                                {
                                    actual[position.X+1,position.Y+1].D3 = actual[position.X+1,position.Y+1].D1;
                                }
                            }
                            else
                            {
                                //90 to 135 1
                                //main diagonal
                                double l = (actual[position].D2 - 90) / 45;
                                double r = 1 - l;
                                double a = XMath.Max(actual[position.X + 1, position.Y + 2].D1 * l + actual[position.X + 2, position.Y + 2].D1 * r, actual[position.X + 1, position.Y].D1 * l + actual[position.X, position.Y].D1 * r);
                                if (actual[position.X + 1, position.Y + 1].D1 >= a)
                                {
                                    actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                                }
                            }
                        }
                        else
                        {
                            //45 to 90 1
                            //side diagonal
                            double l = (angle - 45) / 45;
                            double r = 1 - l;
                            double a = XMath.Max(actual[position.X + 1 - 1, position.Y + 2].D1 * l + actual[position.X + 1, position.Y + 2].D1 * r, actual[position.X + 2, position.Y].D1 * l + actual[position.X + 1, position.Y].D1 * r);
                            if (actual[position.X + 1, position.Y + 1].D1 >= a)
                            {
                                actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                            }
                        }
                    }
                    else
                    {
                        //0 to 45 1
                        //vertical right
                        double l = angle / 45;
                        double r = 1 - l;
                        double a = XMath.Max(actual[position.X, position.Y + 1].D1 * l + actual[position.X, position.Y + 2].D1 * r, actual[position.X + 2, position.Y + 1].D1 * l + actual[position.X + 2, position.Y].D1 * r);
                        if (actual[position.X + 1, position.Y + 1].D1 >= a)
                        {
                            actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                        }
                    }
                }
                else
                {
                    //-45 to 0 1
                    //vertical left
                    double l = (-45 - angle) / -45;
                    double r = 1 - l;
                    double a = XMath.Max(actual[position.X, position.Y].D1 * l + actual[position.X, position.Y + 1].D1 * r, actual[position.X + 2, position.Y + 2].D1 * l + actual[position.X + 2, position.Y + 1].D1 * r);
                    if (actual[position.X + 1, position.Y + 1].D1 >= a)
                    {
                        actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                    }
                }
            }
            else
            {
                if (angle < -90)
                {
                    if (angle < -135)
                    {
                        //-180 to -135 1
                        double l = (angle + 135) / -45;
                        double r = 1 - l;
                        double a = XMath.Max(actual[position.X, position.Y + 1].D1 * r + actual[position.X, position.Y + 2].D1 * l, actual[position.X + 2, position.Y + 1].D1 * r + actual[position.X + 2, position.Y].D1 * l);
                        if (actual[position.X + 1, position.Y + 1].D1 >= a)
                        {
                            actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                        }
                    }
                    else
                    {
                        //-135 to -90 1
                        double l = (angle + 90) / -45;
                        double r = 1 - l;
                        double a = XMath.Max(actual[position.X + 1, position.Y + 2].D1 * l + actual[position.X, position.Y + 2].D1 * r, actual[position.X + 1, position.Y].D1 * l + actual[position.X + 2, position.Y].D1 * r);
                        if (actual[position.X + 1, position.Y + 1].D1 >= a)
                        {
                            actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                        }
                    }
                }
                else
                {
                    //-90 to -45
                    double l = (angle + 45) / -45;
                    double r = 1 - l;
                    double a = XMath.Max(actual[position.X + 1 - 1, position.Y].D1 * l + actual[position.X + 1, position.Y].D1 * r, actual[position.X + 2, position.Y + 2].D1 * l + actual[position.X + 1, position.Y + 2].D1 * r);
                    if (actual[position.X + 1, position.Y + 1].D1 >= a)
                    {
                        actual[position.X + 1, position.Y + 1].D3 = actual[position.X + 1, position.Y + 1].D1;
                    }
                }
            }
        }
        static void TresholdKernel(Index2 position, ArrayView2D<TroupleDouble> actual, ArrayView2D<byte> final, double max, double lowerT, double upperT)
        {
            if (actual[position].D3 < lowerT * max)
            {
                //weak edges, which are disgarded
                final[position] = 0;
            }
            else if (actual[position].D3 < upperT * max)
            {
                //edges for which we use hysteresis
                final[position] = 1;
            }
            else
            {
                //strong edges, which are always kept
                final[position] = 2;
            }
        }
        public static byte[,] ApplyAll(int selectedKernel, TroupleDouble[,] actual, double tLow, double tHigh)
        {
            switch (selectedKernel)
            {
                case 1:
                    Xfilter = prewittX;
                    Yfilter = prewittY;
                    break;
                case 2:
                    Xfilter = robertsX;
                    Yfilter = robertsY;
                    break;
                case 3:
                    Xfilter = scharrX;
                    Yfilter = scharrY;
                    break;
                default:
                    Xfilter = sobelX;
                    Yfilter = sobelY;
                    break;
            }
            n = actual.GetLength(0);
            m = actual.GetLength(1);
            lowerT = tLow;
            upperT = tHigh;
            //init
            using var context = new Context();
            context.EnableAlgorithms();
            using var cuda = new CudaAccelerator(context);
            //grayscale
            //saves values to D1
            var k1 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>>(GrayscaleKernel);
            using var bufferMain = cuda.Allocate(actual);
            k1(new Index2(n, m), bufferMain.View);
            cuda.Synchronize();
            k1 = null;
            //gauss
            //saves values to D2
            var k2 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>, ArrayView2D<double>, int, int>(GaussKernel);
            var bufferHelper = cuda.Allocate(gaussianMatrix);
            k2(new Index2(n, m), bufferMain.View, bufferHelper.View, n, m);
            cuda.Synchronize();
            k2 = null;
            //gradient
            //Gx is on D1, Gy is on D3
            var k3 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>, ArrayView2D<double>, bool, int, int>(GradientConvolutionKernel);
            //Gx
            bufferHelper = cuda.Allocate(Xfilter);
            k3(new Index2(n, m), bufferMain.View, bufferHelper.View, true, n, m);
            cuda.Synchronize();
            //Gy
            using var bufferHelper1 = cuda.Allocate(Yfilter);
            k3(new Index2(n, m), bufferMain.View, bufferHelper1.View, false, n, m);
            cuda.Synchronize();
            k3 = null;
            //computeGradient
            //Gradient intensity is on D1, gradient angle is on D2
            k1 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>>(ComputeGradientKernel);
            k1(new Index2(n, m), bufferMain.View);
            cuda.Synchronize();
            //k1 = null;
            //Non-Max suppression
            //Gradient Intensity on D3
            k1 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>>(NMSKernel);
            k1(new Index2(n - 2, m - 2), bufferMain.View);
            cuda.Synchronize();
            k1 = null;
            //Double Tresholding, output to final
            var k4 = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<TroupleDouble>, ArrayView2D<byte>, double, double, double>(TresholdKernel);
            using var bufferFinal = cuda.Allocate(new byte[n, m]);
            //determine max
            actual = bufferMain.GetAs2DArray();
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (actual[i, j].D3 > max)
                    {
                        max = actual[i, j].D3;
                    }
                }
            }
            k4(new Index2(n, m), bufferMain.View, bufferFinal.View, max, lowerT, upperT);
            cuda.Synchronize();
            //return final
            byte[,] r = bufferFinal.GetAs2DArray();
            return r;
        }
        static void ConvolutionKernelInt(Index2 position, ArrayView2D<int> kernel, ArrayView2D<int> target, ArrayView2D<int> final)
        {
            int t = 0;
            for (long i = kernel.Height - 1; i > -1; i--)
            {
                for (long j = kernel.Width - 1; j > -1; j--)
                {
                    long x = position.X + (kernel.Width - j - 1);
                    if (x >= target.Width)
                    {
                        continue;
                    }
                    long y = position.Y + (kernel.Height - i - 1);
                    if (y >= target.Height)
                    {
                        continue;
                    }
                    t += kernel[j, i] * target[x, y];
                }
            }
            final[position] = t;
        }
        static void ConvolutionKernelDouble(Index2 position, ArrayView2D<double> kernel, ArrayView2D<double> target, ArrayView2D<double> final)
        {
            double t = 0;
            long offset = (kernel.Height - 1) / -2;
            for (long i = kernel.Height - 1; i > -1; i--)
            {
                for (long j = kernel.Width - 1; j > -1; j--)
                {
                    long x = position.X + i + offset;
                    double tt;
                    bool u = true;
                    if (x >= target.Width || x < 0)
                    {
                        u = false;
                    }
                    long y = position.Y + j + offset;
                    if (y >= target.Height || y < 0)
                    {
                        u = false;
                    }
                    if (u)
                    {
                        tt = target[x, y];
                    }
                    else
                    {
                        tt = target[position];
                    }
                    t += kernel[i, j] * tt;
                }
            }
            final[position] = t;
        }
        public static void ApplyInt(int[,] kernel, ref int[,] target)
        {
            //init
            using var context = new Context();
            using var cuda = new CudaAccelerator(context);
            var k = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<int>, ArrayView2D<int>, ArrayView2D<int>>(ConvolutionKernelInt);
            using var bufferTarget = cuda.Allocate<int>(target);
            using var bufferFinal = cuda.Allocate<int>(new int[target.GetLength(0), target.GetLength(1)]);
            using var bufferKernel = cuda.Allocate<int>(kernel);
            k(new Index2(target.GetLength(0), target.GetLength(1)), bufferKernel.View, bufferTarget.View, bufferFinal.View);
            cuda.Synchronize();
            target = bufferFinal.GetAs2DArray();
        }
        public static void ApplyDouble(double[,] kernel, ref double[,] target)
        {
            //init
            using var context = new Context();
            using var cuda = new CudaAccelerator(context);
            var k = cuda.LoadAutoGroupedStreamKernel<Index2, ArrayView2D<double>, ArrayView2D<double>, ArrayView2D<double>>(ConvolutionKernelDouble);
            using var bufferTarget = cuda.Allocate<double>(target);
            using var bufferFinal = cuda.Allocate<double>(new double[target.GetLength(0), target.GetLength(1)]);
            using var bufferKernel = cuda.Allocate<double>(kernel);
            k(new Index2(target.GetLength(0), target.GetLength(1)), bufferKernel.View, bufferTarget.View, bufferFinal.View);
            cuda.Synchronize();
            target = bufferFinal.GetAs2DArray();
        }
    }
}
