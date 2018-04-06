using System;
using System.Numerics;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

using MathNet.Numerics.IntegralTransforms;

namespace Convolution
{
    class Program
    {
        public static void TestPERT()
        {
            // compare with graph at https://en.wikipedia.org/wiki/PERT_distribution
            var min = 0.0;
            var mod = 10.0;
            var max = 100.0;

            for(int k = 0; k != 101; ++k) {
                var x = (double)k * 1.0;
                var f = PERT.PDF(x, min, mod, max);
                Console.WriteLine($"   {x}   {f}");
            }
        }

        public static void TestMdN()
        {
            Matrix<double> A = DenseMatrix.OfArray(new double[,] {
                                    {1,1,1,1},
                                    {1,2,3,4},
                                    {4,3,2,1}});
            MathNet.Numerics.LinearAlgebra.Vector<double>[] nullspace = A.Kernel();

            // verify: the following should be approximately (0,0,0)
            var q = (A * (2*nullspace[0] - 3*nullspace[1]));
            Console.WriteLine($"{q[0]}   {q[1]}   {q[2]}");
        }

        public static void TestU01()
        {
            for(int k = 0; k != 512; ++k)
            {
                var x = -2.06 + (double)k * 0.01;
                var f = U01.PDF(x);
                Console.WriteLine($"   {x}   {f}");
            }
        }

        public static void TestFourierU01()
        {
            Complex[] data = new Complex[512];

            for(int k = 0; k != 512; ++k)
            {
                var x = -2.06 + (double)k * 0.01;
                var f = U01.PDF(x);
                data[k] = new Complex(f, 0.0);
            }

            Fourier.Forward(data);
            Fourier.Inverse(data);

            for(int k = 0; k != 512; ++k)
            {
                Console.WriteLine($"   {data[k].Real}   {data[k].Imaginary}");
            }
        }

        public static void ConvolveU01_U01(int N)
        {
            double[]  x    = new double[N];
            Complex[] data = new Complex[N];

            double min = U01.min();
            double max = 2.0*U01.max();

            double step = (max - min)/(double)(N-1);

            for(int k = 0; k != N; ++k)
            {
                x[k]    = 0.0 + (double)k * step;
                data[k] = new Complex(U01.PDF(x[k]), 0.0);
            }

            Fourier.Forward(data, FourierOptions.Default);

            // var invN = 1.0 / Math.Sqrt(N);
            var invN = 1.0;
            // FT of convolution of the same PDSs is just a multiplication
            for(int k = 0; k != N; ++k)
            {
                var v = data[k];
                data[k] = v * v * invN;
            }

            Fourier.Inverse(data, FourierOptions.Default);

            // something is wrong with normalization
            for(int k = 0; k != N; ++k)
            {
                Console.WriteLine($"  {x[k]}  {data[k].Real}");
            }
        }

        public static void ConvolvePERT_PERT(int N)
        {
            double[]  x     = new double[N];

            var min = 0.0;
            var mod = 20.0; // peak at the left
            var max = 100.0;
            Complex[] data1 = new Complex[N];
            for(int k = 0; k != N; ++k)
            {
                x[k]     = 0.0 + (double)k * 1.0;
                data1[k] = new Complex(PERT.PDF(x[k], min, mod, max), 0.0);
            }

            min = 0.0;
            mod = 80.0; // peak at the right
            max = 100.0;
            Complex[] data2 = new Complex[N];
            for(int k = 0; k != N; ++k)
            {
                data2[k] = new Complex(PERT.PDF(x[k], min, mod, max), 0.0);
            }

            Fourier.Forward(data1, FourierOptions.Default);
            Fourier.Forward(data2, FourierOptions.Default);

            // FT of convolution of the PDSs is just a multiplication
            Complex[] data = new Complex[N];
            for(int k = 0; k != N; ++k)
            {
                data[k] = data1[k] * data2[k];
            }

            Fourier.Inverse(data, FourierOptions.Default);

            // something is wrong with normalization, mmm...
            for(int k = 0; k != N; ++k)
            {
                Console.WriteLine($"  {x[k]}  {data[k].Real}");
            }
        }

        public static Complex[] makeDelta(int N)
        {
            var r = new Complex[N];
            for(int k = 0; k != N; ++k)
            {
                r[k] = new Complex(0.0, 0.0);
            }
            r[0] = new Complex(1.0, 0.0);
            return r;
        }

        public static void TestDelat(int N)
        {
            var r = makeDelta(32);
            Fourier.Forward(r, FourierOptions.Default);
            Fourier.Inverse(r, FourierOptions.Default);

            foreach(var v in r) {
                Console.WriteLine($"  {v.Real}  {v.Imaginary}");
            }
        }

        static void Main(string[] args)
        {
            ConvolveU01_U01(401);
            // ConvolvePERT_PERT(201);
        }
    }
}
