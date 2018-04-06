using System;
using System.Diagnostics;

namespace Convolution
{
    // Uniform in [0...1] implemented according to https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    class U01
    {
        public static double PDF(double x)
        {
            if (x < 0.0)
                return 0.0;

            if (x > 1.0)
                return 0.0;

            return 1.0;
        }

        public static double CDF(double x)
        {
            if (x < 0.0)
                return 0.0;

            if (x > 1.0)
                return 0.0;

            return x;
        }

        public static double min()
        {
            return 0.0;
        }

        public static double max()
        {
            return 1.0;
        }
    }
}
