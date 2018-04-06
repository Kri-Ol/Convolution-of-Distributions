using System;
using System.Diagnostics;

namespace Convolution
{
    // implemented according to https://en.wikipedia.org/wiki/PERT_distribution
    class PERT
    {
        public static double PDF(double x, double min, double mod, double max)
        {
            Debug.Assert(min < max);
            Debug.Assert(mod <= max);
            Debug.Assert(mod >= min);

            if (x < min)
                return 0.0;

            if (x > max)
                return 0.0;

            var range = max - min;

            var alpha = (4.0*mod + max - 5.0*min)/range;
            var beta  = (5.0*max - min - 4.0*mod)/range;

            var v = Math.Pow(x - min, alpha-1.0) * Math.Pow(max - x, beta-1.0)/(Utils.Beta(alpha, beta)*Math.Pow(range, alpha+beta-1.0));
            return v;
        }

        public static double CDF(double x, double min, double mod, double max)
        {
            throw new System.NotImplementedException("CDF");
        }
    }
}
