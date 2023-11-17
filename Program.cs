double[] doubles = new double[] {1, -1.5 };
double[] doubles2 = new double[] { 1, -2, 1 };
double[] doubles3 = new double[] { 1, 2, 2 };
double[] doubles4 = new double[] { 1,-4,-4,16 };

Console.WriteLine(string.Join(",", PolyRoots(doubles)));
Console.WriteLine(string.Join(",", PolyRoots(doubles2)));
Console.WriteLine(string.Join(",", PolyRoots(doubles3)));
Console.WriteLine(string.Join(",", PolyRoots(doubles4)));

double[] PolyRoots(double[] coefficients)
{
    List<double> roots = new List<double>();
    double[] quotient = coefficients;
    while (true)
    {
        double root = PolyRoot(quotient);
        if(root == double.NegativeInfinity)
        {
            return roots.ToArray();
        }
        roots.Add(root);
        quotient = PolyDiv(quotient, root);
    }

    double EvaluatingFunction(double[] coefficients, double x)
    {
        int degree = coefficients.Length - 1;
        double res = 0;
        for (int i = 0; i < coefficients.Length; i++)
        {
            res += coefficients[i] * Math.Pow(x, degree - i);
        }
        return res;
    }

    double[] EvaluatingDerivative(double[] coefficients)
    {
        int degree = coefficients.Length - 1;
        double[] derivativeCoefficients = new double[degree]; // initialize new array for derivative coefficients
        for (int i = 0; i < degree; i++)
        {
            derivativeCoefficients[i] = (degree - i) * coefficients[i]; // compute derivative coefficients
        }
        return derivativeCoefficients;
    }

    double PolyRoot(double[] coefficients)
    {
        double x0 = 0.151231;
        //double root = double.NegativeInfinity;
        double error = 0.000001;
        int maxIterations = 1000;

        for (int i = 0; i < maxIterations; i++)
        {
            double x1 = x0 - EvaluatingFunction(coefficients, x0) / EvaluatingFunction(EvaluatingDerivative(coefficients), x0);
            if (Math.Abs(x1 - x0) <= error) return x1;
            x0 = x1;
        }
        return double.NegativeInfinity;
    }

    double[] PolyDiv(double[] coefficients, double xi)
    {
        double[] quotient = new double[coefficients.Length - 1];
        double prevCoeff = 0;
        double result;
        for (int i = 0; i < quotient.Length; i++)
        {
            result = coefficients[i] + prevCoeff;
            quotient[i] = result;
            prevCoeff = xi * result;
        }
        return quotient;
    }
}
