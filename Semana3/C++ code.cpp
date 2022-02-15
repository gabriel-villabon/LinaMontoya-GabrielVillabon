#include<math.h>

double function(double x)
{
    x = -x*x;
    double result = exp(x);
    return result;
}

double derivate(double x, double h)
{
    double derivate = (function(x+h) - function(x-h))/(2*h);
    return derivate;
}

int main()
{
    double x = 0;
    double h = 0.01;
    cout << "Escribe un numero: \n";
    cin >> x;
    cout << "El resultado es: " << derivate(x, h);
}