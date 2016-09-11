#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

// Try to make a function

double function(double x)
{
  double result = atan(x);
  return result;
}

//double first_order(double (*function)(double x),double x, double h)
double first_order(double x, double h)
{
  double result = (function(x+h)-function(x))/h;
  //double result = (*function(x+h)-*function(x))/h;
  return result;
}

double second_order(double x, double h)
{
  double result = (function(x+h)-function(x-h))/(2*h);
  return result; 
}

int main()
{
  double d10;
  double const exact=1./3;
  double const x1 = sqrt (2);

  int const  N = 15;
  double f1[N+1];
  double f2[N+1];

  double error_f1[N+1];
  double error_f2[N+1];
  
  ofstream myfile;
  myfile.open("warmup.csv");
  for(int i = 0 ; i <= N ; i++)
  {
    f1[i] = first_order(x1,pow(10,-i));
    f2[i] = second_order(x1,pow(10,-i));
    error_f1[i] = abs(f1[i]-exact);
    error_f2[i] = abs(f2[i]-exact);
    myfile << error_f1[i] << ", " << error_f2[i] << "\n";
  }

  myfile.close();

  return 0;
}

