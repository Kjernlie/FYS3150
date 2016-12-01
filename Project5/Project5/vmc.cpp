#include "vmc.h"
#include <iostream>
#include <armadillo>
double ran2(long *);


VMC::VMC() :
    nDimensions(3),
    omega(1),
    stepLength(1.0),
    nParticles(2),
    idum(-1),
    alpha(0.5*omega),
    nCycles(100)
{
}


void VMC::runMCIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double waveFuncOld = 0;
    double waveFuncNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE;

    for (int i = 0; i < nParticles; i++){
        for (int j=0; j < nDimensions, j++){
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // Create a loop over the Monte Carlo cycles here
}

double VMC::localEnergy(const mat &r)
{
    double R1_sqrd = r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2);
    double R2_sqrd = r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2);
    double R12 = sqrt(abs(R1_sqrd-R2_sqrd));
    local_energy = 0.5*omega*omega*(R1_sqrd+R2_sqrd)*(1-alpha*alpha) + 3*alpha*omega + 1/R12;
    return local_energy
}

double VMC::waveFunc(const mat &r)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++){
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            rSingleParticle += r(i,j)*r(i,j);
        }
    }
}




// Random number generator


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
