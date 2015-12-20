#ifndef _FFTW1D_H_
#define _FFTW1D_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "fftw3.h"
class FFTW1D:public FFT1D
{
public:
  FFTW1D(const unsigned int& a_M);
  virtual ~FFTW1D();
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat,
                          const vector<complex<double> >& f) const;
  virtual void inverseFFTCC(vector<complex<double> > & a_f,
                          const vector<complex<double> > & a_fHat) const;
};
#endif
