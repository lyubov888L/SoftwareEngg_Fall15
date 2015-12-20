#ifndef _FFTCTRecursive_H_
#define _FFTCTRecursive_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
void applyFFT(vector<complex<double> >& a_fhat, const vector<complex<double> >& a_f, int a_level);
void applyFFTInv(vector<complex<double> >& a_f, const vector<complex<double> >& a_fhat, int a_level);
vector<complex<double> > evenPoints(const vector<complex<double> >& a_vec);
vector<complex<double> > oddPoints(const vector<complex<double> >& a_vec);
#endif
