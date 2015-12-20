#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <cassert>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFTW1D.H"

FFTW1D::FFTW1D(const unsigned int& a_M):FFT1D(a_M)
{
};
FFTW1D::~FFTW1D(){};
void FFTW1D::forwardFFTCC(vector<complex<double> > & a_fHat,
                            const vector<complex<double> >& a_f) const
{
  assert (a_fHat.size()==m_N);
  fftw_complex *in, *out;
  fftw_plan my_plan;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_N);

  for (unsigned int j = 0; j < m_N; j++)
    {
      in[j][0] = real(a_f[j]);
      in[j][1] = imag(a_f[j]);
    }

  my_plan = fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(my_plan);

  for (unsigned int j = 0; j < m_N; j++)
    {
      a_fHat[j] = complex<double>(out[j][0], out[j][1]);
    }

  fftw_destroy_plan(my_plan);
  fftw_free(in);
  fftw_free(out);

}
void FFTW1D::inverseFFTCC(vector<complex<double> >& a_f ,
                                  const vector<complex<double> > & a_fHat ) const
{
  assert (a_fHat.size()==m_N);
  fftw_complex *in, *out;
  fftw_plan my_plan;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_N);

  for (unsigned int j = 0; j < m_N; j++)
    {
      in[j][0] = real(a_fHat[j]);
      in[j][1] = imag(a_fHat[j]);
    }

  my_plan = fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(my_plan);

  for (unsigned int j = 0; j < m_N; j++)
    {
      a_f[j] = complex<double>(out[j][0], out[j][1]);
    }

  fftw_destroy_plan(my_plan);
  fftw_free(in);
  fftw_free(out);

}
