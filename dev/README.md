# Development version of PSE

This is the development version of the PSE code.

This version uses FFTs to compute nonlinear PSE but the FFT routines that
are supported are either CRAY or SGI and need to be re-implemented.

Note that this uses SGI FFT routines as well as outdated LINPACK routines
that need to be replaced with FFTW and LAPACK.

S. Scott Collis\
flow.physics.simulation@gmail.com