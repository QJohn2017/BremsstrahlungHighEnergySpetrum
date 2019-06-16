# Bremsstrahlung

Set of routines to compute the Bremsstrahlung emissivity from relativistic electrons.
The model is described in Appendix A of Strong, Moskalenko and Reimer, ApJ, 537, 763 (2000)

## To compile it:

> g++ -std=c++11 -I$GSL_DIR/include -L$GSL_DIR/lib -lgsl -lgslcblas bremsstrahlung_sigma.cpp 
