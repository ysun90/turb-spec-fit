function [y] = voigtNorm( wavenumberArray,centerLine,widthGauss, widthLorentz )
% The function calculates VOIGT profile using 
% an algorithm kindly provided by Dr. F. Schreier. (in FORTRAN)
% and rewritten by Dr. N. Cherkasov for MATLAB
% see the function complexErrorFunction for details
%
%  intensity at the maximum of this function is equal to 1
% 
% the reference:
% F. Schreier: Optimized Implementations of Rational Approximations for the Voigt ane Complex Error Function. 
% J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 1010–1025, doi 10.1016/j.jqsrt.2010.12.010, 2011. 
%
% INPUT arguments
% wavenumberArray - array 1xN of wavenumbers 
% centerLine - position of the center line of the band
% widthGauss - parameter of the width of the Gaussian component (HWHM)
% widthLorentz - parameter of the width of the Lorentzian component (HWHM)

% converting to dimensionless coordinates
x=sqrt(log(2)).*(wavenumberArray-centerLine)./(widthGauss);
y=sqrt(log(2)).*(widthLorentz/widthGauss);

w=complexErrorFunction(x,y);
y=real(w);
y=y/max(y);

end

