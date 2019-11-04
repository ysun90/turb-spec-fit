function [y] = voigtNormGLform(lorentzShare,width,centerLine,wavenumberArray)  
% voigtNormGLform   Calculation of the VOIGT profile 
%
%   [y] = voigtNormGLform(lorentzShare,width,centerLine,wavenumberArray)  
%   The function calculates the Voight profile using the algorithm 
%   kindly provided by Dr. F. Schreier in FORTRAN and rewritten to MATLAB
%   by Dr. N. Cherkasov
% 
%   This function is more convinient than VOIGT for deconvolution, 
%   because its intensity is 1 at the maximum. Moreover, input arguments
%   are different.
%
%   For more details on algorithm see the publication:
%   F. Schreier: Optimized Implementations of Rational Approximations for the Voigt ane Complex Error Function. 
%   J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 1010–1025, doi 10.1016/j.jqsrt.2010.12.010, 2011. 
%
%   The function was used for the deconvolution of IR spectra see the
%   publication and itsSupplementary
%
%
%   INPUT ARGUMENTS
%       lorentzShare - fraction of the Lorentzian component in the
%          deconvolution in PERCENT 0-100
%       width        - width parameter. It is close to the resulting
%          half-width at half maximum of the resulting band, but it is not equal
%       centerLine   - position of the band center
%       wavenumberArray - array 1*N of wavenumbers 
%
% 	OUTPUT
%       y - array 1*N of intensities
%
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com


% converting to dimensionless coordinates
if lorentzShare<=0, lorentzShare=1e-6; end;
if lorentzShare>=100, lorentzShare=100-1e-6; end;

widthGauss=(100-lorentzShare)/100*width;
widthLorentz=(lorentzShare)/100*width;

x=sqrt(log(2)).*(wavenumberArray-centerLine)./(widthGauss);
y=sqrt(log(2)).*(widthLorentz/widthGauss);

w=complexErrorFunction(x,y);
y=real(w);
y=y/max(y);

end


