function [f,pxx,r,R,F,numc,channel,spec] = spectrum_original(varargin)

global nfft

% SPECTRUM(I) plots the spectrum of index i. 
% SPECTRUM(N,C,S) plots the spectrum of numc,channel,spec 
% [F,PXX] = SPECTRUM(I) returns spectrum of index i. 
% [F,PXX] = SPECTRUM(N,C,S) returns spectrum of numc,channel,spec 

if nargin == 1
    [numc,channel,spec] = ind2tab(varargin{:});
elseif nargin == 3
    [numc,channel,spec] = varargin{:};
else
    error('Number of input arguments can only be 1 or 3!!');
end

% FFT parameters
window = 1024;
noverlap = window/2;
nfft = 1025;
fs = 1e6;

% Load spctrum data for one shot
S = load(['#',num2str(numc),'drefluc.mat']);

% Spectrum from channel 1 or 2
if channel==1
  [pxx,f] = pwelch(S.xc(:,spec)-mean(S.xc(:,spec)),window,noverlap,nfft,fs);
  r=S.r1(spec);R=S.R1(spec);F=S.F1(spec);
elseif channel==2
  [pxx,f] = pwelch(S.xc_2(:,spec)-mean(S.xc_2(:,spec)),window,noverlap,nfft,fs);
  r=S.r2(spec);S.R=S.R2(spec);S.F=S.F2(spec);
end
f = linspace(-max(f)/2,max(f)/2,nfft);
f = f'/1e3; % KHz
pxx = fftshift(pxx);

end