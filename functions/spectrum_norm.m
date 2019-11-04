function [f,pxx] = spectrum_norm(varargin)

global nfft

% NORMALIZATION SPECTRUM
% [F,PXX] = SPECTRUM_NORM(I) returns spectrum 
% according to index i. 
% [F,PXX] = SPECTRUM_NORM(N,C,S) returns spectrum 
% according to tab of numchoc,channel,spec 

if nargin == 1
    [numc,channel,spec] = ind2tab(varargin{:});
elseif nargin == 3
    [numc,channel,spec] = varargin{:};
elseif nargin ~= 1 && nargin ~= 3
    error('Number of input argument can only be 1 or 3!!');
end

% FFT parameters
window = 1024;
noverlap = window/2;
nfft = 1025;
fs = 1e6;

% Load original spctrum data
S = load(['#',num2str(numc),'drefluc.mat']);

if size(S.txc,1)<window
    window = size(S.txc,1);
end

% Spectrum from channel 1 or 2
if channel==1
  [pxx,f] = pwelch(S.xc(:,spec)-mean(S.xc(:,spec)),window,noverlap,nfft,fs);
elseif channel==2
  [pxx,f] = pwelch(S.xc_2(:,spec)-mean(S.xc_2(:,spec)),window,noverlap,nfft,fs);
end
f = linspace(-max(f)/2,max(f)/2,nfft);
f = f'/1e3; % KHz
pxx = fftshift(pxx);

% Normalization
AreaSpec = trapz(f,pxx);
pxx = pxx/AreaSpec;

end