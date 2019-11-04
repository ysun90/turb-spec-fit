function [f,pxx] = spectrum_plot(varargin)

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
nfft = window+1;
fs = 1e6;

% Load spctrum data for one shot
S = load(['#',num2str(numc),'drefluc.mat']);

% Spectrum from channel 1 or 2
if channel==1
  [pxx,f] = pwelch(S.xc(:,spec)-mean(S.xc(:,spec)),window,noverlap,nfft,fs);
elseif channel==2
  [pxx,f] = pwelch(S.xc_2(:,spec)-mean(S.xc_2(:,spec)),window,noverlap,nfft,fs);
end
f = linspace(-max(f)/2,max(f)/2,nfft);
f = f'/1e3; % KHz
pxx = fftshift(pxx);

% Plot spectrum and log spectrum
load('database_r.mat');
figure('Position',[28         165        1265         487],'color','w');
grid on
subplot(121);
plot(f,pxx,'b','LineWidth',2);
grid on
grid minor
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
xlabel('Frequency [kHz]');
ylabel('Power spectral density');
set(gca,'FontSize',20);
idx_spec = tab2ind(numc,channel,spec);
title(sprintf('\\rho = %3.2f',database_r.r(idx_spec)));
subplot(122);
plot(f,10*log10(pxx/trapz(f,pxx)),'r','LineWidth',2);
grid on
grid minor
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
xlabel('Frequency [kHz]');
ylabel('Power spectral density [dB]');
title(sprintf('\\rho = %3.2f, normalized',database_r.r(idx_spec)));
% title(['@',num2str(idx_spec),', \rho = ',num2str(database_r.r(idx_spec),2)]);

end