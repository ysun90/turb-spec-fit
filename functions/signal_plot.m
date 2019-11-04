function [t,x,f,pxx] = signal_plot(varargin)

% Check input parameters
if nargin == 1
    [numc,channel,spec] = ind2tab(varargin{:});
elseif nargin == 3
    [numc,channel,spec] = varargin{:};
else
    error('Number of input arguments can only be 1 or 3!!');
end

% Load spctrum data for one shot
S = load(['#',num2str(numc),'drefluc.mat']);

% Get time and signal
t = S.txc(:,spec);
if channel==1
    x = S.xc(:,spec);
elseif channel==2
    x = S.xc_2(:,spec);
end

% FFT by Welch method
window = 1024;
noverlap = window/2;
nfft = 1024;
fs = 1e6;
[pxx,f] = pwelch((x-mean(x)),window,noverlap,nfft,fs,'centered');
f = f/1e3; % KHz
pxx = pxx/max(pxx); % Normalized to maximam power

% Plot the signal and spectrum
figure('Position',[95         116        1181         550]);
% Amplitude
ax1 = subplot(4,2,1);
plot(t,abs(x))
ax1.XLim = [t(1) t(end)];
title('Amplitude')
% Phase
ax1 = subplot(4,2,3);
plot(t,angle(x))
ax1.XLim = [t(1) t(end)];
title('Phase')
% Real
ax1 = subplot(4,2,5);
plot(t,real(x))
ax1.XLim = [t(1) t(end)];
title('Real')
% Imaginary
ax1 = subplot(4,2,7);
plot(t,imag(x))
ax1.XLim = [t(1) t(end)];
title('Imaginary')
% FFT Spectrum
subplot(4,2,[2,4]);
grid on
plot(f,10*log10(pxx));
ax = gca;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
%xlabel('f [KHz])');
title('FFT Spectrum');
% CWT (Continuous wavelet transform) for amplitude
subplot(4,2,[6,8]);
cwt(abs(x),'amor',1e6) % 'bump'

end