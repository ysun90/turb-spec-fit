function spectra(numc,varargin)

% SPECTRA(N) plots all the spectra in shot N. 
% SPECTRA(N,1) plots spectra in ch1 of shot N.
% SPECTRA(N,2) plots spectra in ch2 of shot N.

if nargin < 1 || nargin > 2
    error('Number of Input arguments can only be 1 or 2');
end

% FFT parameters
window = 1024;
noverlap = window/2;
nfft = 1025;
fs = 1e6;

% Load data 
data = load(['#',num2str(numc),'drefluc.mat']);

% Plot spectra for ch1
if nargin==1 || varargin{1}==1    
    for ii = 1:length(data.R1)        
        [pxx,f] = pwelch(data.xc(:,ii)-mean(data.xc(:,ii)),window,noverlap,nfft,fs);
        f = linspace(-max(f)/2,max(f)/2,nfft);
        f = f'/1e3;
        pxx = fftshift(pxx);
        h = figure;
        set(h,'Position',[10   138   668   500]);
        plot(f,10*log10(pxx/max(pxx)));
        xlabel('Frequency (kHz)');
        ylabel('PSD (dB)');
        title(['#',num2str(data.numc),', rc = ',num2str(data.r1(ii),3),...
            ', t = ',num2str(data.txc(1,ii))]);
        grid on;
    end
end
  
% Plot spectra for ch2
if nargin==1 || varargin{1}==2   
    if isempty(data.xc_2)
        errordlg('Sorry, there is no ch2 in this shot!','Input Wrong');
    end    
    for ii = 1:length(data.R2)       
        [pxx,f] = pwelch(data.xc_2(:,ii)-mean(data.xc_2(:,ii)),window,noverlap,nfft,fs);
        f = linspace(-max(f)/2,max(f)/2,nfft);
        f = f'/1e3;
        pxx = fftshift(pxx);
        h = figure;
        set(h,'Position',[698   138   668   500]);
        plot(f,10*log10(pxx/max(pxx)));
        xlabel('Frequency (kHz)');
        ylabel('PSD (dB)');
        title(['#',num2str(data.numc),', rc = ',num2str(data.r2(ii),3),...
            ', t = ',num2str(data.txc(1,ii))]);
        grid on;
    end
end

end