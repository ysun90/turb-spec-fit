function F = myfun_CS_LF_BB_N_12p_Complex_backup(x,f,S)

global nfft
persistent winlf winlin

if isempty(winlin)
    winlin=zeros(nfft,1);
    winlin(1:100) = 0;%linspace(0,1,100);
    winlin(101:490) = 1;
    winlin(491:535) = 0;%linspace(1,0.5,52)
    winlin(536:925) = 1;
    winlin(926:1025) = 0;%linspace(1,0,100);
end

if isempty(winlf)
    winlf=zeros(nfft,1);
    winlf(1:490) = 0;
    winlf(491:510) = 1;
    winlf(511:515) = 0;
    winlf(516:550) = 1;
    winlf(551:end) = 0;
end

% % Fitting model with 12 parameters
% F = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ) + ... % Central Spike (CS) 
%     x(4)*exp( -( ( abs( f-x(5) )./x(6) ).^x(7) ) ) + ... % Low Frequency (LF)
%     x(8)*exp( -( ( abs( f-x(9) )./x(10) ).^x(11) ) ) + ... % Broadband (BB)
%     x(12); % Noise (N)

% Fitting model with 12 parameters
F = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ) + ... % Central Spike (CS) 
    x(4).*exp( -( ( abs( f-x(5) )./(x(6).*sqrt(gamma(1/x(7))/gamma(3/x(7)))) ).^x(7) ) ) + ... % Low Frequency (LF)
    x(8).*exp( -( ( abs( f-x(9) )./(x(10).*sqrt(gamma(1/x(11))/gamma(3/x(11)))) ).^x(11) ) ) + ... % Broadband (BB)
    x(12); % Noise (N)

% F = sqrt((log10(F)-log10(S)).^2./sum(log10(S).^2) + winlin.*(F-S).^2); % default
% F = sqrt((log10(F)-log10(S)).^2./sum(log10(S).^2) + 0*(F-S).^2); % default
F = sqrt(((log10(F)-log10(S)).^2./sum(log10(S).^2)) + 0*(F-S).^2); % default

end