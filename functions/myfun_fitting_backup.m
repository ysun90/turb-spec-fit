function F = myfun_fitting_16p_QC(x,f,S)

% global nfft
% persistent winlin
% if isempty(winlin)
%     winlin=zeros(nfft,1);
%     winlin(1:100) = 0;%linspace(0,1,100);
%     winlin(101:490) = 1;
%     winlin(491:535) = 0;%linspace(1,0.5,52)
%     winlin(536:925) = 1;
%     winlin(926:1023) = 0;%linspace(1,0,100);
% end

% GD + GGD + VD
F = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ) + ... % Central Spike (CS)
    x(4).*exp( -( ( abs( f-x(5) )./(x(6).*sqrt(gamma(1/x(7))/gamma(3/x(7)))) ).^x(7) ) ) + ... % Low Frequency (LF)
    x(8)*exp( -0.5*( ( ( f+x(9) )/x(10) ).^2 ) ) +... % Quasi-coherent Negetive (QCN)
    x(11)*exp( -0.5*( ( ( f-x(9) )/x(10) ).^2 ) ) +... % Quasi-coherent Positive (QCP)
    x(12)*voigt(f,x(13),x(14),x(15)) + ... % Broadband (BB)
    x(16); % Noise (N)

% Cost Function
F = sqrt(((10*log10(F)-10*log10(S)).^2./sum(10*log10(S).^2)) + 0.05*(F-S).^2);

end