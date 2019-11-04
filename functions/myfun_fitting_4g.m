function F = myfun_fitting_4g(x,f,S)

% GD + GD + TD
F = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ) + ... % Central Spike (CS)
    x(4)*exp( -0.5*( ( ( f-x(5) )/x(6) ).^2 ) ) + ... % Low Frequency (LF)
    x(7).*exp( -( ( abs( f-x(8) )./(x(9).*sqrt(gamma(1/x(10))/gamma(3/x(10)))) ).^x(10) ) ) + ... % Broadband (BB)
    x(11); % Noise (N)

% Cost Function
F = sqrt(((10*log10(F)-10*log10(S)).^2./sum(10*log10(S).^2)) + 0.5*(F-S).^2); 

%F = sum(F.^2);

end
%x(4).*exp( -( ( abs( f-x(5) )./(x(6).*sqrt(gamma(1/x(7))/gamma(3/x(7)))) ).^x(7) ) )