function HWHM = hwhm(F,S)
% F - frequency
% S - Spectrum
% HWHM - Half width at half maximum
S=S/trapz(F,S);
a = find(S>=(max(S)/2));
HWHM = (F(a(end))-F(a(1)))/2;
end