
function sdS = sd_spec(F,S)
% F - frequency
% S - Spectrum
% SD - Standard deviation of the spectrum

% Normalization
S_norm = S/sum(S);

% Mean
muS = trapz(F,F.*S_norm);

% Variance
varS = trapz(F,(F-muS).^2.*S_norm);

% Standard deviation 
sdS = sqrt(varS);
end