function S_dB = lg(S)
    % S - Power spectrum in linear scale
    % S_dB - Power spectrum in logarithmic scale
    % S_dB = 10*log10(S);
    S_dB = 10*log10(S);
end