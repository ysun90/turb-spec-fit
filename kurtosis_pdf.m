function kts = kurtosis_pdf(F,S)
% kts = kurtosis_pdf(F,S)
% The kurtosis is the fourth standardized moment,
    
    %S=S/sum(S); % Normalization
    mu = mean_pdf(F,S); % Mean
    sigma = std_pdf(F,S); % Standard deviation
    kts = sum(((F-mu)/sigma).^4.*S); % Kurtosis

end