function std = std_pdf_norm(F,S)
    
    S=S/sum(S); % Normalization
    mu = mean_pdf(F,S); % Mean
    var = sum((F-mu).^2.*S); % Variance
    std = sqrt(var); % Standard deviation
    
end