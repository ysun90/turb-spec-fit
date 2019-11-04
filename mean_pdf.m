function mu = mean_pdf(F,S)
    
    S=S/sum(S); % Normalization
    mu = sum(F.*S); % Mean
end