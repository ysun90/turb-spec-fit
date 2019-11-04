function [mu,sigma,kts ]= moment124_pdf(F,S)
% kts = kurtosis_pdf(F,S)
% The kurtosis is the fourth standardized moment,
    
    S=S/sum(S); % Normalization
%     mu = mean_pdf(F,S); % Mean
%     sigma = std_pdf(F,S); % Standard deviation
%     kts = sum(((F-mu)).^4.*S)/sigma.^4; % Kurtosis
    mu=sum(F.*S);
    variance=sum((F-mu).^2.*S);
    sigma=sqrt(variance);
    kts=sum((F-mu).^4.*S)/ variance.^2;
   
end