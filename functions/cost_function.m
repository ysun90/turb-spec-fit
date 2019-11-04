function Fcost = cost_function(F,S,Sfit,method)
% Sfit - fitting function
% S = spectrum
% w - weight in linear scale

%global wlinear wcenter

% Weight for linear scale
wlinear=0.5;
wcenter = 0.2;

% Remove central point (513)
Weight = ones(length(F),1); Weight(513) = wcenter;

% Cost function
switch method
    case {'SL','SG','SM','VLS','vector'} % For least squared vector
        Fcost = sqrt(  ( (1-wlinear)*(lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)*1000/1025) + wlinear*(Sfit-S).^2 ).*Weight  ); % default
        %Fcost = sqrt(  (lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)*1000/1025) + w*(Sfit-S).^2  ); % default
    case {'AL','AG'} % For absolute vector
        Fcost = sqrt( abs(lg(Sfit)-lg(S))/(sum(abs(lg(S)))*1000/1025) + wlinear*abs(Sfit-S));
    case {'local','multi','global'} % for fmincon
        Fcost = ( (1-wlinear)*(lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)*1000/1025) + wlinear*(Sfit-S).^2 ).*Weight  ;
        %Fcost = sum(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2); % default
end
% Cost Function
%Fcost = sqrt(((10*log10(Sfit)-10*log10(S)).^2./sum(10*log10(S).^2)) + 0.5*(Sfit-S).^2); 
%Fcost = sqrt(((lg(Sfit)-lg(S)).^2./trapz(lg(S).^2)) + w*(Sfit-S).^2); 
%Fcost = sqrt(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2);
%Fcost = sqrt(sqrt(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2));
%Fcost = sqrt((abs(lg(Sfit)-lg(S))./sum(lg(S))) + w*abs(Sfit-S));
%Fcost = sum(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2);
end