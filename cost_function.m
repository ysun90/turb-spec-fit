function Fcost = cost_function(F,S0,S,Sfit,method,wlinear,wcenter)
% Sfit - fitting function
% S = spectrum
% w - weight in linear scale

% Weight for linear part and center point
%global wlinear wcenter
% wlinear=0.5;
% wcenter = 0.2;

% Remove central point (513)
Weight = ones(length(F),1); Weight(513) = wcenter;

% Cost function
switch method
    case {'SL','SM','VLS','vector',} % For least squared vector
        Fcost = sqrt(  ( (1-wlinear)*(lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)) + wlinear*(Sfit-S).^2 ).*Weight  ); % default
        %Fcost = sqrt(  (lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)*1000/1025) + w*(Sfit-S).^2  ); % default
    case {'AL','AG','A'} % For absolute vector
        Fcost = sqrt( abs(lg(Sfit)-lg(S))/(sum(abs(lg(S)))) + wlinear*abs(Sfit-S));
    case {'local','multi','global','SG'} % for fmincon
        Fcost = ( (1-wlinear)*(lg(Sfit)-lg(S)).^2./(sum(lg(S).^2)) + wlinear*(Sfit-S).^2 ).*Weight;
        %Fcost = ( (1-wlinear)*(lg(sum(S0)*Sfit)-lg(S0)).^2./(sum(lg(S0).^2)) + wlinear*(Sfit-S).^2 ).*Weight;
        %Fcost = sum(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2); % default
        Fcost = sum(Fcost);
end
% Cost Function
%Fcost = sqrt(((10*log10(Sfit)-10*log10(S)).^2./sum(10*log10(S).^2)) + 0.5*(Sfit-S).^2); 
%Fcost = sqrt(((lg(Sfit)-lg(S)).^2./trapz(lg(S).^2)) + w*(Sfit-S).^2); 
%Fcost = sqrt(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2);
%Fcost = sqrt(sqrt(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2));
%Fcost = sqrt((abs(lg(Sfit)-lg(S))./sum(lg(S))) + w*abs(Sfit-S));
%Fcost = sum(((lg(Sfit)-lg(S)).^2./sum(lg(S).^2)) + w*(Sfit-S).^2);
end