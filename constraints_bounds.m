function [LB,UB] = constraints_bounds(S,model)
% [LB,UB] = constraints_bounds(S,model)
LB = zeros(1,11);
UB = zeros(1,11);

% Same constraints 
ampCS_lb = 0; ampCS_ub = inf; LB(1) = ampCS_lb; UB(1) = ampCS_ub; % ampCS 
ampLF_lb = 0; ampLF_ub = inf; LB(4) = ampLF_lb; UB(4) = ampLF_ub; % ampLF
ampBB_lb = 0; ampBB_ub = inf; LB(7) = ampBB_lb; UB(7) = ampBB_ub; % ampBB
muCS_lb =   -1; muCS_ub =   1;LB(2) = muCS_lb; UB(2) = muCS_ub; % muCS
muLF_lb =  -10; muLF_ub =  10;LB(5) = muLF_lb; UB(5) = muLF_ub;  % muLF
muBB_lb = -250; muBB_ub = 250;LB(8) = muBB_lb; UB(8) = muBB_ub; % muBB
sigmaCS_lb = 0; sigmaCS_ub =  2.5;LB(3) = sigmaCS_lb; UB(3) = sigmaCS_ub; % sigmaCS
sigmaLF_lb = 1; sigmaLF_ub = 20;LB(6) = sigmaLF_lb; UB(6) = sigmaLF_ub; % sigmaLF
Noise_lb = 0; Noise_ub = inf; LB(11) = Noise_lb; UB(11) = Noise_ub; % Noise

% Different constraints
switch model   
    case 'GGD'
        sigmaBB_lb = 10; sigmaBB_ub = 300;LB(9) = sigmaBB_lb; UB(9) = sigmaBB_ub; % alphaBB
        betaBB_lb = 0.5; betaBB_ub = 8; LB(10) = betaBB_lb; UB(10) = betaBB_ub; % betaBB    
    case {'VD','PVD','PVD2VD'}
        sigmaBBG_lb = 10; sigmaBBG_ub = 300; LB(9) = sigmaBBG_lb; UB(9) = sigmaBBG_ub; % sigmaBBG
        gammaBBL_lb = 0.001; gammaBBL_ub = 300; LB(10) = gammaBBL_lb; UB(10) = gammaBBL_ub; % gammaBBL
    case 'bisTD'
        k2DBB_lb = 0.01; k2DBB_ub = 100; LB(9) = k2DBB_lb; UB(9) = k2DBB_ub; % nuBB
        tauBB_lb = 0.01; tauBB_ub = 100; LB(10) = tauBB_lb; UB(10) = tauBB_ub; % tauBB        
    case 'TD'
        kappaBB_lb = 0.1; kappaBB_ub = 5; LB(9) = kappaBB_lb; UB(9) = kappaBB_ub; % kappaBB
        tauBB_lb = 0.001; tauBB_ub = 10; LB(10) = tauBB_lb; UB(10) = tauBB_ub; % tauBB        
end

end

%% backup
% ampCS_lb = 0; ampCS_ub = inf; ampCS_x0 = max(S); LB(1) = ampCS_lb; UB(1) = ampCS_ub; X0(1) = ampCS_x0;    % ampCS     
% ampLF_lb = 0; ampLF_ub = inf; ampLF_x0 = max(S); LB(4) = ampLF_lb; UB(4) = ampLF_ub; X0(4) = ampLF_x0;    % ampLF
% ampBB_lb = 0; ampBB_ub = inf; ampBB_x0 = max(S); LB(7) = ampBB_lb; UB(7) = ampBB_ub; X0(7) = ampBB_x0;    % ampBB
% muCS_lb =   -2; muCS_ub =   2; muCS_x0 = 0; LB(2) = muCS_lb; UB(2) = muCS_ub; X0(2) = muCS_x0;  % muCS
% muLF_lb =  -15; muLF_ub =  15; muLF_x0 = 0; LB(5) = muLF_lb; UB(5) = muLF_ub; X0(5) = muLF_x0;  % muLF
% muBB_lb = -300; muBB_ub = 300; muBB_x0 = 0; LB(8) = muBB_lb; UB(8) = muBB_ub; X0(8) = muBB_x0;  % muBB
% sigmaCS_lb = 0; sigmaCS_ub =  3; sigmaCS_x0 =  1; LB(3) = sigmaCS_lb; UB(3) = sigmaCS_ub; X0(3) = sigmaCS_x0; % sigmaCS
% sigmaLF_lb = 3; sigmaLF_ub = 20; sigmaLF_x0 = 10; LB(6) = sigmaLF_lb; UB(6) = sigmaLF_ub; X0(6) = sigmaLF_x0; % sigmaLF
% switch model
%     case 'GD'
%         sigmaBB_lb = 20; sigmaBB_ub = 300; sigmaBB_x0 = 80; LB(9) = sigmaBB_lb; UB(9) = sigmaBB_ub; X0(9) = sigmaBB_x0; % sigmaBB       
%     case 'GGD'
%         alphaBB_lb = 20; alphaBB_ub = 300; alphaBB_x0 = 50; LB(9) = alphaBB_lb; UB(9) = alphaBB_ub; X0(9) = alphaBB_x0; % alphaBB
%         betaBB_lb = 0.5; betaBB_ub = 8; betaBB_x0 = 1.05; LB(10) = betaBB_lb; UB(10) = betaBB_ub; X0(10) = betaBB_x0;   % betaBB    
%     case 'VD'
%         sigmaBBG_lb = 30; sigmaBBG_ub = 400; sigmaBBG_x0 = 80; LB(9) = sigmaBBG_lb; UB(9) = sigmaBBG_ub; X0(9) = sigmaBBG_x0;% sigmaBBG
%         gammaBBL_lb = 1e-4; gammaBBL_ub = 50; gammaBBL_x0 = 2; LB(10) = gammaBBL_lb; UB(10) = gammaBBL_ub; X0(10) = gammaBBL_x0;% gammaBBL
%     case 'TD'
%         tauBB_lb = 0; tauBB_ub = 10; tauBB_x0 = 1; LB(9) = tauBB_lb; UB(9) = tauBB_ub; X0(9) = tauBB_x0;   % tauBB
%         kappaBB_lb = 0.2; kappaBB_ub = 8; kappaBB_x0 = 1; LB(10) = kappaBB_lb; UB(10) = kappaBB_ub; X0(10) = kappaBB_x0;   % kappaBB
% end
% Noise_lb = 0; Noise_ub = max(S); Noise_x0 = min(S);
% switch model
%     case 'GD'
%         LB(10) = Noise_lb; UB(10) = Noise_ub; X0(10) = Noise_x0;
%     case {'GGD','VD','TD'}
%         LB(11) = Noise_lb; UB(11) = Noise_ub; X0(11) = Noise_x0;
% end