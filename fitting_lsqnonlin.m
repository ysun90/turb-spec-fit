
function [X, Spec, Err, Turb] = fitting_lsqnonlin(si, model, method)
% si - Spectrum index
% model - 'GD': The Gaussian distribution
%         'GGD': The generalized Gaussian distribution
%         'VD': The Voigt distribution
%         'TD': FFT of the Taylor distribution
% method - 'SL': Squared Local
%          'AL': Absolute Local
%          'SG': Squared Global
%          'AG': Absolute Global
% X - model optimal parameters
% Err - error, BIC
% Spec - S, Sfit, each component
% Turb - Energy and width of the fitting components

% Number of frequency
global nfft

% Spectrum in kHz
[F,S0] = spectrum_original(si);
% Normalization
S=S0/trapz(F,S0);

%% Optimization
tic;
% Solver and set algorithm
OPTIONS = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

% The objective function
switch model
    case 'GD'
        FUN = @(x)cost_function_gaussian_model(x,F,S,method);
    case 'GGD'
        FUN = @(x)cost_function_ggaussian_model(x,F,S,method);
    case 'VD'
        FUN = @(x)cost_function_voigt_model(x,F,S,method);
    case 'TD'
        FUN = @(x)cost_function_taylor_model(x,F,S,method);
end

% Constraints
ampCS_lb = 0; ampCS_ub = inf; ampCS_x0 = max(S); LB(1) = ampCS_lb; UB(1) = ampCS_ub; X0(1) = ampCS_x0;    % ampCS     
ampLF_lb = 0; ampLF_ub = inf; ampLF_x0 = max(S); LB(4) = ampLF_lb; UB(4) = ampLF_ub; X0(4) = ampLF_x0;    % ampLF
ampBB_lb = 0; ampBB_ub = inf; ampBB_x0 = max(S); LB(7) = ampBB_lb; UB(7) = ampBB_ub; X0(7) = ampBB_x0;    % ampBB
muCS_lb =   -2; muCS_ub =   2; muCS_x0 = 0; LB(2) = muCS_lb; UB(2) = muCS_ub; X0(2) = muCS_x0;  % muCS
muLF_lb =  -15; muLF_ub =  15; muLF_x0 = 0; LB(5) = muLF_lb; UB(5) = muLF_ub; X0(5) = muLF_x0;  % muLF
muBB_lb = -500; muBB_ub = 500; muBB_x0 = 0; LB(8) = muBB_lb; UB(8) = muBB_ub; X0(8) = muBB_x0;  % muBB
sigmaCS_lb = 0; sigmaCS_ub =  3; sigmaCS_x0 =  1; LB(3) = sigmaCS_lb; UB(3) = sigmaCS_ub; X0(3) = sigmaCS_x0; % sigmaCS
sigmaLF_lb = 3; sigmaLF_ub = 20; sigmaLF_x0 = 10; LB(6) = sigmaLF_lb; UB(6) = sigmaLF_ub; X0(6) = sigmaLF_x0; % sigmaLF
switch model
    case 'GD'
        sigmaBB_lb = 20; sigmaBB_ub = 300; sigmaBB_x0 = 80; LB(9) = sigmaBB_lb; UB(9) = sigmaBB_ub; X0(9) = sigmaBB_x0; % sigmaBB       
    case 'GGD'
        alphaBB_lb = 20; alphaBB_ub = 300; alphaBB_x0 = 50; LB(9) = alphaBB_lb; UB(9) = alphaBB_ub; X0(9) = alphaBB_x0; % alphaBB
        betaBB_lb = 0.5; betaBB_ub = 8; betaBB_x0 = 1.05; LB(10) = betaBB_lb; UB(10) = betaBB_ub; X0(10) = betaBB_x0;   % betaBB    
    case 'VD'
        sigmaBBG_lb = 20; sigmaBBG_ub = 400; sigmaBBG_x0 = 80; LB(9) = sigmaBBG_lb; UB(9) = sigmaBBG_ub; X0(9) = sigmaBBG_x0;% sigmaBBG
        gammaBBL_lb = 1e-4; gammaBBL_ub = 50; gammaBBL_x0 = 2; LB(10) = gammaBBL_lb; UB(10) = gammaBBL_ub; X0(10) = gammaBBL_x0;% gammaBBL
    case 'TD'
        tauBB_lb = 0; tauBB_ub = 10; tauBB_x0 = 1; LB(9) = tauBB_lb; UB(9) = tauBB_ub; X0(9) = tauBB_x0;   % tauBB
        kappaBB_lb = 0.2; kappaBB_ub = 10; kappaBB_x0 = 1; LB(10) = kappaBB_lb; UB(10) = kappaBB_ub; X0(10) = kappaBB_x0;   % kappaBB
end
Noise_lb = 0; Noise_ub = max(S); Noise_x0 = min(S);
switch model
    case 'GD'
        LB(10) = Noise_lb; UB(10) = Noise_ub; X0(10) = Noise_x0;
    case {'GGD','VD','TD'}
        LB(11) = Noise_lb; UB(11) = Noise_ub; X0(11) = Noise_x0;
end

% Set problem
problem = createOptimProblem('lsqnonlin','objective',FUN,'x0',X0,'lb',LB,'ub',UB,'options',OPTIONS);
%problem = createOptimProblem('lsqnonlin','objective',FUN,'x0',X0,'lb',LB,'ub',UB,'options',optimset('OutputFcn',@PlotIterates));
% Optimize
switch method
    case {'SL','AL'}
        [x,normFUN,resFUN] = lsqnonlin(problem);
    case {'SG','AG'}
        % Global optimization by MultiStart
        ms =  MultiStart; 
        ms.StartPointsToRun = 'bounds';
        k = 20;
        x = run(ms,problem,k);
end
time = toc;
fprintf('Optimization time = %3.2f [s]\n',time);

%% Save the fitting parameters
% Components
CCS = gaussian_function(x(1),x(2),x(3),F); % Central Spike (CS)
CLF = gaussian_function(x(4),x(5),x(6),F); % Low Frequency (LF)
switch model  % Broadband (BB)
    case 'GD'
        CBB = gaussian_function(x(7),x(8),x(9),F);
    case 'GGD'
        CBB = ggaussian_function(x(7),x(8),x(9),x(10),F);
    case 'VD'
        CBB = x(7)*voigt(F,x(8),x(9),x(10));
    case 'TD'
        CBB = x(7)*taylorFFT(F,x(8),x(9),x(10)); 
end
CN = ones(nfft,1).*x(end); % Noise (N)
Sfit = CCS + CLF + CBB + CN;
Spec.F = F; Spec.S = S; Spec.Sfit = Sfit;
Spec.CCS = CCS; Spec.CLF = CLF; Spec.CBB = CBB; 
CBBcut = CBB; CBBcut(abs(F)<=sigmaLF_ub) = min(CBB(abs(F)<=sigmaLF_ub));
Spec.CBBcut = CBBcut;
Spec.CN = CN;
% Parameters
X.ampCS = x(1); X.muCS = x(2); X.sigmaCS = x(3);
X.ampLF = x(4); X.muLF = x(5); X.sigmaLF = x(6);
X.ampBB = x(7); X.muBB = x(8);
switch model
    case 'GD'
        X.sigmaBB = x(9); 
    case 'GGD'
        X.alphaBB = x(9); X.betaBB = x(10);
    case 'VD'
        X.sigmaBBG = x(9); X.gammaBBL = x(10); 
    case 'TD'
        X.tauBB = x(9); X.kappaBB = x(10); 
end
X.noise = x(end);
% Error
Err.res = Sfit-S; 
Err.norm2 = norm(Sfit-S); % residual sum of squares (RSS)
Err.BIC = 2*nfft*log(std(Sfit-S)) + length(x)*log(nfft);
Err.resFUN = resFUN;
Err.norm2FUN = sqrt(normFUN);
Err.BICFUN = 2*nfft*log(std(resFUN)) + length(x)*log(nfft);
% Turbulence properties
Turb.ES0 = trapz(F,S0);
fun_gaussian = @(f,amp,mu,sigma) amp*exp( -0.5*( ( ( f-mu )/sigma ).^2 ) );
Turb.ECS = trapz(F,CCS);Turb.ECSint = integral(@(f)fun_gaussian(f,x(1),x(2),x(3)),-500,500);
Turb.ELF = trapz(F,CLF);Turb.ELFint = integral(@(f)fun_gaussian(f,x(4),x(5),x(6)),-500,500);
Turb.EBB = trapz(F,CBB);Turb.EBBcut = trapz(F,CBBcut); Turb.EN = trapz(F,CN);
Turb.sdCS = x(3); %Turb.hwhmCS = hwhm(F,CCS)/2;
Turb.sdLF = x(6); %Turb.hwhmLF = hwhm(F,CLF)/2;
Turb.sdBB = sd_spec(F,CBB); %Turb.hwhmBB = hwhm(F,CBB);
% if strcmp(model,'VD')
%     fG = 2*x(9)*sqrt(2*log(2)); fL = 2*x(10);
%     Turb.hwhmBBD = ( 0.5346*fL + sqrt(0.2166*fL.^2 + fG.^2) )/2;
% end

%% Plot the fitting results
figure('Position',[18 128 1327 508],'Color','w');
% Linear scale
subplot(121);
hold on
plot(F,S,'c-','LineWidth',1.5);
plot(F,Sfit,'r-','LineWidth',3);
plot(F,CBB,'b-.','LineWidth',2);
plot(F,CLF,'m--','LineWidth',2);
plot(F,CCS,'k:','LineWidth',2);
plot(F,CN,'g-','LineWidth',2);
hold off
xlabel('Frequency [kHz]');
ylabel('Power spectrum [dB]');
title(['@',num2str(si),', #',num2str(numc),', Fx = ',num2str(Facqu,3), ' GHz']);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XMinorTick = 'on';
ax.XLim = [-450 450];
ax.YLim = [min(S) max(S)];
ax.YAxis.Exponent = 0;
lgd = legend('S','CS+LF+BB+N','BB','LF','CS','N');
lgd.FontSize = 12;
if ~strcmp(model,'GD')
    if strcmp(model,'GGD')
        t1 = text(ax.XLim(1)+10,max(Sfit)*0.95,sprintf(' \\alpha_{BB} = %4.3f,\n \\beta_{BB} = %4.3f,',x(9),x(10)));
    elseif strcmp(model,'VD')
        t1 = text(ax.XLim(1)+10,max(Sfit)*0.95,sprintf(' \\sigma_{BBG} = %4.3f,\n \\gamma_{BBL} = %4.3f,',x(9),x(10)));
    elseif strcmp(model,'TD')
        t1 = text(ax.XLim(1)+10,max(Sfit)*0.95,sprintf(' \\tau_{BB} = %4.3f,\n \\kappa_{BB} = %4.3f,',x(9),x(10)));
    end
t1.FontSize = 12;
t1.VerticalAlignment = 'top'; % cap
t1.HorizontalAlignment = 'left';
end
subplot(122);
hold on
plot(F,lg(S),'c-','LineWidth',1.5);
plot(F,lg(Sfit),'r-','LineWidth',3);
plot(F,lg(CBB),'b-.','LineWidth',2);
plot(F,lg(CLF),'m--','LineWidth',2);
plot(F,lg(CCS),'k:','LineWidth',2);
plot(F,lg(CN),'g-','LineWidth',2);
hold off
xlabel('Frequency [kHz]');
ylabel('Power spectrum');
title(['@',num2str(si),', \rho = ',num2str(r,2),', ',model,', ',method]);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XLim = [-450 450];
ax.XMinorTick = 'on';
ax.YLim = [min(lg(S))*1.1 max(lg(S))];
t2 = text(ax.XLim(1)+10,max(lg(S)),sprintf(' norm2 = %.3f,\n BIC = %d,\n norm2FUN = %.3f,\n BICFUN = %d,\n E_{S0} = %4.3f,\n E_{CS} = %4.3f,\n E_{LF} = %4.3f,\n E_{BB} = %4.3f,\n E_{BBcut} = %4.3f',...
            Err.norm2,round(Err.BIC),Err.norm2FUN,round(Err.BICFUN),...
            Turb.ES0,Turb.ECS,Turb.ELF,Turb.EBB,Turb.EBBcut));
t2.FontSize = 12;
t2.VerticalAlignment = 'top'; % cap
t2.HorizontalAlignment = 'left';
t3 = text(ax.XLim(2)-10,max(lg(S)),sprintf(' SD_{CS} = %4.3f,\n SD_{LF} = %4.3f,\n SD_{BB} = %4.3f',Turb.sdCS,Turb.sdLF,Turb.sdBB));
t3.FontSize = 12;
t3.VerticalAlignment = 'top'; % cap
t3.HorizontalAlignment = 'right';