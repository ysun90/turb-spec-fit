    
function [ERROR,X,XI,TURB,SPEC] = fitmincon(INDEX,MODEL,METHOD,WL,W0,RESULTS)
%*************************************************************************
%FITMINCON fits a frequency spectrum by fmincon slover
% [ERROR,X,X0,TURB,SPEC] = fitmincon(INDEX,MODEL,METHOD,WL,W0,RESULTS)
%*************************************************************************
%INPUTS:
% INDEX - Spectrum index [1,358874]
% MODEL - 'GD','GGD','VD','PVD','PVD2VD','TD'
% METHOD - 'local','multi','global'
% WL - Weight for the linear scale of the spectrum [0,1]
% W0 - Weight for the zero frequency component[0,1]
% RESULTS - 'plot','plotx','ploty'
%*************************************************************************
%OUTPUTS:
% ERROR - FVAL, BIC
% X - Optimal solution
% X0 - Estimated Initial guess
% TURB - Turblence properties
% SPEC - Spectrum and components
%*************************************************************************

% Defaults for the inputs
if nargin<6
   RESULTS = [];
   if nargin<5
       W0 = 1;
       if nargin<4
           WL = 0.5;
           if nargin<3
               METHOD = 'local';
               if nargin<2
                   MODEL = 'GGD';
               end
           end
       end
   end
end

%% Fitting
% Spectrum
nfft = 1025; % Number of frequency
loadata; % Load frequency of spectrum and database of index
S0 = ind2spec(INDEX,database_index); % Spectrum in kHz
S=S0/sum(S0); % Normalization

% Solver and optimization options
OPTIONS = optimoptions(@fmincon); % fmincon for constraints and GlobalSearch
OPTIONS.Algorithm = 'interior-point'; % 'trust-region-reflective', 'sqp','active-set'
if ~isempty(RESULTS)
    OPTIONS.FunValCheck = 'on'; % 'off'
    OPTIONS.MaxIterations = 500; % Default 1000
    OPTIONS.Diagnostics = 'on';
    OPTIONS.Display = 'iter-detailed';
    if strcmp(RESULTS,'plotx');OPTIONS.PlotFcn = @optimplotx;end % plot X
    if strcmp(RESULTS,'ploty');OPTIONS.PlotFcn = @optimplotfval;end % plot FVAL
end

% The objective function
FUN = @(x)objective_function(x,F,S0,S,MODEL,METHOD,WL,W0);

% Constraints
[LB,UB] = constraints_bounds(S,MODEL);
if strcmpi(MODEL,'TD')|| strcmpi(MODEL,'bisTD')
    A = []; b=[];
else
    eta=1.5; %sigmaLF>1.5*sigmCS, sigmaBB>1.5*sigmaLF
    A = [0,0,eta,0,0,-1,0,0,0,0,0; 0,0,0,0,0,eta,0,0,-1,0,0];
    b = [0;0];
end

% Initial guess
X0 = initial_values(F,S,MODEL);

% Optimization
tic;
problem = createOptimProblem('fmincon','objective',FUN,...
    'x0', X0,'Aineq',A,'bineq',b,'lb', LB, 'ub', UB,'options', OPTIONS);
switch METHOD
    case 'local' % Local optimization
        [x,fval] = fmincon(problem); 
        
    case 'multi' % Multiple initial guesses
        Niv=3;
        X0Matrix = Multi_initial_guess(X0,MODEL,F,S);
        ms =  MultiStart;            
        tpoints = CustomStartPointSet(X0Matrix(1:Niv,:));
        [x,fval] = run(ms,problem,tpoints);
    case 'global' % Global optimization
        gs =  GlobalSearch;
        [x,fval] = run(gs,problem);
end
disp(['Optimization time: ',num2str(toc)]);
disp(['FVAL at optimum: ',num2str(fval)]);

% Components
CCS = x(1)*gaussian_function(F,x(2),x(3)); % Central Spike (CS)
CLF = x(4)*gaussian_function(F,x(5),x(6)); % Low Frequency (LF)
switch MODEL  % Broadband (BB)
    case 'GD'
        CBB = gaussian_function(x(7),x(8),x(9),F);
    case 'GGD'
        CBB = x(7)*ggaussian_function(F,x(8),x(9),x(10));
    case {'VD','PVD2VD'}
        CBB = x(7)*voigt_function(F,x(8),x(9),x(10));
    case 'PVD'
        CBB = x(7)*pvoigt(F,x(8),x(9),x(10));
    case 'TD'
        CBB = x(7)*taylorFFT(F,x(8),x(9),x(10)); 
    case 'bisTD'
        t=(-512:512)';
        CBB = x(7)*taylorbisFFT(t,x(8),x(9),x(10));
end
CN = ones(1025,1).*x(end); % Noise (N)
Sfit = CCS + CLF + CBB + CN;

%% Save parameters
% Optimum and initial guess
X.ampCS = x(1); X.muCS = x(2); X.sigmaCS = x(3);
X.ampLF = x(4); X.muLF = x(5); X.sigmaLF = x(6);
X.ampBB = x(7); X.muBB = x(8);
XI.ampCS = X0(1); XI.muCS = X0(2); XI.sigmaCS = X0(3);
XI.ampLF = X0(4); XI.muLF = X0(5); XI.sigmaLF = X0(6);
XI.ampBB = X0(7); XI.muBB = X0(8);
switch MODEL
    case 'GGD'
        X.sigmaBB = x(9); X.betaBB = x(10);
        XI.sigmaBB = X0(9); XI.betaBB = X0(10);
    case {'VD','PVD','PVD2VD'}
        X.sigmaBBG = x(9); X.gammaBBL = x(10);
        XI.sigmaBBG = X0(9); XI.gammaBBL = X0(10); 
    case 'TD'
        X.kappaBB = x(9); X.tauBB = x(10); 
        XI.kappaBB = X0(9); XI.tauBB = X0(10);
    case 'bisTD'
        X.k2DBB = x(9); X.tauBB = x(10);
        XI.k2DBB = X0(9); XI.tauBB = X0(10);       
end
X.noise = x(11); XI.noise = X0(11);
% Error
ERROR.FVAL = fval;
ERROR.BIC = 2*nfft*log(std(Sfit-S)) + length(x)*log(nfft);
% Components
SPEC.F = F; SPEC.S = S; SPEC.Sfit = Sfit;
SPEC.CS = CCS; SPEC.LF = CLF; SPEC.BB = CBB; SPEC.N = CN;
% Turbulence properties
TURB.ES0 = sum(S0); %fun_gaussian = @(f,amp,mu,sigma) amp*exp( -0.5*( ( ( f-mu )/sigma ).^2 ) );
TURB.ECS = sum(CCS);%TURB.ECSint = integral(@(f)fun_gaussian(f,x(1),x(2),x(3)),-500,500);
TURB.ELF = sum(CLF);%TURB.ELFint = integral(@(f)fun_gaussian(f,x(4),x(5),x(6)),-500,500);
TURB.EBB = sum(CBB);%Turb.EBBcut = trapz(F,CBBcut); 
TURB.EN = sum(CN);
TURB.SDCS = x(3);
TURB.SDLF = x(6);
TURB.SDBB = std_pdf_norm(F,CBB);
% Turb.sdCS = x(3); Turb.hwhmCS = fwhm(F,CCS)/2;
% Turb.sdLF = x(6); Turb.hwhmLF = fwhm(F,CLF)/2;
% Turb.sdBB = sd_spec(F,CBB); Turb.hwhmBB = fwhm(F,CBB)/2;
% if strcmp(model,'VD')
%     fG = 2*x(9)*sqrt(2*log(2)); fL = 2*x(10);
%     Turb.hwhmBBD = ( 0.5346*fL + sqrt(0.2166*fL.^2 + fG.^2) )/2;
% end

%% Plot the fitting results
if ~isempty(RESULTS)
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
ylabel('Power spectrum');
title('Linear Scale');
%title(['@',num2str(si),', #',num2str(numc),', Fx = ',num2str(Facqu,3), ' GHz']);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XMinorTick = 'on';
ax.XLim = [-35 35];
ax.YLim = [min(S) max(S)];
ax.YAxis.Exponent = 0;
lgd = legend('Spectrum','Fitting','BB','LF','CS','N');
lgd.FontSize = 12;
switch MODEL
    case 'GGD'
        t1 = text(ax.XLim(1)+5,max(Sfit)*0.65,sprintf(' \\sigma_{BB} = %4.3f,\n \\beta_{BB} = %4.3f,',x(9),x(10)));
    case {'VD','PVD','PVD2VD'}
        t1 = text(ax.XLim(1)+5,max(Sfit)*0.65,sprintf(' \\sigma_{BBG} = %4.3f,\n \\gamma_{BBL} = %4.3f,',x(9),x(10)));
    case 'TD'
        t1 = text(ax.XLim(1)+5,max(Sfit)*0.65,sprintf(' \\kappa_{BB} = %4.3f,\n \\tau_{BB} = %4.3f,',x(9),x(10)));
    case 'bisTD'
        t1 = text(ax.XLim(1)+5,max(Sfit)*0.65,...
            sprintf(' k2D_{BB} = %4.3f,\n \\tau_{BB} = %4.3f',x(9),x(10)));
end
t1.FontSize = 12;
t1.VerticalAlignment = 'bottom'; % cap
t1.HorizontalAlignment = 'left';
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
ylabel('Power spectrum [dB]');
title(['Spectrum : ',num2str(INDEX)]);
%title(['@',num2str(si),', \rho = ',num2str(r,2),', ',model,', ',method]);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XLim = [-450 450];
ax.XMinorTick = 'on';
ax.YLim = [min(lg(S))*1.1 max(lg(S))];
t2 = text(ax.XLim(1)+10,max(lg(S)),sprintf(' FVAL = %.5f,\n BIC = %d,\n E_{S0} = %4.3f,\n E_{CS} = %4.3f,\n E_{LF} = %4.3f,\n E_{BB} = %4.3f',...
            ERROR.FVAL,round(ERROR.BIC),TURB.ES0,TURB.ECS,TURB.ELF,TURB.EBB));
t2.FontSize = 12;
t2.VerticalAlignment = 'top'; % cap
t2.HorizontalAlignment = 'left';
t3 = text(ax.XLim(2)-10,max(lg(S)),sprintf(' \\sigma_{CS} = %4.3f,\n \\sigma_{LF} = %4.3f,\n \\sigma_{BB} = %4.3f',TURB.SDCS,TURB.SDLF,TURB.SDBB));
t3.FontSize = 12;
t3.VerticalAlignment = 'top'; % cap
t3.HorizontalAlignment = 'right';

end
