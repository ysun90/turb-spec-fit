
function [fval,X, Spec, Err, Turb] = fitting(si, Niv, model, method, wlinear, wcenter)
% [X, Spec, Err, Turb] = fitting(si, Niv, model, method, wlinear, wcenter)
%
% Inputs:
% si - Spectrum index
% model - 'GD','GGD','VD','TD'
% method - 'SL' : Vector Least Squares (lsqnonlin, 'trust-region-reflective')
%          'SG': Scalar Absolute Global (fmincon, 'interior-point', GlobalSearch)
%
% Outputs:
% X - model optimal parameters
% Err - error, BIC
% Spec - S, Sfit, each component
% Turb - Energy and width of the fitting components

% Defaults for the inputs
if nargin<6
   wcenter = 0.2;
   if nargin<5
       wlinear = 0.5;
       if nargin<4
           method = 'SL';
           if nargin<3
               model='GGD';
               if nargin<2
                   Niv = 1;
               end
           end
       end
   end
end

% Load frequency of spectrum and database of index
loadata;

% Spectrum in kHz
S0 = ind2spec(si,database_index);
S=S0/trapz(F,S0); % Normalization

% Check the length of F and S
if length(F)~=length(S)
    error('Lengthes of F and S are not the same!');
end

%% Optimization
tic;
% Solver and set options
if Niv==inf
    OPTIONS = optimoptions(@fmincon,'Algorithm','interior-point',...
         'Display','iter-detailed','Diagnostics','on','FunValCheck','on','PlotFcn',@optimplotx);
else 
    OPTIONS = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
         'Display','iter-detailed','Diagnostics','on','FunValCheck','on','PlotFcn',@optimplotx);
end

% The objective function
FUN = @(x)objective_function(x,F,S,model,method,wlinear,wcenter);

% Constraints and initial values
[LB,UB] = parameter_constraints(S,model);
X0 = initial_values(F,S,model);
X0Matrix = repmat(X0,5,1);
X0Matrix(2,10) = X0(10)*4;
X0Matrix(3,10) = X0(10)/4;
X0Matrix(4,10) = X0(10)*2;
X0Matrix(5,10) = X0(10)/2;

% Optimize
switch method
    case 'SL'
        problem = createOptimProblem('lsqnonlin','objective',FUN,'x0',X0,'lb',LB,'ub',UB,'options',OPTIONS);
        ms =  MultiStart;       
        tpoints = CustomStartPointSet(X0Matrix(1:2,:));
        %rpts = RandomStartPointSet('NumStartPoints',2);
        %allpts = {tpoints,rpts};
        [x,fval] = run(ms,problem,tpoints);
        %[x,fval] = lsqnonlin(problem);
    case 'SG'
        problem = createOptimProblem('fmincon','objective',FUN,'x0',X0,'lb',LB,'ub',UB,'options',OPTIONS);
        gs = GlobalSearch; 
        [x,fval] = run(gs,problem);
end
time = toc;
fprintf('Optimization time = %3.2f [s]\n',time);
disp(['SSR = ',num2str(fval)]);

%% Save the fitting parameters
% Components
CCS = gaussian_function(x(1),x(2),x(3),F); % Central Spike (CS)
CLF = gaussian_function(x(4),x(5),x(6),F); % Low Frequency (LF)
switch model  % Broadband (BB)
    case 'GD'
        CBB = gaussian_function(x(7),x(8),x(9),F);
    case 'GGD'
        CBB = ggaussian_function(x(7),x(8),x(9),x(10),F);
    case 'PVD'
        CBB = pseudo_voigt_function(F,x(7),x(8),x(9),x(10));
    case 'VD'
        CBB = x(7)*voigt(F,x(8),x(9),x(10));
    case 'TD'
        CBB = x(7)*taylorFFT(F,x(8),x(9),x(10)); 
end
CN = ones(length(F),1).*x(end); % Noise (N)
Sfit = CCS + CLF + CBB + CN;
Spec.F = F; Spec.S = S; Spec.Sfit = Sfit;
Spec.CCS = CCS; Spec.CLF = CLF; Spec.CBB = CBB; 
%CBBcut = CBB; CBBcut(abs(F)<=sigmaLF_ub) = min(CBB(abs(F)<=sigmaLF_ub));
%Spec.CBBcut = CBBcut;
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
Err.resnorm = fval; 
Err.RSS = sum((Sfit-S).^2); % residual sum of squares (RSS)
Err.BIC = 2*length(F)*log(std(Sfit-S)) + length(x)*log(length(F));
% Turbulence properties
Turb.ES0 = trapz(F,S0);
fun_gaussian = @(f,amp,mu,sigma) amp*exp( -0.5*( ( ( f-mu )/sigma ).^2 ) );
Turb.ECS = trapz(F,CCS);Turb.ECSint = integral(@(f)fun_gaussian(f,x(1),x(2),x(3)),-500,500);
Turb.ELF = trapz(F,CLF);Turb.ELFint = integral(@(f)fun_gaussian(f,x(4),x(5),x(6)),-500,500);
Turb.EBB = trapz(F,CBB);Turb.EN = trapz(F,CN);%Turb.EBBcut = trapz(F,CBBcut); 
Turb.sdCS = x(3); Turb.hwhmCS = hwhm(F,CCS)/2;
Turb.sdLF = x(6); Turb.hwhmLF = hwhm(F,CLF)/2;
Turb.sdBB = sd_spec(F,CBB); Turb.hwhmBB = hwhm(F,CBB);
if strcmp(model,'VD')
    fG = 2*x(9)*sqrt(2*log(2)); fL = 2*x(10);
    Turb.hwhmBBD = ( 0.5346*fL + sqrt(0.2166*fL.^2 + fG.^2) )/2;
end

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
%title(['@',num2str(si),', #',num2str(numc),', F = ',num2str(Facqu,3), 'GHz']);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XMinorTick = 'on';
ax.XLim = [-450 450];
ax.YLim = [min(Sfit) max(Sfit)];
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
ylabel('Power spectrum [dB]');
%title(['@',num2str(si),', \rho = ',num2str(r,2),', R = ',num2str(R,2),'m, ',model]);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XLim = [-450 450];
ax.XMinorTick = 'on';
ax.YLim = [min(lg(S))*1.1 max(lg(S))];
% t2 = text(ax.XLim(1)+10,max(lg(S)),...
%     sprintf(' RSS = %4.3e,\n BIC = %4.3e,\n',...
%             ' E_{S0} = %4.3f,\n E_{CS} = %4.3f,\n E_{LF} = %4.3f,\n',...
%             ' E_{BB} = %4.3f,\n',...
%             Err.RSS,Err.BIC,Turb.ES0,Turb.ECS,Turb.ELF,Turb.EBB));
% t2.FontSize = 12;
% t2.VerticalAlignment = 'top'; % cap
% t2.HorizontalAlignment = 'left';
% t3 = text(ax.XLim(2)-10,max(lg(S)),...
%     sprintf(' SD_{CS} = %4.3e,\n  HW_{CS} = %4.3e,\n',...
%             ' SD_{LF} = %4.3f,\n HW_{LF} = %4.3f,\n SD_{BB} = %4.3f,\n HW_{BB} = %4.3f',...
%             Turb.sdCS,Turb.hwhmCS,Turb.sdLF,Turb.hwhmLF,Turb.sdBB,Turb.hwhmBB));
% t3.FontSize = 12;
% t3.VerticalAlignment = 'top'; % cap
% t3.HorizontalAlignment = 'right';