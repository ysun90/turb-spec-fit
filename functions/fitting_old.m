
function [x,resnorm,FBB,FLF,FCS,FN] = fitting(si,method)

global nfft

load('database_r.mat'); % Load position database

[f,S] = spectrum_norm(si); % Get normalization spectrum

numc = ind2tab(si); % Get shot number

if method=='3g'   
%% Fit and get components by GD + GD + GD
%      x(1)     x(2)    x(3)     x(4)     x(5)     x(6)      x(7)        x(8)       x(9)       x(10)  
%    amp_CS , mu_CS, sigma_CS, amp_LF,   mu_LF,  sigma_LF,  amp_BB,     mu_BB    sigma_BB,     Noise
lb = [0,       -2,      0,      0,        -20,     3,         0,        -500       30,           0];
ub = [inf,      2,      3,      inf,       20,     30,        inf,       500,     300,      max(S)];
x0 = [max(S),   0,      1,      max(S),     0,     10,        max(S),      0,      80,       min(S)];
%options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','TolFun',1e-9);
[x,resnorm] = lsqnonlin(@myfun_fitting_3g,x0,lb,ub,[],f,S);
% Get fitting components
FCS = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ); % Central Spike (CS) 
FLF = x(4)*exp( -0.5*( ( ( f-x(5) )/x(6) ).^2 ) ); % Low Frequency (LF)
FBB = x(7)*exp( -0.5*( ( ( f-x(8) )/x(9) ).^2 ) ); % Broadband (BB)
FN = ones(nfft,1).*x(end); % Noise (N)
% Full fitting curve
F = FCS + FLF + FBB + FN;
% Plot the spectrum and fitting
figure('Position',[18 128 1327 508],'Color','w');
subplot(121);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('Spectrum','Fitting');
lgd.FontSize = 12;
% t = text(-490,max(10*log10(S)),...
%     sprintf(' A_{BB} = %3.1e,\n A_{LF} = %3.1e,\n A_{CS} = %3.1e, \n N = %3.1e, \n Res = %4.3f'...
%     ,x(8),x(4),x(1),x(12),resnorm));
% t.FontSize = 14;
% t.VerticalAlignment = 'top'; % cap
% t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]');
ylabel('Power spectrum [dB]');
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off
subplot(122);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(FBB),'r-','LineWidth',1.5);
plot(f,10*log10(FLF),'m--','LineWidth',1.5);
plot(f,10*log10(FCS),'g:','LineWidth',1.5);
plot(f,10*log10(FN),'c-.','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -400:200:400;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('S','BB','LF','CS','N','CS+LF+BB+N');
lgd.FontSize = 16;
% t = text(-490,max(10*log10(S)),...
%     sprintf(' \\sigma_{BBG} = %3.1f, \\gamma_{BBL} = %3.1f\n \\alpha_{LF} = %3.1f, \\beta_{LF} = %3.1f\n \\sigma_{CS} = %3.1f'...
%     ,x(10),x(11),x(6),x(7),x(3)));
% t.FontSize = 14;
% t.VerticalAlignment = 'top'; % cap
% t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]','FontSize',20);
ylabel('Power spectrum [dB]','FontSize',20);
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off  
end

if method=='gg'   
%% Fit and get components by GD + GD + GGD
% %    amp_CS , mu_CS, sigma_CS, amp_LF,  mu_LF, sigma_LF,  amp_BB,   mu_BB   sigma_BB,     Noise
% lb = [0,       -2,      0,      0,      -20,     3,          0,     -500      30,           0];
% ub = [max(S),   2,      3,      inf,     20,     30,        inf,      500,    300,      max(S)];
% x0 = [max(S),   0,      1,      max(S),   0,     10,        max(S),    0,     80,       min(S)];
%      x(1)     x(2)    x(3)     x(4)    x(5)    x(6)      x(7)    x(8)      x(9)     x(10)     x(11) 
%    amp_CS , mu_CS, sigma_CS, amp_LF,  mu_LF, sigma_LF, amp_BB,  mu_BB,  alpha_BB, beta_BB,  Noise
lb = [0,       -2,      0,      0,      -20,     3,         0,     -300,      30,       1.5,       0];
ub = [inf,      2,      3,      inf,     20,     30,       inf,      300,     300,      8,     max(S)];
x0 = [max(S),   0,      1,      max(S),   0,     10,      max(S),     0,      50,       1,      min(S)];
% lb = [min(S),  -2,      0,      min(S),   -20,     3,       0,     min(S),    -400,     30,     1e-4,       0];
% ub = [max(S),   2,      3,      max(S),    20,     20,      8,     max(S),     400,    450,     200,    max(S)];
% x0 = [max(S),   0,      1,      max(S),     0,     5,       2,     max(S),     0,      100,     100,    min(S)];
%options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','TolFun',1e-9);
[x,resnorm] = lsqnonlin(@myfun_fitting_4g,x0,lb,ub,[],f,S);
% Get fitting components
FCS = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ); % Central Spike (CS) 
FLF = x(4)*exp( -0.5*( ( ( f-x(5) )/x(6) ).^2 ) ); % Low Frequency (LF)
FBB = x(7)*exp( -( ( abs( f-x(8) )./(x(9).*sqrt(gamma(1/x(10))/gamma(3/x(10)))) ).^x(10) ) ); % Broadband (BB)
FN = ones(nfft,1).*x(end); % Noise (N)
% Full fitting curve
F = FCS + FLF + FBB + FN;
% Plot the spectrum and fitting
figure('Position',[18 128 1327 508],'Color','w');
subplot(121);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('Spectrum','Fitting');
lgd.FontSize = 12;
t = text(-490,max(10*log10(S)),...
    sprintf('\\alpha_{BB} = %3.1f, \\beta_{LF} = %3.1f',x(9),x(10)));
t.FontSize = 18;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]');
ylabel('Power spectrum [dB]');
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off
subplot(122);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(FBB),'r-','LineWidth',1.5);
plot(f,10*log10(FLF),'m--','LineWidth',1.5);
plot(f,10*log10(FCS),'g:','LineWidth',1.5);
plot(f,10*log10(FN),'c-.','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -400:200:400;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('S','BB','LF','CS','N','CS+LF+BB+N');
lgd.FontSize = 16;
t = text(-490,max(10*log10(S)),...
    sprintf('\\alpha_{BB} = %3.1f, \n \\beta_{LF} = %3.1f',x(9),x(10)));
t.FontSize = 18;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]','FontSize',20);
ylabel('Power spectrum [dB]','FontSize',20);
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off  
end

if method=='v'   
%% Fit and get components by GD + GD + VD
%      x(1)     x(2)    x(3)     x(4)    x(5)     x(6)      x(7)    x(8)      x(9)     x(10)     x(11)    
%    amp_CS , mu_CS, sigma_CS, amp_LF,  mu_LF,  sigma_LF,   amp_BB,  mu_BB,  sigma_BBG, sigma_BBL,  Noise
lb = [0,       -2,      0,      0,      -20,      3,          0,     -300,     30,       1e-4,       0];
ub = [max(S),   2,      3,      inf,     20,      30,        inf,      300,     400,      50,     max(S)];
x0 = [max(S),   0,      1,      max(S),   0,      10,        max(S),     0,       80,      2,      min(S)];
% lb = [min(S),  -2,      0,      min(S),   -20,     3,       0,     min(S),    -400,     30,     1e-4,       0];
% ub = [max(S),   2,      3,      max(S),    20,     20,      8,     max(S),     400,    450,     200,    max(S)];
% x0 = [max(S),   0,      1,      max(S),     0,     5,       2,     max(S),     0,      100,     100,    min(S)];
%options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','TolFun',1e-9);
[x,resnorm] = lsqnonlin(@myfun_fitting_voigt_GD,x0,lb,ub,[],f,S);
% Get fitting components
FCS = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ); % Central Spike (CS) 
FLF = x(4)*exp( -0.5*( ( ( f-x(5) )/x(6) ).^2 ) ); % Low Frequency (LF)
FBB = x(7)*voigt(f,x(8),x(9),x(10)); % Broadband (BB)
FN = ones(nfft,1).*x(end); % Noise (N)
% Full fitting curve
F = FCS + FLF + FBB + FN;
% Plot the spectrum and fitting
figure('Position',[18 128 1327 508],'Color','w');
subplot(121);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('Spectrum','Fitting');
lgd.FontSize = 12;
t = text(-490,max(10*log10(S)),...
    sprintf(' A_{BB} = %3.1e,\n A_{LF} = %3.1e,\n A_{CS} = %3.1e, \n N = %3.1e, \n Res = %4.3f'...
    ,x(7),x(4),x(1),x(end),resnorm));
t.FontSize = 14;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]');
ylabel('Power spectrum [dB]');
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off
subplot(122);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(FBB),'r-','LineWidth',1.5);
plot(f,10*log10(FLF),'m--','LineWidth',1.5);
plot(f,10*log10(FCS),'g:','LineWidth',1.5);
plot(f,10*log10(FN),'c-.','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -400:200:400;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('S','BB','LF','CS','N','CS+LF+BB+N');
lgd.FontSize = 16;
t = text(-490,max(10*log10(S)),...
    sprintf(' \\sigma_{BBG} = %3.1f, \\gamma_{BBL} = %3.1f\n \\sigma_{LF} = %3.1f, \\beta_{LF} = %3.1f\n \\sigma_{CS} = %3.1f'...
    ,x(9),x(10),x(6),x(7),x(3)));
t.FontSize = 14;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]','FontSize',20);
ylabel('Power spectrum [dB]','FontSize',20);
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off  
end


if method=='t'   
%% Fit and get components by GD + GD + TD
%      x(1)     x(2)    x(3)     x(4)    x(5)    x(6)      x(7)    x(8)      x(9)     x(10)     x(11)  
%    amp_CS , mu_CS, sigma_CS, amp_LF,  mu_LF, sigma_LF,  amp_BB,   mu_BB     tL_BB,    K_BB,     Noise
lb = [0,       -2,      0,      0,      -20,     3,          0,     -500      0,        0.1,        0];
ub = [max(S),   2,      3,      inf,     20,     30,        inf,      500,     10,      10,     max(S)];
x0 = [max(S),   0,      1,      max(S),   0,     10,        max(S),     0,     1,       1,      min(S)];
[x,resnorm] = lsqnonlin(@myfun_fitting_taylor_GD,x0,lb,ub,[],f,S);
% Get fitting components
FCS = x(1)*exp( -0.5*( ( ( f-x(2) )/x(3) ).^2 ) ); % Central Spike (CS) 
FLF = x(4)*exp( -0.5*( ( ( f-x(5) )/x(6) ).^2 ) ); % Low Frequency (LF)
FBB = x(7)*taylorFFT(f,x(8),x(9),x(10)); % Broadband (BB)
FN = ones(nfft,1).*x(11); % Noise (N)
% Energy of BB and LF
EBB = trapz(f,FBB);
ELF = trapz(f,FLF);
% Full fitting curve
F = FCS + FLF + FBB + FN;
% Plot the spectrum and fitting
figure('Position',[18 128 1327 508],'Color','w');
subplot(121);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -500:250:500;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('Spectrum','Fitting');
lgd.FontSize = 12;
t = text(-490,max(10*log10(S)),...
    sprintf(' A_{BB} = %3.1e,\n A_{LF} = %3.1e,\n A_{CS} = %3.1e, \n N = %3.1e, \n Res = %4.3f'...
    ,x(7),x(4),x(1),x(11),resnorm));
t.FontSize = 14;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]');
ylabel('Power spectrum [dB]');
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off
subplot(122);
grid on
grid minor
hold on
plot(f,10*log10(S),'b-','LineWidth',1.5);
plot(f,10*log10(FBB),'r-','LineWidth',1.5);
plot(f,10*log10(FLF),'m--','LineWidth',1.5);
plot(f,10*log10(FCS),'g:','LineWidth',1.5);
plot(f,10*log10(FN),'c-.','LineWidth',1.5);
plot(f,10*log10(F),'k-','LineWidth',1.5);
ax = gca;
ax.Box = 'on';
ax.FontSize = 20;
ax.XLim = [-500 500];
ax.XTick = -400:200:400;
ax.YLim = [min(10*log10(S))*1.1 max(10*log10(S))];
lgd = legend('S','BB','LF','CS','N','CS+LF+BB+N');
lgd.FontSize = 16;
t = text(-490,max(10*log10(S)),...
    sprintf(' tL_{BB} = %3.1f, K_{BB} = %3.1f\n \\sigma_{LF} = %3.1f,\\sigma_{CS} = %3.1f\n E_{BB} = %4.3f\n E_{LF} = %4.3f',...
    x(9),x(10),x(6),x(3),EBB,ELF));
t.FontSize = 14;
t.VerticalAlignment = 'top'; % cap
t.HorizontalAlignment = 'left';
xlabel('Frequency [kHz]','FontSize',20);
ylabel('Power spectrum [dB]','FontSize',20);
title(['@',num2str(si),', #',num2str(numc),', \rho = ',num2str(database_r.r(si),3)]);
hold off  
end

end