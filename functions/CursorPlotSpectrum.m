
function  CursorPlotSpectrum(p1,p2,indexing,numVar)

global ind method
ind = indexing;
method = numVar;
figure;set(gcf,'Pos',[16   211   560   420]);
plot(p1(indexing),p2(indexing),'b.');

dcm = datacursormode(gcf);
datacursormode toggle
set(dcm,'updatefcn',@myfunction);

end

function output_txt = myfunction(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global ind method nfft
index_logic = ind;
pos = get(event_obj,'Position');

dataObj=get(gca,'children');
xData=get(dataObj,'Xdata');
yData=get(dataObj,'Ydata');
a=[];c=[];
%a=get(dataObj,'Zdata');
%c=get(dataObj,'Cdata');

index=find(xData ==  pos(1) & yData ==  pos(2));

% [row,~] for index[n*1]
% [~,col] for index[1*n]
[row,~] = find(index_logic,index); 
idx = row(end);

if ~isempty(a) && ~isempty(c)
    output_txt = {['x: ',num2str(pos(1),4)],...
        ['y: ',num2str(pos(2),4)],...
        ['a: ',num2str(a(index),4)],...
        ['c: ',num2str(c(index),6)],...
        ['index: ',num2str(idx,6)]};
    
elseif  isempty(a) && ~isempty(c) && length(xData)==length(c)
    output_txt = {['x: ',num2str(pos(1),4)],...
        ['y: ',num2str(pos(2),4)],...
        ['c: ',num2str(c(index),6)],...
        ['index: ',num2str(index,6)]};
else
    output_txt = {['x: ',num2str(pos(1),4)],...
        ['y: ',num2str(pos(2),4)],...
        ['index: ',num2str(idx,6)]};
end

%[choc,channel,plateau] = ind2tab(idx);
[f,S] = spectrum_norm(idx);


    fitting(idx,method);
    

% Fitting with generalized Gaussion + Gaussion
if method==7
% 7 parameters
x0 = [max(S),0,100,1,max(S),2,min(S)];
lb = [0,-300,0,0,0,0,0];
ub = [max(S),300,500,10,max(S),10,1];
[x,resnorm] = lsqnonlin(@myfun_fitting_7p,x0,lb,ub,[],f,S);
FB = x(1)*exp( - 0.5*( abs((f-x(2)))./x(3) ).^x(4));
FC = x(5)*exp( - 0.5*( abs((f))./x(6)).^2);
FN = ones(nfft,1).*x(7);
elseif method==8
%8 parameters
        x0 = [max(S),0,100,1,max(S),0,2,min(S)];
        lb = [0,-300,0,0,0,-20,0,0];
        ub = [max(S),300,500,10,max(S),20,10,1];
        [x,resnorm] = lsqnonlin(@myfun_fitting_8p,x0,lb,ub,[],f,S);
FB = x(1)*exp( - 0.5*( abs((f-x(2)))./x(3) ).^x(4));
FC = x(5)*exp( - 0.5*( abs((f-x(6)))./x(7)).^2);
FN = ones(nfft,1).*x(8);
end

% sum of components
F = FB + FC + FN;

% Power of ocmponents
PowerFB = trapz(f,FB);
PowerFC = trapz(f,FC);
PowerFN = trapz(f,FN);

% plots
figure;
set(gcf,'Pos',[602   143   741   537]);
plot(f,10*log10(S),'LineWidth',1.5);hold on;
plot(f,10*log10(F),'k','LineWidth',1.5); hold on;
plot(f,10*log10(FB),'r','LineWidth',1.5); hold on;
plot(f,10*log10(FC),'g','LineWidth',1.5); hold on;
plot(f,10*log10(FN),'y'); 
xlabel('f (KHz)');ylabel('PSD');
title(['@',num2str(idx)]);
axis([-500 500 min(10*log10(S))*1.1 max(10*log10(S))]);
legend('Spectrum','Fitting','Broadband','Central Peak','Noise');
% load('database_r.mat');
% title(['@',num2str(idx),', r = ',num2str(database_r.r(idx),3)]);  

[~,~] = spectrum(idx);
% Print 

disp(['Fitting parameters: ',num2str(x,3)]);
disp(['Spectrum Index: ',num2str(idx)]);
disp(['Resnorm: ',num2str(resnorm)]);
disp(['PowerFB: ',num2str(PowerFB)]);
disp(['PowerFC: ',num2str(PowerFC)]);
disp(['PowerFN: ',num2str(PowerFN)]);

end