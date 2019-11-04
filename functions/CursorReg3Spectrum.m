function  CursorReg3Spectrum(p1,p2,p3,indexing,model)

global ind
ind = indexing;
h = figure;
set(h,'Pos',[16   211   560   420]);

x = linspace(0,8,100);
y = linspace(0,4,100);
[xNl,yTe0] = meshgrid(x,y);
zPower = 10.^(model.Coefficients.Estimate(1)).*...
    power(xNl,model.Coefficients.Estimate(2)).*...
    power(yTe0,model.Coefficients.Estimate(3));
surf(xNl,yTe0,zPower);hold on;
shading interp
colormap hsv;
alpha(.3)
zlim([0 1]);

scatter3(p1(indexing),p2(indexing),p3(indexing),10,'filled');

dcm = datacursormode(gcf);
datacursormode toggle
set(dcm,'updatefcn',@myfunction);

end

function output_txt = myfunction(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global ind nfft
index_logic = ind;

pos = get(event_obj,'Position');

dataObj=get(gca,'children');
xData=get(dataObj,'Xdata');
yData=get(dataObj,'Ydata');
zData=get(dataObj,'Zdata');
a=[];c=[];
%a=get(dataObj,'Zdata');
%c=get(dataObj,'Cdata');

index=find(xData{1,1} ==  pos(1) & yData{1,1} ==  pos(2) & zData{1,1} == pos(3));

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

% Fitting with generalized Gaussion + Gaussion
% 7 parameters
x0 = [max(S),0,100,1,max(S),2,min(S)];
lb = [0,-300,0,0,0,0,0];
ub = [max(S),300,500,10,max(S),10,1];
% %8 parameters
% x0 = [max(S),0,100,1,max(S),2,2,min(S)];
% lb = [0,-300,0,0,0,0,0,0];
% ub = [max(S),300,500,10,max(S),10,10,1];
[x,resnorm] = lsqnonlin(@myfun_fitting_7p,x0,lb,ub,[],f,S);
FB = x(1)*exp( - 0.5*( abs((f-x(2)))./x(3) ).^x(4));
FC = x(5)*exp( - 0.5*( abs((f))./x(6)).^2);
FN = ones(nfft,1).*x(7);
F = FB + FC + FN;

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

%[~,~] = spectrum(idx);

% Print 

disp(['Fitting parameters: ',num2str(x,3)]);
disp(['Resnorm: ',num2str(resnorm)]);
disp(['Spectrum Index: ',num2str(idx)]);

end