function [fittingError areaBands plots]=deconvolute3Vt(x,y,nIterations)
% deconvolute3Vt Deconvolution of the Mordenite spectrum with 3 Voigt bands
%   OUTPUT
%     fittingError - array 2*1 - ranges of the fitting error    
%     areaBands - array 6*1 of the band rfactiond achivable using the
%       deconvolution in % of total area. [LFmin LFmax HFmin HFmax TFmin TFmax]
%     plots - array N*8 of the points corresponding to the minimum and
%       maximum adsorbances of the components and their sum.
%       [LFmin HFmin TFmin SUMmin LFmax HFmax TFmax SUMmax]


%% initialisation
options = optimset('Display','off');
rng('shuffle', 'twister');
multiplier=0.04/max(y);
y=y.*multiplier;yMax=max(y);
areaExp=trapz(x,y);

% boundaries of allowed parameter variation - b1=lower, b2=upper
z=[19.1	3587.2	13.8	3607.1	11.7	3619.3	19.4	3588.2	14.4	3607.4	11.9	3622.2	70	30	15];

% green [16.9	3587	13.6	3607	10.7	3619.6	18.2	3587.9	14.4	3607.7	11.2	3622.6	75	45	20];

% yellow [19.1	3587.2	13.8	3607.1	11.7	3619.3	19.4	3588.2	14.4	3607.4	11.9	3622.2	70	30	15]

lorentzian=z(13:15);
b1=[z(1:2)	0       z(3:4)      0       z(5:6)  0   ];  
b2=[z(7:8)  yMax    z(9:10)     yMax    z(11:12)  yMax];

derivatives=zeros(nIterations,4);
constants=zeros(nIterations,10);

% the fitting function - 3 Voigt profiles
fSum = @(c,x)  (c(3).*voigtNormGLform(lorentzian(1),c(1),c(2),x)+...
                c(6).*voigtNormGLform(lorentzian(2),c(4),c(5),x)+ ...
                c(9).*voigtNormGLform(lorentzian(3),c(7),c(8),x)) ;




% Spectrum deconvolution
a_band_range=zeros(10,6);
a_band=zeros(10,3);
errorFit=zeros(10,1);totalUncertainty=0;




%% calculation

%deconvolution USE PARRALEL COMPUTATIONAL TOOLBOX if you can !
parfor i=1:nIterations
    Cinit=rand(1,9).*(b2-b1)+b1;    %random initial parameters
    [a error]=lsqcurvefit(fSum,Cinit, x, y,b1,b2,options);
    constants(i,:)=[a error];
end;
constants=(sortrows(constants,10));  % sorted according to residual

% calculating the fitting error and fractions of the bands
% and preparing the plots
yBandMin=ones(length(x),4).*Inf;yBandMax=zeros(length(x),4);

for i=1:nIterations
    yLF=constants(i,3) *voigtNormGLform(lorentzian(1),constants(i,1),constants(i,2),x);
    yHF=constants(i,6) *voigtNormGLform(lorentzian(2),constants(i,4),constants(i,5),x);
    yTF=constants(i,9) *voigtNormGLform(lorentzian(3),constants(i,7),constants(i,8),x);
    
    areaError=100*trapz(x,abs(y-yLF(:)-yHF(:)-yTF(:))) / areaExp; %fitting error
    
    aLF=100*trapz(x,yLF)/areaExp;   % area fraction of the LF component
    aHF=100*trapz(x,yHF)/areaExp;
    aTF=100*trapz(x,yTF)/areaExp;
    
    derivatives(i,1)=areaError;
    derivatives(i,2)=aLF;
    derivatives(i,3)=aHF;
    derivatives(i,4)=aTF;    
    
    % checking if current components are very large or small
    % to include them into the plots
    
    index=yBandMin(:,1)>yLF; 
    yBandMin(index,1)= yLF(index);
    index=yBandMax(:,1)<yLF(:); 
    yBandMax(index,1)= yLF(index);   
    
    index=yBandMin(:,2)>yHF; 
    yBandMin(index,2)= yHF(index);    
    index=yBandMax(:,2)<yHF; 
    yBandMax(index,2)= yHF(index);  
    
    index=yBandMin(:,3)>yTF; 
    yBandMin(index,3)= yTF(index);    
    index=yBandMax(:,3)<yTF; 
    yBandMax(index,3)= yTF(index);  
    
    index=yBandMin(:,4)>(yLF(:)+      yHF(:)+      yTF(:)); 
    yBandMin(index,4)=  yLF(index)+  yHF(index)+  yTF(index); 
    index=yBandMax(:,4)<(yLF(:)+      yHF(:)+      yTF(:)); 
    yBandMax(index,4)=  yLF(index)+  yHF(index)+  yTF(index);     
end;
plots=[yBandMin, yBandMax]./multiplier;
y=y./multiplier;
fittingError=[min(derivatives(:,1)) max(derivatives(:,1))];

numMarginal=nIterations;
%{
% removing the values, which are too far - with the fitting error more that
margin=2; %2 percent higher that the best deconvolution
minErr=min(derivatives(:,1)); 
derivatives=sortrows(derivatives,1);
numMarginal=find(derivatives(1,:)>margin+minErr,1, 'first');    %index of the first value with the error 2% higher
if isempty(numMarginal), numMarginal=nIterations; else numMarginal=numMarginal-1; end;
%}
% preparing the output on the fractions of the components
areaBands(1)=min(derivatives(1:numMarginal,2)); 	%LF min - the minimal ALF value achievable within the accepted error range
areaBands(2)=max(derivatives(1:numMarginal,2)); 	%LF max - the max ALF value achievable within the accepted error range
areaBands(3)=min(derivatives(1:numMarginal,3)); 	%HF min - the minimal AHF value achievable within the accepted error range
areaBands(4)=max(derivatives(1:numMarginal,3)); 	%HF max - the max AHF value achievable within the accepted error range
areaBands(5)=min(derivatives(1:numMarginal,4)); 	%TF min - the minimal ATF value achievable within the accepted error range
areaBands(6)=max(derivatives(1:numMarginal,4)); 	%TF max - the max ATF value achievable within the accepted error range


end
    