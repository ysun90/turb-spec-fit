
%% Open parallel computing
parpool(16)

%% Parameters
N = 358874; 
FVAL = zeros(1,N);
BIC = zeros(1,N);

AmpCS = zeros(1,N);
AmpLF = zeros(1,N);
AmpBB = zeros(1,N);
muCS = zeros(1,N);
muLF = zeros(1,N);
muBB = zeros(1,N);
sigmaCS = zeros(1,N);
sigmaLF = zeros(1,N);
k2DBB = zeros(1,N);
tauBB = zeros(1,N);
Noise = zeros(1,N);

ES0 = zeros(1,N);
ECS = zeros(1,N);
ELF = zeros(1,N);
EBB = zeros(1,N);
EN = zeros(1,N);
SDCS = zeros(1,N);
SDLF = zeros(1,N);
SDBB = zeros(1,N);

%% Fitting
tic;       
parfor ii=170001:190000
    
    [E, X, X0, T, S] = fitmincon(ii,'bisTD','multi',0.5,1)
    
    FVAL(ii) = E.FVAL; 
    BIC(ii) = E.BIC;
    
    AmpCS(ii) = X.ampBB;
    AmpLF(ii) = X.ampLF;
    AmpBB(ii) = X.ampBB;
    muCS(ii) = X.muCS;
    muLF(ii) = X.muLF;
    muBB(ii) = X.muBB;
    sigmaCS(ii) = X.sigmaCS;
    sigmaLF(ii) = X.sigmaLF;
    k2DBB(ii) = X.k2DBB;
    tauBB(ii) = X.tauBB;
    Noise(ii) = X.noise;
    
    ES0(ii) = T.ES0;
    ECS(ii) = T.ECS;
    ELF(ii) = T.ELF;
    EBB(ii) = T.EBB;
    EN(ii) = T.EN;
    SDCS(ii) = T.SDCS;
    SDLF(ii) = T.SDLF;
    SDBB(ii) = T.SDBB;
    
end
time = toc; 
save('T170001to358874.mat')