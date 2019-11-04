

load('T1to170000.mat');
TF=load('T170001to358874.mat');

FVAL(170001:end)=TF.FVAL(170001:end);
BIC(170001:end)=TF.BIC(170001:end);

AmpCS(170001:end)=TF.AmpCS(170001:end);
AmpLF(170001:end)=TF.AmpLF(170001:end);
AmpBB(170001:end)=TF.AmpBB(170001:end);
muCS(170001:end)=TF.muCS(170001:end);
muLF(170001:end)=TF.muLF(170001:end);
muBB(170001:end)=TF.muBB(170001:end);
sigmaCS(170001:end)=TF.sigmaCS(170001:end);
sigmaLF(170001:end)=TF.sigmaLF(170001:end);
k2DBB(170001:end)=TF.k2DBB(170001:end);
tauBB(170001:end)=TF.tauBB(170001:end);
Noise(170001:end)=TF.Noise(170001:end);

ES0(170001:end)=TF.ES0(170001:end);
ECS(170001:end)=TF.ECS(170001:end);
ELF(170001:end)=TF.ELF(170001:end);
EBB(170001:end)=TF.EBB(170001:end);
EN(170001:end)=TF.EN(170001:end);
SDCS(170001:end)=TF.SDCS(170001:end);
SDLF(170001:end)=TF.SDLF(170001:end);
SDBB(170001:end)=TF.SDBB(170001:end);

save('T358874.mat','FVAL','BIC','AmpCS','AmpLF','AmpBB',...
'muCS','muLF','muBB','sigmaCS','sigmaLF','k2DBB','tauBB','Noise',...
'ES0','ECS','ELF','EBB','EN','SDCS','SDLF','SDBB','r');