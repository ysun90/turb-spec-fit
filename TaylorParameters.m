%% BisT
clear SigmaMatrix KtsMatrix
nu=exp(log(0.002):0.2:log(5)); 
tau=exp(log(0.1):0.1:log(2));
%Sigma = zeros(length(tau),length(kappa));
[Nu,Tau]=meshgrid(nu,tau);
t=0:1024;
F=-0.5:1/1024:0.5;
for ii=1:length(nu)
    for jj=1:length(tau)
%         SigmaMatrix(jj,ii)=1e3*std_pdf_norm(F,taylorFFT(t,0,tau(jj),kappa(ii)));
%         KtsMatrix(jj,ii)=kurtosis_pdf_norm(F,taylorFFT(t,0,tau(jj),kappa(ii)));
%        [dummy,SigmaMatrix(jj,ii),KtsMatrix(jj,ii)]=moment124_pdf(F*1e3,taylorFFT(t,0,kappa(ii),tau(jj)));
       TT=taylorbisFFT(t,0,nu(ii),tau(jj));
       [dummy,SigmaMatrix(jj,ii),KtsMatrix(jj,ii)]=moment124_pdf(F*1e3,TT);
  %     plot(F,lg(TT));hold on;pause(0.05);
    end
 %   title(sprintf('kappa=%f',kappa(ii)));pause;hold off;
end
figure(1);surf(Nu,Tau,SigmaMatrix);
xlabel('\nu=\kappa*\tau');
ylabel('\tau');
zlabel('\sigma_{spec}');
figure(2);surf(Nu,Tau,KtsMatrix);
xlabel('\nu=\kappa*\tau');
ylabel('\tau');
%xlim([0 0.05]);ylim([0 0.1])
zlabel('k_{spec}');

%% T
clear SigmaMatrix KtsMatrix
kappa=exp(log(0.002):0.2:log(5)); 
tau=exp(log(0.1):0.1:log(2));
%Sigma = zeros(length(tau),length(kappa));
[Kappa,Tau]=meshgrid(kappa,tau);
t=-512:512;
F=-0.5:1/1024:0.5;
for ii=1:length(kappa)
    for jj=1:length(tau)
%         SigmaMatrix(jj,ii)=1e3*std_pdf_norm(F,taylorFFT(t,0,tau(jj),kappa(ii)));
%         KtsMatrix(jj,ii)=kurtosis_pdf_norm(F,taylorFFT(t,0,tau(jj),kappa(ii)));
%        [dummy,SigmaMatrix(jj,ii),KtsMatrix(jj,ii)]=moment124_pdf(F*1e3,taylorFFT(t,0,kappa(ii),tau(jj)));
       TT=taylorFFT(t,0,kappa(ii),tau(jj));
       [dummy,SigmaMatrix(jj,ii),KtsMatrix(jj,ii)]=moment124_pdf(F*1e3,TT);
  %     plot(F,lg(TT));hold on;pause(0.05);
    end
 %   title(sprintf('kappa=%f',kappa(ii)));pause;hold off;
end
figure(1);surf(Kappa,Tau,SigmaMatrix);
xlabel('\kappa');
ylabel('\tau');
zlabel('\sigma_{spec}');
figure(2);surf(Kappa,Tau,KtsMatrix);
xlabel('\kappa');
ylabel('\tau');
%xlim([0 0.05]);ylim([0 0.1])
zlabel('k_{spec}');