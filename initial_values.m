function  X0 = initial_values(F,S,model)
% X0 = initial_values(F,S,model)  
if strcmpi(model,'PVD2VD')
    [EP,XP,X0P,TP,SP] = fitmincon(2296,'PVD','local',0.5,1);
    X0 = struct2array(XP);
else

    % Central spike
    
    ind_CS = (abs(F)<3);
    ampCS_x0 = S(F==0); X0(1) = ampCS_x0;
    muCS_x0 = mean_pdf(F(ind_CS),S(ind_CS)); X0(2) = muCS_x0;
    sigmaCS_x0 = std_pdf(F(ind_CS),S(ind_CS)); X0(3) = sigmaCS_x0;
    
    % Low frequency
    ind_LF = (abs(F)>3 & abs(F)<20);
    ampLF_x0 = max(S(ind_LF)); X0(4) = ampLF_x0;
    muLF_x0 = mean_pdf(F(ind_LF),S(ind_LF)); X0(5) = muLF_x0;
    sigmaLF_x0 = std_pdf(F(ind_LF),S(ind_LF)); X0(6) = sigmaLF_x0;
    
    % Broadband
    ind_BB = (abs(F)>20 & abs(F)<=300);
    ampBB_x0 = max(S(ind_BB)); X0(7) = ampBB_x0;
    muBB_x0 = mean_pdf(F(ind_BB),S(ind_BB)); X0(8) = muBB_x0;
      
    switch model
        case 'GGD'
            sigmaBB_x0 = std_pdf(F(ind_BB),S(ind_BB)); X0(9) = sigmaBB_x0;
            kts = kurtosis_pdf(F(ind_BB),S(ind_BB));           
            betaBB_x0 =  fsolve(@(x)gamma(5/x)*gamma(1/x)/(gamma(3/x)^2)-kts,2);
            X0(10) = betaBB_x0;
        case {'VD','PVD'}          
            sigmaBBG_x0 = 100; X0(9) = sigmaBBG_x0;
            gammaBBL_x0 = 0.1; X0(10)= gammaBBL_x0;            
        case 'TD'
            tau=0.2;
            kappaTAB=exp(log(0.0001):0.1:log(10));
            t=-512:512;
            Fsample=(-0.5:1/1024:0.5);
            for ii=1:length(kappaTAB)
                sigmaTAB(ii)=std_pdf_norm(Fsample*1e3,taylorFFT(t,0,kappaTAB(ii),tau));
            end
            indBB = (abs(F)>20 & abs(F)<=300);
            sigmaSPEC=std_pdf(F(indBB),S(indBB))/sum(S(indBB));
            [~,iU1]=unique(sigmaTAB);[~,iU2]=unique(kappaTAB);
            iU=intersect(iU1,iU2);
            kappa_estimated = interp1(sigmaTAB(iU),kappaTAB(iU),sigmaSPEC,'linear','extrap');
            kappaBB_x0 = kappa_estimated; X0(9) = kappaBB_x0;
            tauBB_x0 = tau; X0(10) = tauBB_x0;
        case 'bisTD'
            tau=0.1;   
            k2Dtab=exp(log(0.01):0.2:log(100));
            t=-512:512;
            Fnorm=-0.5:1/1024:0.5;
            for ii=1:length(k2Dtab)
                [m1,m2,m4]=moment124_pdf(Fnorm*1e3,taylorbisFFT(t,0,k2Dtab(ii),tau));
                SIGMAtab(ii)=m2;
            end
            indBB = (abs(F)>20 & abs(F)<=300);
            sigmaSPEC=std_pdf(F(indBB),S(indBB))/sum(S(indBB));
            [~,iU1]=unique(SIGMAtab);[~,iU2]=unique(k2Dtab);
            iU=intersect(iU1,iU2);
            k2D_estimated = interp1(SIGMAtab(iU),k2Dtab(iU),sigmaSPEC,'linear','extrap');
            X0(9)=k2D_estimated;
            X0(10)= tau;
            
    end
    
    % Noise level
    ind_noise = (abs(F)>400 & abs(F)<450);
    Noise_x0 = mean(S(ind_noise)); X0(11) = Noise_x0;
end
    
end