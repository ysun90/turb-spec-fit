function X0Matrix = Multi_initial_guess(X0,MODEL,F,S)

        X0Matrix = repmat(X0,5,1);
        
switch MODEL        
    case 'GGD'
        X0Matrix(2,10) = X0(10)*4;
        X0Matrix(3,10) = X0(10)/4;
        X0Matrix(4,10) = X0(10)*2;
        X0Matrix(5,10) = X0(10)/2;
    case {'VD','PVD','PVD2VD'}
        X0Matrix(2,9) = 10; X0Matrix(2,10) = 100;
        X0Matrix(3,9) = 50; X0Matrix(3,10) = 50;
        X0Matrix(4,9) = 75; X0Matrix(4,10) = 25;
        X0Matrix(5,9) = 25; X0Matrix(5,10) = 75;

    case 'TD'
        X0Matrix(2,9) = X0(9)*16;
        X0Matrix(3,9) = X0(9)/1.8;
        X0Matrix(4,9) = X0(9)*8;
        X0Matrix(5,9) = X0(9)/1.6;
    case 'bisTD'
%         tau=1;
%         
%             k2Dtab=exp(log(0.01):0.2:log(100));
%             t=-512:512;
%             Fnorm=-0.5:1/1024:0.5;
%             for ii=1:length(k2Dtab)
%                 [m1,m2,m4]=moment124_pdf(Fnorm*1e3,taylorbisFFT(t,0,k2Dtab(ii),tau));
%                 SIGMAtab(ii)=m2;
%             end
%             indBB = (abs(F)>20 & abs(F)<=300);
%             sigmaSPEC=std_pdf(F(indBB),S(indBB))/sum(S(indBB));
%             [~,iU1]=unique(SIGMAtab);[~,iU2]=unique(k2Dtab);
%             iU=intersect(iU1,iU2);
%             k2D_estimated = interp1(SIGMAtab(iU),k2Dtab(iU),sigmaSPEC,'linear','extrap');
%             k2D=k2D_estimated;
        
            
           X0Matrix(2,10)=X0(10)/10;
           X0Matrix(3,10)=X0(10)*10;
           X0Matrix(4,10)=X0(10)/5;
           X0Matrix(5,10)=X0(10)*5;
end

end

%         tau=0.1:0.1:1;
%         for jj=1:length(tau)
%             k2Dtab=exp(log(0.01):0.2:log(100));
%             t=-512:512;
%             Fnorm=-0.5:1/1024:0.5;
%             for ii=1:length(k2Dtab)
%                 [m1,m2,m4]=moment124_pdf(Fnorm*1e3,taylorbisFFT(t,0,k2Dtab(ii),tau(jj)));
%                 SIGMAtab(ii)=m2;
%             end
%             indBB = (abs(F)>20 & abs(F)<=300);
%             sigmaSPEC=std_pdf(F(indBB),S(indBB))/sum(S(indBB));
%             [~,iU1]=unique(SIGMAtab);[~,iU2]=unique(k2Dtab);
%             iU=intersect(iU1,iU2);
%             k2D_estimated = interp1(SIGMAtab(iU),k2Dtab(iU),sigmaSPEC,'linear','extrap');
%             k2D(jj)=k2D_estimated;
%         end