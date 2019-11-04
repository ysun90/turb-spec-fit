function [fpV] = pvoigt(F,mu,sigma,gamma)
        
        fG=2*sigma*sqrt(2*log(2));
        fL=2*gamma;
        
        f = (fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + ...
                4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5).^(1/5);
        eta = 1.36603*(fL/f) - 0.47719*(fL/f)^2 + 0.11116*(fL/f)^3;
        
        fpV = (2*eta/pi/f)*((1+4*((F-mu)./f).^2).^(-1)) + ...
                          (1-eta)*(2/f)*((log(2)/pi).^(0.5))*exp(-4*log(2)*(((F-mu)./f).^2));
        fpV=fpV/sum(fpV);

end