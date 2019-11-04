function VF = voigt_function(F,mu,sigma,gamma)

        fG=sigma*sqrt(2*log(2));
        fL=gamma;
        
        VF = voigt(F,mu,fG,fL);
        VF=VF/sum(VF);

end