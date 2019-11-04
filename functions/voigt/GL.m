function output=GL(lorentzShare,width,x0,x)
    %% defining pseudo-Voigt Function
    % designations are according to http://www.casaxps.com/help_manual/line_shapes.htm
    % INPUT:
    %  Lorentzian share
    %  width
    %  centerline
    %  wavenumber
    % 01.02.13 CherkasovN
    
    E=x0;
    F=width;
    m=lorentzShare/100;
    output=exp(-4*log(2)*(1-m)/F^2*(x-E).^2)./(1+4*m*(x-E).^2/F^2);
end
