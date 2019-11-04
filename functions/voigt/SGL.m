function output=SGL(lorentzShare,width,x0,x)
    %%defining Gauss + Lorentz Function
    % designations are according to http://www.casaxps.com/help_manual/line_shapes.htm
    % 01.02.13 CherkasovN
    % INPUT
    % Lorantzian fraction
    % Width
    % Position
    % X
    
    E=x0;
    F=width;
    m=lorentzShare/100;
    output=(1-m)*exp(-4*log(2)*(x-E).^2./F^2) + m ./(1+4*(x-E).^2./F^2);    
end
