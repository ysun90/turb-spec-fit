function output=pearson(constant,x0,x)
    % defining Pearson VII Function
    % designations are according to http://pd.chem.ucl.ac.uk/pdnn/peaks/pvii.htm
    % 
    % output=pearson(constant,x0,x) 
    %   constants=[Ymax w m]
    %
    % in Iglesia form: Y=a./(p.*(1+(x-x0).^2.*(2.^(1/q)-1)./p).^q);
    % recalculation to Iglesia: q=m, p=w^2, a=Ymax*p
    % 27.01.13 CherkasovN
    
    Ymax=constant(1);
    w=constant(2);
    m=constant(3);
    Y=Ymax./(1+(2^(1/m)-1).*(x-x0).^2./w.^2).^m;
    output=Y;
end
