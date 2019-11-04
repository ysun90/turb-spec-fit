function Fobj = objective_function(x,F,S0,S,model,method,wlinear,wcenter)

switch model
    case 'GGD'
        Sfit = x(1)*gaussian_function(F,x(2),x(3)) + ... % Central Spike (CS)
            x(4)*gaussian_function(F,x(5),x(6)) + ... % Low Frequency (LF)
            x(7)*ggaussian_function(F,x(8),x(9),x(10)) + ... % Broadband (BB)
            x(11);
    case {'VD','PVD2VD'}
        Sfit = x(1)*gaussian_function(F,x(2),x(3)) + ... % Central Spike (CS)
            x(4)*gaussian_function(F,x(5),x(6)) + ... % Low Frequency (LF)
            x(7)*voigt(F,x(8),x(9),x(10)) + ...       % Broadband (BB)
            x(11);                                    % Noise (N)
        
    case 'PVD'
        Sfit = x(1)*gaussian_function(F,x(2),x(3)) + ... % Central Spike (CS)
            x(4)*gaussian_function(F,x(5),x(6)) + ... % Low Frequency (LF)
            x(7)*pvoigt(F,x(8),x(9),x(10)) + ...       % Broadband (BB)
            x(11);                                    % Noise (N)
        
    case 'TD'
        Sfit = x(1)*gaussian_function(F,x(2),x(3)) + ... % Central Spike (CS)
            x(4)*gaussian_function(F,x(5),x(6)) + ... % Low Frequency (LF)
            x(7)*taylorFFT(F,x(8),x(9),x(10)) + ...       % Broadband (BB)
            x(11);                                    % Noise (N)
    case 'bisTD'
        t=(-512:512)';
        Sfit = x(1)*gaussian_function(F,x(2),x(3)) + ... % Central Spike (CS)
            x(4)*gaussian_function(F,x(5),x(6)) + ... % Low Frequency (LF)
            x(7)*taylorbisFFT(t,x(8),x(9),x(10)) + ...       % Broadband (BB)
            x(11);
end

% Cost function
    Fobj = cost_function(F,S0,S,Sfit,method,wlinear,wcenter);
end