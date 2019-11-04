function GD = gaussian_function(F,mu,sigma)
    % amp - amplitude
    % mu - mean (central position)
    % sigma - stanfard deviation
    % F - frequency
    GD = 1/sqrt(2*pi)/sigma*exp( -0.5*( ( ( F-mu )/sigma ).^2 ) );
end