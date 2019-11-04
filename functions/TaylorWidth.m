function HWHM = TaylorWidth(f,A,mu,tauL,K)
    
    F = A.*taylorFFT(f,mu,tauL,K);
    a = find(F>=(max(F)/2));
    HWHM = (f(a(end))-f(a(1)))/2;
    
end
