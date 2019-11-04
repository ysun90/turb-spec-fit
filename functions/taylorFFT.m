function S = taylorFFT( tau,wc,K,tL )
% Taylor function as described by Hennequin 1999
%   Detailed explanation goes here
t=tau/tL;
K2=K*K;
tL2=tL*tL;
Rtt=exp( - K2*tL2*( abs(t) -1 + exp(-abs(t)) ) ).*exp(-1i*2*pi*wc*tau);
%plot(tau,Rtt);pause
S = fftshift(abs((fft(Rtt))));
S=S/sum(S);

end