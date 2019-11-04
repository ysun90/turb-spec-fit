function S = taylorbisFFT( tau,wc,k2D,tL )
% Taylor function as described by Hennequin 1999
%   Detailed explanation goes here
t=tau/tL;
% K2=K*K;
% tL2=tL*tL;
% nu=K*tL;
Rtt=exp( - (k2D*tL)*( abs(t) -1 + exp(-abs(t)) ) ).*exp(1i*2*pi*wc*tau);
%figure(10);plot(tau,Rtt);pause;
S = fftshift(abs((fft(Rtt))));
S=S/sum(S);

end