function Rtt = taylorfit( tau,tL,K )
% Taylor function as described by Hennequin 1999
%   Detailed explanation goes here
tau=tau;
t=tau/tL;
K2=K*K;
Rtt=exp(-K2*(t-1+exp(-t)));
end

%% backup
% function Rtt = taylorfit( tau,tL,K )
% % Taylor function as described by Hennequin 1999
% %   Detailed explanation goes here
% t=tau/tL;
% K2=K*K;
% Rtt=exp(-K2*(t-1+exp(-t)));
% end
