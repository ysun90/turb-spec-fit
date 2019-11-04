function [R,rc] = PlotCorrelation(d,P1,P2,index)

rMin = -1;
rMax = 1;
rInterval = 0.1;
rStep = 0.02;

rL = rMin:rStep:(rMax-rInterval);
rU = rL + rInterval;

R = zeros(1,length(rL));

for gi = 1:length(rL)
    % clean database
    rl = rL(gi);
    ru = rU(gi);
    indr = d.rhoc>rl & d.rhoc<ru;  

    % calculate the correlation coefficients 
    Rcc = corrcoef(P1(indr & index),P2(indr & index),'rows','complete');
    R(gi) = Rcc(1,2);

end

figure();
rc = (rL+rU)./2;
plot(rc,R);

end
