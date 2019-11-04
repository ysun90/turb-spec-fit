function [para_temp,tps_temp] = time_unique(S,t)

if length(unique(t))~=length(t)
    [~,ind] = unique(t);
    tps_temp = t(ind);
    para_temp = S(ind);
else
    tps_temp = t;
    para_temp = S;
end

end
