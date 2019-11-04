function [x,t] = signal(varargin)

% Check input parameters
if nargin == 1
    [numc,channel,spec] = ind2tab(varargin{:});
elseif nargin == 3
    [numc,channel,spec] = varargin{:};
else
    error('Number of input arguments can only be 1 or 3!!');
end

% Load spctrum data for one shot
S = load(['#',num2str(numc),'drefluc.mat']);

% Get time and signal
t = S.txc(:,spec);
if channel==1
    x = S.xc(:,spec);
elseif channel==2
    x = S.xc_2(:,spec);
end

end