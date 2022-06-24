%% Helper function for implementation of DBM

%  generate a (discretized) distribution for the bernoulli rate parater $r$,
%  as a broad beta distribution with mean $r0$, and scale parameter bscale. 
%  used as a prior distribution.

%  Returns: discretized distribution and its variance.



function [betaprior, betavar] = getprior(rzero, N, bscale)

if nargin < 3; bscale = 2; end

%% make sure not degenerate mean.
if rzero == 0 || rzero == 1 
    error('getprior:: extreme value for prior mean');
end

%% Another check:: make sure palph pbet >= 1;
if min(rzero, (1-rzero)) < 1/bscale
     bscale = 1./min(rzero, (1-rzero));
     warning('getprior_r: changing scale parameter');
end
%% Compute prior dist

param_alph = rzero*bscale;
param_bet = (1-rzero)*bscale;

betaprior = pdf('beta', [0:N]/N, param_alph, param_bet);
betaprior = betaprior/sum(betaprior);

x = [0:N]./N;
betavar =  sum(betaprior.*((x-rzero).^2));


end
