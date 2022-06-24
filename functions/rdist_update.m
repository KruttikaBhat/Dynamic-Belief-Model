%%% code for updating the distribution over the bernoulli parameter $r$. 

% Input: distribution before current trial, current binary obs ($obs)
% a "reset probability" 1-alpha, and a "prior" distribution to reset to.
% (both distributions are discretized representations.)

% Output: Bayes estimate of parameter $r$ (i.e., probability-weighted average)
%       Distribution over parameter estimate
%       Variance of parameter estimate

function [bayes_r, rdist_post, varest] = rdist_update(rdist_pre, obs_t, betaprior, alph)

global params; % model-specific code; ignore.
if nargin < 4; alph = params.rdist.alpha; end

rdist_post = compute_rpost(rdist_pre, obs_t, betaprior,alph);
[bayes_r, varest] = getBayesEst(rdist_post);

end

function map_r = getMAPest(rdist)
    map_r = find(max(rdist) == rdist);
    map_r = map_r(1)/length(rdist);

end
function [bayes_r, varest] = getBayesEst(rdist)
    N = length(rdist)-1;
    prob = [0:N]/N;
    bayes_r = sum(prob.*rdist);

    varest = sum(rdist.* ((prob-bayes_r).^2));

end


%% sequential nonstationary update (DBM, yu/cohen08).
% straightforward numerical calculations:
% 1. discretized representation of estimated pdf for bernoulli parameter $r$
% 2. 0/1 observation on current trial
% 3. probability of "reset" $1-alph$, to a distribution $betaprior$
% normalize and return discretized/approximate posterior distribution.

function rdist_post = compute_rpost(rdist_pre, obs_t, betaprior, alph)

rdist_pre = rdist_pre(:)';

Nr = length(rdist_pre)-1;

newp_r = alph*rdist_pre + (1-alph)*betaprior;
newp_r = newp_r/sum(newp_r); % normalize 


x=[0:Nr]/Nr;
obslik = obs_t*x + (1-obs_t)*(1-x);

%% update with current data point.
% p(rt|x1..t) = p(xt|rt).*p(rt|x1..t-1); 
% p(rt|x1..t) = p(rt|x1..t)/sum(p(rt|x1..t));

rdist_post = newp_r.*obslik;  
rdist_post = rdist_post./sum(rdist_post); % normalize.

end