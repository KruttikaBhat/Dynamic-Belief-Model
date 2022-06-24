% Estimate the probability of an event (the bernoulli parameter $r$ (i.e. gamma in Yu,Cohen 2008)) 
% as the mean over a distribution updated with observations.

% return value is the PRIOR expectation at the onset of each trial.

function [rtraj, postdist, vartraj] = est_prob(stopseq, alph, priormean, scale)

N = 1000;

[rprior, varzero] = getprior(priormean,N, scale); % computes prior distribution
rpost = rprior;                                   % Initially sets posterior = prior   

rtraj(1) = priormean; vartraj(1) = varzero;  

postdist = zeros(length(rprior), length(stopseq));

for i = 2:length(stopseq)
   [rtraj(i), rpost, vartraj(i)] = rdist_update(rpost, stopseq(i-1), rprior, alph); % updates the posterior distribution on each trial based on the likelihood of the observed event
   postdist(:,i) = rpost;
end

end

