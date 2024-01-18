function ind = sample(w)
%% SAMPLE - Obtains an index as a sample from a discrete probability distribution
%
% Syntax:
%   ind = sample(w)
%
% In:
%   w   - A vector of weights representing a discrete probability distribution
%
% Out:
%   ind - The index of the sampled element
%
% Description:
%   This function samples an index from a discrete probability distribution 
%   defined by the input vector 'w'. Each element of 'w' represents the 
%   probability of selecting the corresponding index. The function uses 
%   cumulative summation and uniform random sampling to select an index based 
%   on the provided probabilities.
%
% Examples:
%   ind = sample([0.1, 0.2, 0.3, 0.4]);
%   ind = sample([0.25, 0.25, 0.25, 0.25]);
%
% See also:
%   rand, cumsum
%
% Copyright:
%   2018-   Manon Kok

wc = cumsum(w);
u = rand;
ind = sum(wc < u) + 1;

end

% %% Test functionality of function
% % Make some weights
% N_P = 100;
% w = rand(N_P,1);
% w = w ./sum(w);
% 
% % Compute cumsum
% wc = cumsum(w);
% 
% % Draw samples
% N_S = 100000;
% ind = zeros(N_S,1);
% for i = 1:N_S
%     u = rand;
%     ind(i) = sum(wc < u) + 1;
% end
% 
% % Compute ratio of how often each sample has been drawn
% nw = zeros(N_P,1);
% for i = 1:N_P
%     nw(i) = sum(ind == i);
%     nw(i) = nw(i) / N_S;
% end
% 
% % Compare ratio with weights
% close all
% plot(w)
% hold all
% plot(nw)