function y = fastrandsample(p,n)
    
    % Fast multinomially-distributed random numbers.
    %
    % USAGE: y = fastrandsample(p,[n])
    %
    % INPUTS:
    %   p - [1 x K] probability distribution
    %   n (optional) - number of samples (default: 1)
    %
    % OUTPUTS:
    %   y - [1 x N] random numbers
    %
    % Sam Gershman, June 2015
    
    if nargin < 2; n=1; end
    [~, y] = histc(rand(1,n),[0 cumsum(p)]);