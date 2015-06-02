function X = reflecting_random_walk(N,D,s,b)
    
    % Generate samples from a Gaussian random walk with reflecting boundaries.
    %
    % USAGE: X = reflecting_random_walk(N,D,s,b)
    %
    % INPUTS:
    %   N - number of samples
    %   D - number of variables
    %   s - standard deviation of the random walk
    %   b - [1 x 2] lower and upper boundaries
    %
    % OUTPUTS:
    %   X - [N x D] samples from the random walk
    %
    % Sam Gershman, June 2015
    
    X = zeros(N,D);
    X(1,:) = unifrnd(b(1),b(2),1,D);    % initialize uniformly within bounds
    
    for n = 2:N
        x = X(n-1,:) + normrnd(0,s,1,D);
        x(x<b(1)) = 2*b(1) - x(x<b(1));
        x(x>b(2)) = 2*b(2) - x(x>b(2));
        X(n,:) = x;
    end