function se = wse(X)
    
    % Within-subject error, following method of Cousineau (2005).
    %
    % USAGE: se = wse(X)
    %
    % INPUTS:
    %   X - [N x D] data with N observations and D subjects
    %
    % OUTPUTS:
    %   se - [1 x D] within-subject standard errors
    %
    % Sam Gershman, June 2015
    
    X = bsxfun(@minus,X,mean(X,2));
    N = sum(~isnan(X));
    se = bsxfun(@rdivide,nanstd(X),sqrt(N));