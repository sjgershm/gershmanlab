function Z = dpmm_bernoulli(X,varargin)
    
    % Nonparametric Bayesian clustering of binary data. Uses Gibbs sampling to approximate
    % the posterior of a Dirichlet process mixture of Bernoullis. Also
    % samples the concentration parameter (alpha) which controls the number
    % of clusters.
    %
    % Generative model:
    %   alpha|a,b ~ inv-gamma(a,b)
    %   z|alpha ~ CRP(alpha)
    %   theta(k,d) ~ Beta(B,B)
    %   X(n,d) ~ Bernoulli(theta(z(n),d))
    %
    % USAGE: Z = dpmm_bernoulli(X,[varargin])
    %
    % EXAMPLES:
    %   >> X = [1 0; 0 1];
    %   >> Z = dpmm_bernoulli(X);
    %   >> Z = dpmm_bernoulli(X,'nIter',200);
    %   >> Z = dpmm_bernoulli(X,'alpha',1);
    %
    % INPUTS:
    %   X - [N x D] binary data
    %
    % OPTIONAL INPUTS (see examples):
    %   'nIter' - number of Gibbs iterations (default: 1000)
    %   'a' - shape parameter for hyperiorprior on concentration parameter (default: 1)
    %   'b' - scale parameter for hyperiorprior on concentration parameter (default: 1)
    %   'B' - parameter for the beta hyperprior (default: 1)
    %   'alpha' - initial concentration parameter (default: 1)
    %
    % OUTPUTS:
    %   Z - [nIter x N x K] cluster assignments, where Z(i,n,k)=1 if
    %   observation n was assigned to cluster k on iteration i
    %
    % Sam Gershman, May 2013
    % updated Oct 2013
    
    % parameters
    nIter = 1000;   % number of sampling iterations
    a = 1; % shape parameter for the hyperprior on alpha
    b = 1; % scale parameter for the hyperprior on alpha
    B = 1; % base distribution with parameter beta
    
    % initialize hyperparameters
    alpha=1;  % concentration parameter
    eta=0.5;  % auxiliary parameter
    
    % process inputs
    if nargin > 1
        for i = 1:2:length(varargin)-1
            eval([varargin{i},'=',num2str(varargin{i+1}),';']);
        end
    end
    
    % initialize cluster assignments
    [N D] = size(X);
    K = N+1;                % max number of clusters
    M = zeros(K,D);         % cluster/feature co-occurence counts
    M(1,:) = sum(X);        % initialize all data points to the first cluster
    k = ones(N,1);          % cluster assignment vector
    Q = [N zeros(1,K-1)];   % number of observations per cluster
    Z = zeros(nIter,N,K);
    
    % run Gibbs sampler
    for i = 1:nIter
        
        if ~mod(i,100)
            disp(['Iteration ',num2str(i)]);
        end
        
        for n = randperm(N)
            % remove current observation from counts
            M(k(n),:) = M(k(n),:) - X(n,:);
            Q(k(n)) = Q(k(n)) - 1;
            
            active = find(Q>0); % active clusters
            new = find(Q==0,1); % new clusters
            C = [active new];
            
            % sample from conditional distribution
            m = bsxfun(@rdivide,M(C,:)'+B,Q(C)+2*B);
            p = log([Q(active) alpha]) + X(n,:)*log(m) + (1-X(n,:))*log(1-m);
            p = exp(p-logsumexp(p,2));
            j = fastrandsample(p);
            k(n) = C(j);        % new cluster assignment
            Z(i,n,k(n)) = 1;
            
            % update counts
            M(k(n),:) = M(k(n),:) + X(n,:);
            Q(k(n)) = Q(k(n)) + 1;
        end
        
        % sample concentration parameter
        K = sum(Q>0);
        alpha = gamrnd(a+K-1,1/(b-log(eta)));
        eta = betarnd(alpha,N);
        
        % make sure these numbers don't get too small
        alpha = max(alpha,10^-10);
        eta = max(eta,10^-10);
    end
    
    % remove inactive clusters
    m = sum(squeeze(mean(Z)));
    Z(:,:,m==0) = [];