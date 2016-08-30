function results = dpkf(Y,opts)
    
    % Dirichlet process Kalman filter algorithm. Uses a local maximum a
    % posteriori estimate of the partition.
    %
    % USAGE: results = dpkf(Y,[opts])
    %
    % INPUTS:
    %   Y - [T x D] observation sequence, where Y(t,d) is dimension d of the observation at time t
    %   opts (optional) - structure with any of the following fields
    %                     (missing fields are set to defaults):
    %                       .R = noise covariance (default: eye(D))
    %                       .Q = diffusion covariance (default: 0.01*eye(D))
    %                       .W = dynamics matrix (default: eye(D))
    %                       .C = prior state covariance (default: 10*eye(D))
    %                       .alpha = concentration parameter (default: 0.1)
    %                       .sticky = stickiness of last mode (default: 0)
    %                       .x0 = prior mean (default: zeros(1,D))
    %                       .Kmax = upper bound on number of state (default: 10)
    %                   Note: if R, Q, W, or C are given as scalars, it is
    %                         assumed that they are the same for all dimensions
    %
    % OUTPUTS:
    %   results - [1 x T] structure with the following fields:
    %             .P = [1 x Kmax] cell array of [D x D] posterior
    %                  covariance matrices (one for each mode)
    %             .x - [K x D] matrix of posterior mean state estimates
    %             .G - [1 x Kmax] cell array of [D x D] Kalman gain
    %             .err - [1 x D] prediction error vector
    %             .pZ - [1 x K] posterior probability of modes
    %
    % Sam Gershman, June 2015
    % Reference: Gershman, Radulescu, Norman & Niv (2014). Statistical
    % computations underlying the dynamics of memory updating. PLOS
    % Computational Biology.
    
    [T,D] = size(Y);
    
    % set missing parameters to their defaults
    def_opts.Q = 0.01;
    def_opts.W = 1;
    def_opts.R = 1;
    def_opts.C = 10;
    def_opts.alpha = 0.1;
    def_opts.sticky = 0;
    def_opts.x0 = zeros(1,D);
    def_opts.Kmax = 10;
    F = fieldnames(def_opts);
    if nargin < 2 || isempty(opts)
        opts = def_opts;
    else
        for i = 1:length(F); if ~isfield(opts,F{i}); opts.(F{i}) = def_opts.(F{i}); end; end
    end
    
    % if scalar parameters are given, assume these are the same for all dimensions
    if isscalar(opts.x0); opts.x0 = zeros(1,D)+opts.x0; end
    if isscalar(opts.Q); opts.Q = diag(zeros(1,D)+opts.Q); end
    if isscalar(opts.W); opts.W = diag(zeros(1,D)+opts.W); end
    if isscalar(opts.R); opts.R = diag(zeros(1,D)+opts.R); end
    if isscalar(opts.C); opts.C = diag(zeros(1,D)+opts.C); end
    
    % initialization
    for k = 1:opts.Kmax
        P{k} = opts.C;
        x(k,:) = opts.x0;
    end
    pZ = [1 zeros(1,opts.Kmax-1)];  % mode 1 starts with probability 1
    M = [1 zeros(1,opts.Kmax-1)];
    lik = zeros(1,opts.Kmax);
    
    % run DPKF
    for t = 1:T
        
        x = x*opts.W;           % predicted (a priori) estimate
        err = Y(t,:) - pZ*x;    % prediction error
        for k = 1:opts.Kmax
            P{k} = opts.W*P{k}*opts.W' + opts.Q;    % predicted (a priori) estimate covariance
        end
        
        if all(~isnan(err))
            
            % compute posterior over modes
            if t > 1 && opts.alpha > 0
                
                % Chinese restaurant prior
                prior = M;
                prior(find(prior==0,1)) = opts.alpha;   % probability of new mode
                prior(k) = prior(k) + opts.sticky;      % make last mode sticky
                prior = prior./sum(prior);
                
                % multivariate Gaussian likelihood
                for k = 1:opts.Kmax
                    lik(k) = mvnpdf(Y(t,:),x(k,:),P{k}+opts.R);
                end
                
                % posterior
                pZ = prior.*lik;
                pZ = pZ./sum(pZ);
                
                % MAP estimate
                [~,k] = max(pZ);
                M(k) = M(k) + 1;
            end
            
            % update estimates
            for k = 1:opts.Kmax
                S{k} = (pZ(k)^2)*P{k} + opts.R;         % error covariance
                G{k} = (P{k}*pZ(k))/S{k};               % Kalman gain
                x(k,:) = x(k,:) + err*G{k};             % updated (a posteriori) estimate
                P{k} = P{k} - pZ(k)*G{k}*P{k};          % updated (a posteriori) estimate covariance
            end
        end
        
        % store results
        results(t).P = P;
        results(t).x = x;
        results(t).G = G;
        results(t).pZ = pZ;
        results(t).err = err;
        
    end