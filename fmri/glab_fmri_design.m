function SPM = glab_fmri_design(EXPT,subj,model)
    
    % extract subject-specific info from EXPT
    S = EXPT.subject(subj);
    
    % make analysis directory if one doesn't already exist
    adir = fullfile(EXPT.analysis_dir,S.name);
    if ~exist(adir,'dir'); mkdir(adir); end
    
    SPM.xBF.UNITS = 'secs';         % set units to seconds
    
    nSess = length(S.functional);   % number of functional sessions
    
    for s = 1:nSess
        
        % get user-defined regressors
        SPM.Sess(s).U = EXPT.make_regressors(subj,model);
        
        % get motion regressors
        rp = dlmread(fullfile(S.functional(1).niftidir,sprintf('rp_
        SPM.Sess(s).C.C
    end