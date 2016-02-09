function matlabbatch = glab_fmri_preproc(EXPT,subj)
    
    % Preprocess fMRI data. This function does the following: (1)
    % realignment, (2) segmentation, (3) coregistration of functionals to
    % structural, (4) normalization, (5) smoothing.
    %
    % USAGE: glab_fmri_preproc(EXPT,subj)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   subj - subject number
    %
    % OUTPUTS:
    %   matlabbatch - batch structure for use in spm_jobman
    %   Realignment and coregistration induce changes in the headers of the
    %   nifti files. Normalization writes out new files with 'w' prefix.
    %   Smoothing writes out new files with 's' prefix.
    %
    % Sam Gershman, June 2015
    % Adapted from preproc_fmri_simplified.m in SPM12
    
    %% Setup
    S = EXPT.subject(subj);
    matlabbatch = {};
    nSess = length(S.functional);   % number of functional sessions
    
    %% Realignment
    iRU = 1;
    for s = 1:nSess
        niftidir = S.functional(s).niftidir;
        sess = S.functional(s).sess;
        scans = get_files(fullfile(niftidir,sprintf('*.%d.*',sess)));
        matlabbatch{iRU}.spm.spatial.realignunwarp.data(s).scans = scans;
        sessdep(s) = cfg_dep(sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', s), ...
            substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sess', '()',{s}, '.','uwrfiles'));
    end
    
    %% Segmentation
    iSeg = length(matlabbatch) + 1;
    matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
    ngaus  = [1 1 2 3 4 2];
    for c = 1:6 % tissue class c
        matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
            fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
        matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
        matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [0 0];
        matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
    end
    matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];
    
    %% Coregister mean functional to bias-corrected structural
    iCR = length(matlabbatch) + 1;
    matlabbatch{iCR}.spm.spatial.coreg.estimate.ref(1) = ...
        cfg_dep('Segment: Bias Corrected (1)', ...
        substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{iCR}.spm.spatial.coreg.estimate.source(1) = ...
        cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
        substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','meanuwr'));
    matlabbatch{iCR}.spm.spatial.coreg.estimate.other = sessdep;
    
    %% Normalise functional images
    iNF = length(matlabbatch) + 1;
    matlabbatch{iNF}.spm.spatial.normalise.write.subj.def(1) = ...
        cfg_dep('Segment: Forward Deformations', ...
        substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','fordef', '()',{':'}));
    matlabbatch{iNF}.spm.spatial.normalise.write.subj.resample = sessdep;
    
    %% Smooth normalised functionals
    iSm = length(matlabbatch) + 1;
    matlabbatch{iSm}.spm.spatial.smooth.data(1) = ...
        cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
        substruct('.','val', '{}',{iNF}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('()',{1}, '.','files'));
    
    %% Run jobs
    spm_jobman('run', matlabbatch)