function files = get_files(filename)
    
    % Get filenames matching a regular expression.
    %
    % USAGE: files = get_files(filename)
    %
    % INPUTS:
    %   filename - regular expression specifying set of file names
    %
    % OUPUTS:
    %   files - character array of file names (including full path)
    %
    % Sam Gershman, June 2015
    
    files = dir(filename);
    filedir = fileparts(filename);
    files = dir2char(files,filedir);