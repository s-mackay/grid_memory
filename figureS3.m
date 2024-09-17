% figureS3(data_dir, read_from_excel)
% INPUTS:
%   data_dir (str):         Path to source data (default: 'source_data').
%
%   read_from_excel (bool): Flag to read from Excel (true) 
%                           or .mat files (false) (default: false).
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function figureS3(data_directory, read_from_excel)

if nargin < 1
    data_directory  = 'source_data';
end
if nargin < 2
    read_from_excel = false;
end

alpha_resp = .001;  plot_visibility = 'on'; narrow_ylims = true;
seed = 12;          no_overlaps = false;    which_half = 0;
binned_test = true; 

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, read_from_excel, data_directory);