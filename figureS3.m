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