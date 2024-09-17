% function [results] = figureS1_item(data_dir, id, data, read_from_xlsx)
% this function creates full rasterplots for item cells
%
% INPUTS:
%       data_dir (str): path to source_data directory
%             id (int): 1 or 2, since there are 2 units of each type
%                       (responding to items or locations)
%
% data (table or cell): the classic option is spiketable, which is loaded
%                       from .mat files. alternatively the struct called
%                       results.
% read_from_xlsx(bool): was the data loaded from .xlsx? (or .mat)
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function [results] = figureS1_item(data_dir, id, data, read_from_xlsx)

% set source directory for images and density data
if nargin < 1; data_dir = 'source_data'; end
if nargin < 2
    id = 1;
end
% if spiketimes are not passed as argument, load data
if nargin < 3
    load(sprintf('%s/Fig2_data.mat', data_dir), 'spiketable_item');
    spiketable = spiketable_item;
end
if nargin < 4
    read_from_xlsx = false;
end

% set up figure
fs = 14;
paper_size = [5 21];
h1 = figure('visible', 'on');
set(h1,'PaperUnits','centimeters');
set(h1, 'PaperSize', paper_size, 'PaperPosition',...
    [0 0 paper_size(1) paper_size(2)]);
set(h1, 'Units', 'pixels', 'Position', [0 0 300 1024]);

% set some variables depending on whether we're loading from .mat or .xlsx
if read_from_xlsx
    results = data;
    reg_str = results.brainreg;
else
    spiketable = data;
    reg = spiketable.brain_reg(id);
    reg_strs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
    reg_str = reg_strs{reg};
    results = struct; % for output
end

raster_bottom = .07;
raster_height = .7;
raster_width = .5;

%% create the giant rasterplot
mainAx = axes('position',[.4,raster_bottom,raster_width,raster_height]);
mainAx.YAxis.Visible = 'off';
if read_from_xlsx
    % get spike times
    fn = fieldnames(results);
    sp_fn = fn(startsWith(fn, 'rel_ts'));
    sp_img = cell(size(sp_fn));
    for i = 1:length(sp_fn)
        sp_img{i} = results.(sp_fn{i});
    end
    n_stims = length(sp_img);
    n_trials = length([sp_img{:}]); % number of trials
    jpg_paths = results.all_jpg;
else
    theSpikes = spiketable.rel_ts{id};
    stim_order = spiketable.item_index{id};
    valid_trials = spiketable.successful_tap{id};
    n_stims = length(unique(stim_order));
    sp_img = cell(1, n_stims);
    jpg_paths = cell(1, n_stims);
    image_folder = strrep(data_dir, 'source_data', 'images');
    for i = 1:n_stims
        fieldname = sprintf('rel_ts_stim%i', i);
        sp_img{i} = theSpikes(stim_order==i & valid_trials);
        results.(fieldname) = sp_img{i}; 
        image_fname = spiketable.fname{id}{find(stim_order==i, 1)};
        jpg_paths{i} = sprintf('%s/%s', image_folder, image_fname);
        jpg_paths{i} = strrep(jpg_paths{i}, 'patient_person',...
            'patient_person_written');
    end
    n_trials = sum(valid_trials);
end
hold on
counter = n_trials;
plot([0 0], [0 counter], 'k')
for i = 1:n_stims
    plot([-500 1500], [counter, counter], 'color', [.5 .5 .5]);
    for tr = 1:length(sp_img{i})
        ypos = [counter; counter-1];
        xpos = sp_img{i}{tr};
        line([xpos; xpos], repmat(ypos,1, length(xpos)), 'color', 'k');
        counter = counter-1;
    end

end
assert(counter==0);
xlim([-500 1500])
ylim([0, n_trials])
xline(1000, 'color', [.7 .7 .7]);
xlabel('ms')
set(gca, 'yticklabel', [], 'fontsize', fs-1)

%% put  in image thumbnails
thumb_padding = raster_height*.2/(n_stims);
thumb_height = raster_height*.8/n_stims;
for i=1:n_stims
    axes('position', [.05, raster_bottom+thumb_padding*(n_stims-i+.5)+...
        thumb_height*(n_stims-i), .25, thumb_height])
    axis off
    imshow(jpg_paths{i})
end

%% create density plot
axes('position', [.4, .84, raster_width, .09])
if read_from_xlsx
    load(results.spikefile, 'spikeshapes');
else
    spikeshapes = spiketable.spikeshapes{id};
end
density_plot(spikeshapes)
set(gca, 'fontsize', fs)
ylims = double(ylim);
text(.5, ylims(1)-(ylims(2)-ylims(1))*.4, 'ms', 'fontsize', fs,...
    'horizontalalignment', 'center')
xlabel('')

%% text: brain region
axes('position', [.05, .945, .9, .05], 'fontweight', 'bold')
xlim([-1 1]); ylim([-1 1]);
text(0,0, reg_str, 'Fontsize', fs, 'horizontalalignment', 'center');
axis off

%% print figures
fname = sprintf('%s/figures/FigureS1_item_%i.mat', data_dir, id);
print(gcf, '-dpng', fname, '-r600');
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);

%% complete output 
%fill output struct results
spike_f = sprintf('%s/density_item_%i.mat', data_dir, id);
results.all_jpg       = jpg_paths;
results.brainreg      = reg_str;
results.spikefile     = spike_f;