%function gm_raster_vertical(data_dir, id, item_or_loc, spiketable)
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function gm_raster_vertical(data_dir, id, item_or_loc, spiketable)

%load data
if nargin < 1; data_dir = 'source_data'; end

load(fullfile(data_dir, 'Fig2_data.mat'),...
    'spiketable_item', 'spiketable_loc');

% cond_strs = {'enc', 'rec'};
% root_dir = '/media/data';
% if ~exist('spikecell', 'var')
%     load(sprintf('%s/datacell_%s.mat', root_dir, cond_strs{which_part}))
% end

if nargin < 1
    id = 1;
end
if nargin < 2
    item_or_loc = 1;
end
% if spiketimes are not passed as argument, load data
if nargin < 3
    load(sprintf('%s/Fig2_data.mat', data_dir),...
            'spiketable_item', 'spiketable_loc');
    if item_or_loc == 1 % load item neuron data
        spiketable = spiketable_item;
    else % load location neuron data
        spiketable = spiketable_loc;
    end
end

% set up figure
fs = 14;
paper_size = [5 21];
h1 = figure('visible', 'off');
set(h1,'PaperUnits','centimeters');
set(h1, 'PaperSize', paper_size, 'PaperPosition', [0 0 paper_size(1) paper_size(2)]);

subject = spiketable.subject(id);
sess = ids(id, 2);
chan = ids(id, 3);
clus = ids(id, 4);
reg = spikecell{id, 1};
% reg_strs = {'Amy', 'Hipp', 'EC', 'PHC'};
reg_strs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
reg_str = reg_strs{reg};
raster_bottom = .07;
raster_height = .7;
raster_width = .5;
%% create the giant rasterplot
mainAx = axes('position',[.4,raster_bottom,raster_width,raster_height]);
mainAx.YAxis.Visible = 'off';
theSpikes = spikecell{id, 13};
img_order = spikecell{id, 8};
valid_trials = spikecell{id, 9};
n_img = length(unique(img_order));
sp_img = cell(1, n_img);
jpg_path = cell(1, n_img);
for i = 1:n_img
    sp_img{i} = theSpikes(img_order==i & valid_trials);
    the_img = spikecell{id,12}{find(img_order==i, 1)};
    garbage = dir(sprintf('%s/%.3i*gm%i/response_img/%s',...
        root_dir, subject, sess, the_img));
    jpg_path{i} = sprintf('%s/%s', garbage.folder, garbage.name);
end
counter = sum(valid_trials);
hold on
plot([0 0], [0 counter], 'k')
for i = 1:n_img
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
ylim([0 sum(valid_trials)])
xline(1000, 'color', [.7 .7 .7]);
xlabel('ms')
%going with fs-1 here because somehow that font got printed bigger
set(gca, 'yticklabel', [], 'fontsize', fs-1)
%% put  in img thumbnails
thumb_padding = raster_height*.2/(n_img);
thumb_height = raster_height*.8/n_img;
for i=1:n_img
    axes('position', [.05, raster_bottom+thumb_padding*(n_img-i+.5)+thumb_height*(n_img-i),...
        .25, thumb_height])
    axis off
    imshow(jpg_path{i})
end
spikes_dir = dir(sprintf('%s/%.3i*gm%i/times_CSC%i.mat',...
    root_dir, subject, sess, chan));
load(fullfile(spikes_dir.folder, spikes_dir.name))
%% create density plot
axes('position', [.4, .84, raster_width, .09])
plot_spikes = spikes(cluster_class(:,1)==clus, :);
density_plot(plot_spikes)
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
fname = sprintf('%s/gm_resps/rasterCompact_%.3igm%i_ch%i_cl%i_%i', root_dir,...
    subject, sess, chan, clus, id);
print(gcf, '-dpng', fname, '-r600');
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);
