function figureS1_item(data_dir, id, spiketable)

% this function creates full rasterplots
% parameters
% data_dir (str):     path to source_data directory
% id (int):           1 or 2, since there are 2 units of each type
%                     (responding to items or locations)
% spiketable (table): contains spiketimes and trial information. Different
%                     tables have to be used for item and location
%                     responses.


if nargin < 1; data_dir = 'source_data'; end
if nargin < 2
    id = 1;
end
% if spiketimes are not passed as argument, load data
if nargin < 3
    load(sprintf('%s/Fig2_data.mat', data_dir), 'spiketable_item');
    spiketable = spiketable_item;
end

% set up figure
fs = 14;
paper_size = [5 21];
h1 = figure('visible', 'on');
set(h1,'PaperUnits','centimeters');
set(h1, 'PaperSize', paper_size, 'PaperPosition', [0 0 paper_size(1) paper_size(2)]);

reg = spiketable.brain_reg(id);
reg_strs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
reg_str = reg_strs{reg};
raster_bottom = .07;
raster_height = .7;
raster_width = .5;

%% create the giant rasterplot
mainAx = axes('position',[.4,raster_bottom,raster_width,raster_height]);
mainAx.YAxis.Visible = 'off';
theSpikes = spiketable.rel_ts{id};
stim_order = spiketable.item_index{id};
valid_trials = spiketable.successful_tap{id};
n_stims = length(unique(stim_order));
sp_img = cell(1, n_stims);
jpg_path = cell(1, n_stims);
image_folder = strrep(data_dir, 'source_data', 'images');
for i = 1:n_stims
    sp_img{i} = theSpikes(stim_order==i & valid_trials);
    image_fname = spiketable.fname{id}{find(stim_order==i, 1)};
    jpg_path{i} = sprintf('%s/%s', image_folder, image_fname);
end
counter = sum(valid_trials);
hold on
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
ylim([0 sum(valid_trials)])
xline(1000, 'color', [.7 .7 .7]);
xlabel('ms')
set(gca, 'yticklabel', [], 'fontsize', fs-1)

%% put  in image thumbnails
thumb_padding = raster_height*.2/(n_stims);
thumb_height = raster_height*.8/n_stims;
for i=1:n_stims
    axes('position', [.05, raster_bottom+thumb_padding*(n_stims-i+.5)+thumb_height*(n_stims-i),...
        .25, thumb_height])
    axis off
    imshow(jpg_path{i})
end

%% create density plot
axes('position', [.4, .84, raster_width, .09])
plot_spikes = spiketable.spikeshapes{id};
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

%% print figures
fname = sprintf('%s/figures/FigureS1_item_%i.mat', data_dir, id);
print(gcf, '-dpng', fname, '-r600');
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);
