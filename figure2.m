
function [results] = figure2(data_dir, id, item_or_loc, spiketable)

% [results] = figure2(data_dir, id, item_or_loc, spiketable)
%
% this function will plot one of the 4 panels as seen in figure 2
% 
% INPUTS:
%      data_dir (str): path to source data
%            id (int): row number within spiketable of unit to plot
%    img_or_pos (int): 1: plot image response, 2: plot location response
%  spiketable (table): optoinal. Contains the experimental data and spike times    

image_directory = fullfile(strrep(data_dir, 'source_data', ''), 'images');

if nargin < 1
    data_dir  = 'source_data';
end
if nargin < 2
    id = 1;
end
if nargin < 3
    item_or_loc = 1;
end

results = struct;

% if spiketimes are not passed as argument, load data
if nargin < 4
    if item_or_loc == 1 % load item neuron data
        spiketable = load(sprintf('%s/Fig2_data.mat', data_dir),...
            'spiketable_item');
    else % load location neuron data
        spiketable = load(sprintf('%s/Fig2_data.mat', data_dir),...
            'spiketable_loc');
    end
end

stim_types = {'item', 'loc'};
stim_type = stim_types{item_or_loc};

%% configure figure
fs = 16; %font size
h1 = figure('visible', 'on');
set(h1, 'PaperUnits','centimeters');
set(h1, 'PaperSize', [10 13], 'PaperPosition', [0 0 10 13]);

if item_or_loc ==1
    % find all response eliciting items
    prefersThis = spiketable.pref_items{id};
    stimuli = spiketable.item_index{id};
    % thumbnail
    f_struct = dir(sprintf('%s/%s', image_directory,...
        spiketable.fname{id}{find(stimuli==prefersThis, 1)}));
    nstims   = length(unique(stimuli));
    st_c = 'item_index'; % spiketable column for image indices
else
    % find all response eliciting locations
    prefersThis = spiketable.pref_locs{id};
    sign_str = strrep(strjoin(string(prefersThis)), ' ', '');
    f_struct = dir(sprintf('%s/pos%s.png', image_directory, sign_str));
    if isempty(f_struct)
        fprintf('please create %s/pos%s.png\n', image_directory, sign_str) 
    end
    st_c = 'locs'; % spiketable column for position
    nstims = 9;
end

% get stimulus jpg
pref_jpg = fullfile(f_struct.folder, f_struct.name);

% some parameters for figure setup
padding     = .11; % white frame around everything
pl_h        = (1 - 4 * padding) / 3; % we want 3 rows of plots
full_w      = 1 - 2.8 * padding; % for raster plot and convolved plot
thumbnail_w = (full_w - 1.5 * padding) * (1/3); %thumbnail width
dens_w      = (full_w - 1.5 * padding) * (2/3); %density plot width

% stimulus image
axes('Position', [padding, 3 * padding + 2 * pl_h, thumbnail_w, pl_h]);
imshow(pref_jpg)
% spike shape
axes('Position', [1.5 * padding + 1.5 * padding + thumbnail_w,...
    3 * padding + 2 * pl_h, dens_w, pl_h]);

%% density plot
spike_f = sprintf('%s/density_%s_%i.mat', data_dir, stim_type, id);
load(spike_f, 'spikeshapes')
density_plot(spikeshapes)
xticks([-.5, 0 .5, 1, 1.34])
xticklabels({'-.5', '0', '.5', '1', 'ms'})
xlabel('')
set(gca, 'Fontsize', fs)

%% raster plot
axes('Position', [1.5 * padding, 2* padding + pl_h, full_w, pl_h])
[rel_ts_r, rel_ts_f] = plot_resp_byID(id, item_or_loc, spiketable);
ylabel('trials'); %xlabel('ms');
set(gca, 'Fontsize', fs)
xlim([-500 1500])


%% now the lines / spiking rates
axes('Position', [1.5*padding, padding, full_w, pl_h]);
hold on
binsize     = 100;
edgesEven   = -700:binsize:1700;
edgesOdd    = -650:binsize:1650;
allEdges    = -700:binsize/2:1700;
maxval = 0;
[solidBlue, solidRed, dashedBlue, dashedRed] = deal([]);

%for subsequent memory test
mem_test_win = [0 1500];
n_spikes_c = [];
n_spikes_f = [];
% Define a function to filter the values within each array
filter_and_length_func = ...
    @(x) length(x(x >= mem_test_win(1) & x <= mem_test_win(2)));

for m = 1:nstims % e.g. 1:8 different stimuli or 1:9 positions
    hold on
    rmbrd_all{2,m} = [];
    false_all{2,m} = [];
    rmbrd_all{1,m} = find(spiketable.(st_c){id} == m & ...
            spiketable.successful_tap{id} &  spiketable.subsequent_memory{id});
    false_all{1,m} = find(spiketable.(st_c){id} == m & ...
            spiketable.successful_tap{id} & ~spiketable.subsequent_memory{id});
    ts_relC = spiketable.rel_ts{id}(rmbrd_all{1,m});
    ts_relF = spiketable.rel_ts{id}(false_all{1,m});
    rmbrd_all{2,m} = ts_relC;
    false_all{2,m} = ts_relF;
    N_c_even = histc(cell2mat(rmbrd_all{2,m}), edgesEven);
    N_c_odd  = histc(cell2mat(rmbrd_all{2,m}), edgesOdd);
    N_f_even = histc(cell2mat(false_all{2,m}), edgesEven);
    N_f_odd  = histc(cell2mat(false_all{2,m}), edgesOdd);
    N_c = zeros(1,length(allEdges));
    N_c(1:2:end) = N_c_even; N_c(2:2:end) = N_c_odd;
    %we have to do this binning for individual traces now, to be able to do
    %a cluster permutation test.
    if size(N_f_even, 2)==0
        N_f_even = zeros(size(N_f_even, 1), 1);
    end
    if size(N_f_odd, 2) == 0
        N_f_odd = zeros(size(N_f_odd, 1), 1);
    end
    N_f = zeros(1,length(allEdges));
    N_f(1:2:end) = N_f_even; N_f(2:2:end) = N_f_odd;
    %if m==prefersThis{id} || ismember(m, sign_stims)
    if any(m==prefersThis)% || ismember(m, sign_stims)
        %change this from indiviual lines per stimulus to one avg line for
        %all preferred and one (dashed) line for all non pref
        solidBlue = [solidBlue; N_c/(binsize/1000)/length(rmbrd_all{1,m})];
        solidRed = [solidRed; N_f/(binsize/1000)/length(false_all{1,m})];
        % keep track of spike counts for subsequent memory test
        % subsequent memory is compared in a time window of mem_test_win,
        % which is 0-1500ms post stim onset
        % Apply the function to each array in ts_relC using cellfun
        filtered_lengths_c = cellfun(filter_and_length_func, ts_relC);
        filtered_lengths_f = cellfun(filter_and_length_func, ts_relF);
        n_spikes_c = [n_spikes_c, filtered_lengths_c];
        n_spikes_f = [n_spikes_f, filtered_lengths_f];
    else
        dashedBlue = [dashedBlue; N_c/(binsize/1000)/length(rmbrd_all{1,m})];
        dashedRed  = [dashedRed; N_f/(binsize/1000)/length(false_all{1,m})];
    end
    maxval = max([maxval, max(N_c), max(N_f)]);
end

%perform memory test: ranksum test
%set method to approximate in order to obtain stats.zval every time. in one
%case, one n is 6 such that by default it would be set to 'exact'
[p, ~, stats] = ranksum(n_spikes_c, n_spikes_f,...
    'method', 'approximate','tail', 'right');
stim_types = {'item', 'loc'};
stim_type  = stim_types{item_or_loc};

fprintf('%s unit %i: subs rem vs forg, p=%.3g, z=%.3g\n',...
    stim_type, id, p, stats.zval);

fprintf('    n spikes fired within 1500ms: remembered %i, forgotten %i\n',...
    length(n_spikes_c), length(n_spikes_f));
hclnp= plot(allEdges, mean(dashedBlue,1), '--b');
hflnp= plot(allEdges, mean(dashedRed ,1), '--r');
hclp = plot(allEdges, mean(solidBlue, 1), 'b', 'linewidth', 1.5);
hflp = plot(allEdges, mean(solidRed,  1), 'r', 'linewidth', 1.5);
line([0 0],[-1.7 maxval*1.2], 'Color', 'k'); % stim presentation time
% ylim([0 maxval]);
xlim([-500 1500])
if id == 1 && item_or_loc == 1
    yls = ylim; ylim([-.05*yls(2) 55])
elseif id == 2 && item_or_loc == 1
    yls = ylim; ylim([-.05*yls(2) 20])
elseif id == 1 && item_or_loc == 2
    yls = ylim; ylim([-.05*yls(2) 18])
elseif id == 2 && item_or_loc == 2
    yls = ylim; ylim([-.05*yls(2) 20])
else
    yls = ylim; ylim([-.05*yls(2) yls(2)])
end

set(gca, 'Fontsize', fs)
lgnd = legend([hclp, hflp], {'remembered', 'forgotten'}, 'Location', 'northeast',...
    'numcolumns', 1, 'fontsize', fs-2);
legend('boxoff')
lgnd.Position =[.475 .2 .5 .1];
xlabel('ms')
ylabel('Hz')

%fill output struct results
results.stim_type     = stim_type;
results.pref_jpg      = pref_jpg;
results.x             = allEdges;
results.rel_ts_rem    = rel_ts_r;
results.rel_ts_forg   = rel_ts_f;
results.pref_rem      = mean(solidBlue, 1);
results.pref_forg     = mean(solidRed, 1);
results.non_pref_rem  = mean(dashedBlue, 1);
results.non_pref_forg = mean(dashedRed, 1);
results.spikefile     = spike_f;

% save output
fname = sprintf('%s/figures/fig2_%s_unit%i.png', ...
        data_dir, stim_type, id);
print(gcf, '-dpng', fname,'-r300')
set(gcf, 'renderer', 'painters');
print(gcf, '-depsc', '-r600', [fname(1:end-3), 'eps']);
