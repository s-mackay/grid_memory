function gm_response_permutation(nperm, seed)

if ~exist('nperm', 'var')
    nperm = 20;
end

randn('state', seed);

bsl_period = [-500 0];
response_period = [0 1000];
binsize = 100;
pcrit = .001; % alpha for binw ranksum
rootdir = '/media/data';
load(sprintf('%s/spiketable_enc.mat', rootdir))
continuous_sess_index = nan(size(spiketable.subject));
sess_counter = 0;
%add a session index to spiketable, iterate over cells to add it
for i =1:size(spiketable, 1)
    if i==1 || spiketable.subject(i) ~= spiketable.subject(i-1) || ...
               spiketable.session(i) ~= spiketable.session(i-1)
        sess_counter = sess_counter+1;
    end
    continuous_sess_index(i) = sess_counter;
end

% for total responsive cell counts, columns are
% amy, hipp, ec, phc, total
resp_unit_counts = zeros(nperm, 5);
has_item_response = false(size(spiketable, 1), nperm);
has_loc_response  = false(size(spiketable, 1), nperm);

% create an array with the number of trials in each session
n_trials_sess = zeros(max(continuous_sess_index), 1);
for i = 1:length(n_trials_sess)
    table_index = find(continuous_sess_index==i, 1);
    n_trials_sess(i) = length(spiketable.locs{table_index});
end
tic
for i = 1:nperm
    % create set of shuffled indices
    shuff_cell = cell(max(continuous_sess_index), 1);
    for sess = 1:length(shuff_cell)
        shuff_cell{sess} = randperm(n_trials_sess(sess));
    end   
    % iterate over cells
    for ii = 1:size(spiketable, 1)
        sess_index = continuous_sess_index(ii);
        %deal with items
        n_items = length(unique(spiketable.item_index{ii}));
        n_locs  = length(unique(spiketable.locs{ii})); % always 9
        trial_order = shuff_cell{continuous_sess_index(ii)};
        items       = spiketable.item_index{ii}(trial_order);
        locs        = spiketable.locs{ii}(trial_order);
        %valid_tap   = spiketable.successful_tap{ii}(trial_order);
        %subs_mem    = spiketable.subsequent_memory{ii}(trial_order);
        spikes4bsl = spiketable.rel_ts{ii};
        distribution = [];
        %iterate over items, calculate responses
        for iii = 1:n_items
            trials = items == iii;
            theSpikes = spiketable.rel_ts{ii}(trials);
            %trials_corr = trial_order == iii & valid_tap & subs_mem;
            %trials_false = trial_order == iii & valid_tap & ~subs_mem;
            [pval, consider, distribution] = binwise_ranksum(theSpikes, ...
                spikes4bsl, bsl_period, response_period, binsize,...
                0, 1, distribution, 0);
            if pval < pcrit && consider > 0
                has_item_response(ii, i) = true;
                %if only one item elicits a response, it's enough for us to
                %call this unit responsive and we don't have to go on
                %testing with other items
                continue
            end
        end
        for iii = 1:n_locs
            trials = locs == iii;
            theSpikes = spiketable.rel_ts{ii}(trials);
            %trials_corr = trial_order == iii & valid_tap & subs_mem;
            %trials_false = trial_order == iii & valid_tap & ~subs_mem;
            [pval, consider, distribution] = binwise_ranksum(theSpikes, ...
                spikes4bsl, bsl_period, response_period, binsize,...
                0, 1, distribution, 0);
            if pval < pcrit && consider > 0
                has_loc_response(ii, i) = true;
                %if only one location elicits a response, it's enough for 
                %us to call this unit responsive and we don't have to go on
                %testing with other locations
                continue
            end
        end
    end
    if mod(i, 10)==0
        disp(i)
        toc
    end
end

fname = sprintf("labelshuf%iperm_seed%i.mat", nperm, seed);
save(fname, 'has_item_response', 'has_loc_response');
keyboard
%{
pos_n = [24, 26, 1, 53];
stim_n = [80, 81, 32, 96];
brainregs = {'Amy', 'Hipp', 'EC', 'PHC'};
figure
for i =1:4
    units_reg = spiketable.brain_reg == i;
    %item
    distr = sum(has_item_response(units_reg, :), 1);
    n_larger = sum(distr > stim_n(i));
    n_equal = sum(distr== stim_n(i));
    emp_size = (n_larger + .5*n_equal) / length(distr);
    subplot(2,4,i)
    edges = 0:max(distr)+1;
    bar((edges(1:end-1) + .5)/sum(units_reg), histcounts(distr, edges));
    hold on
    xline(stim_n(i)/sum(units_reg), 'r')
    disp([i, emp_size]);
    units_reg = spiketable.brain_reg == i;
    title(sprintf('%s item cells shuffled', brainregs{i}));
    if i==1
        legend({'1220 shuffles', 'measured value'}, 'location', 'southeast')
    end
    xl = xlim;
    yl = ylim;
    text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n= %.2g',...
        emp_size), 'horizontalalignment', 'right');
    xlim(xl); ylim(yl);
    %pos
    distr = sum(has_loc_response(units_reg, :), 1);
    n_larger = sum(distr > pos_n(i));
    n_equal = sum(distr== pos_n(i));
    emp_size = (n_larger + .5*n_equal) / length(distr);
    subplot(2,4,i+4)
    edges = 0:max(distr)+1;
    bar((edges(1:end-1)+ .5)/sum(units_reg), histcounts(distr, edges));
    hold on
    xline(pos_n(i)/sum(units_reg), 'r')
    disp([i, emp_size]);
    title(sprintf('%s location cells shuffled', brainregs{i}));
    xlim(xl);
    ylim(yl);
    text(xl(2)*.5, yl(2)*.9, sprintf('empirical size\n= %.2g',...
        emp_size));
end
%}
            
        