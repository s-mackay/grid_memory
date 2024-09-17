% [relSpikesC, relSpikesF] = plot_resp_byID(id, ...
%     img_or_pos, spiketable)
%
% plots simple raster plot of preferred stimulus only
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function [relSpikesC, relSpikesF] = plot_resp_byID(id, ...
    img_or_pos, spiketable)

if img_or_pos ==1
    %prefThings = spiketable.pref_items{id};
    st_c   = 'item_index'; % spikecell column for image indices
    pref_c = 'pref_items'; % spikecell column for preferred images
else
    %prefThings = spiketable.pref_locs{id};
    st_c   = 'locs';      % spikecell column for position
    pref_c = 'pref_locs'; % spikecell column for preferred locs
end

prefThings = spiketable.(pref_c){id};

[relSpikesC, relSpikesF] = deal(cell(0,1));
%will increase this in the loop 
n_trials_c_f = zeros(length(prefThings), 2);
for i = 1:length(prefThings)
    %get indices of valid and subsequently correctly remembered trials
    tr_indsC = spiketable.(st_c){id}==prefThings(i) & ...
        spiketable.successful_tap{id} &  spiketable.subsequent_memory{id};
    %get indices of valid and subsequently not remembered (False) trials
    tr_indsF = spiketable.(st_c){id}==prefThings(i) & ...
        spiketable.successful_tap{id} & ~spiketable.subsequent_memory{id};
    relSpikesC = [relSpikesC, spiketable.rel_ts{id}(tr_indsC)];
    relSpikesF = [relSpikesF, spiketable.rel_ts{id}(tr_indsF)];
    n_trials_c_f(i,:) = [sum(tr_indsC), sum(tr_indsF)];
end

hold on

% plot red lines, one trial, i.e. row, at a time
for trF = 1:length(relSpikesF)
    ypos = length(relSpikesF)- trF + 1;
    xvals = [relSpikesF{trF}; relSpikesF{trF}];
    yvals = repmat([ypos-0.5; ypos+0.5], 1, length(relSpikesF{trF}));
    line(xvals, yvals, 'LineWidth', 0.8, 'Color', 'r');
end

% plot horizontal grey lines across the whole panel
lines_ypos = [length(relSpikesF)+length(relSpikesC)-[cumsum(n_trials_c_f(:,1));...
    sum(n_trials_c_f(:,1))+cumsum(n_trials_c_f(:,2))]]+.5;
xvals = repmat([-500; 1500], [1, numel(lines_ypos)-1]);
line(xvals, repmat(lines_ypos(1:end-1)', 2, 1), 'color', [.5,.5,.5])

% plot blue lines, one trial, i.e. row, at a time
for trC = 1:length(relSpikesC)
    ypos = trF+length(relSpikesC)-trC+1;
    xvals = [relSpikesC{trC}; relSpikesC{trC}];
    yvals = repmat([ypos-0.5; ypos+0.5], 1, length(relSpikesC{trC}));
    line(xvals, yvals, 'LineWidth', 0.8, 'Color', 'b');
end

%plot black vertical line at x=0
line([0 0],[-1.7 length(relSpikesC)+length(relSpikesF)+1.5], 'Color', 'k');
ylim([0, length(relSpikesC)+length(relSpikesF)+1])

xlim([-500 1500])


