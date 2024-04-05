function figureS1_loc(data_dir, id, spiketable)

%load data
if nargin < 1; data_dir = 'source_data'; end
if nargin < 2; id = 1; end
if nargin < 3
    load(sprintf('%s/Fig2_data.mat', data_dir), 'spiketable_loc');
    spiketable = spiketable_loc;
end

fs = 14;
factor_aspect_ratio = 1.32;

paper_size = [8,10.5];
h1 = figure('visible', 'on');
set(h1,'PaperUnits','centimeters');
set(h1, 'PaperSize', paper_size, 'PaperPosition', [0 0 paper_size(1) paper_size(2)]);
% 
% space_height = .01;
% raster_width = .5;
% thumbnail_width = .3;
% subject = spiketable.subject(id);
% sess = ids(id, 2);
% chan = ids(id, 3);
% clus = ids(id, 4);
reg = spiketable.brain_reg(id);
% reg_strs = {'Amy', 'Hipp', 'EC', 'PHC'};
reg_strs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
reg_str = reg_strs{reg};
% raster_bottom = .07;
% raster_height = .7;
% raster_width = .5;
%% create the background
width_square = .9;
height_square = width_square/factor_aspect_ratio;
bapos = [.05 .05/factor_aspect_ratio, .9, height_square];
blackAx = axes('position',bapos);

fill([0,1,1,0], [0,0,1,1], [.8, .8, .8], 'linestyle', 'none');
set(blackAx, 'yticklabel', []); 
set(blackAx, 'xticklabel', [])
blackAx.YAxis.Visible = 'off';
blackAx.XAxis.Visible = 'off';
theSpikes = spiketable.rel_ts{id};
valid_tr  = spiketable.successful_tap{id};
loc_order = spiketable.locs{id};
n_pos = length(unique(loc_order)); %always 9
sp_loc = cell(1, n_pos);
for i = 1:n_pos
    sp_loc{i} = theSpikes(loc_order==i & valid_tr);
end
% counter = sum(valid_tr);
hold on
%squares padding horizontal and vertical
s_pad_hor = .1/4;
s_pad_ver = s_pad_hor / factor_aspect_ratio;
%squares size horizontal and vertical
s_s_hor = .8/3;
s_s_ver = s_s_hor / factor_aspect_ratio;
histboxes =...
    [.05+s_pad_hor,               bapos(2)+s_pad_ver,               s_s_hor, s_s_ver/2;...%1 bottom left
     .05+s_pad_hor*2 + s_s_hor,   bapos(2)+s_pad_ver,               s_s_hor, s_s_ver/2;...%2
     .05+s_pad_hor*3 + s_s_hor*2, bapos(2)+s_pad_ver,               s_s_hor, s_s_ver/2;...%3
     .05+s_pad_hor,               bapos(2)+s_pad_ver*2 + s_s_ver,   s_s_hor, s_s_ver/2;...%4 center left
     .05+s_pad_hor*2 + s_s_hor,   bapos(2)+s_pad_ver*2 + s_s_ver,   s_s_hor, s_s_ver/2;...%5 
     .05+s_pad_hor*3 + s_s_hor*2, bapos(2)+s_pad_ver*2 + s_s_ver,   s_s_hor, s_s_ver/2;...%6 
     .05+s_pad_hor,               bapos(2)+s_pad_ver*3 + s_s_ver*2, s_s_hor, s_s_ver/2;...%7 top left
     .05+s_pad_hor*2 + s_s_hor,   bapos(2)+s_pad_ver*3 + s_s_ver*2, s_s_hor, s_s_ver/2;...%8
     .05+s_pad_hor*3 + s_s_hor*2, bapos(2)+s_pad_ver*3 + s_s_ver*2, s_s_hor, s_s_ver/2]; %9
rasterboxes = histboxes + repmat([0, s_s_ver/2, 0, 0], size(histboxes, 1), 1);
edges = -500:100:1500;
centers = edges(1:end-1)+mean(diff(edges))/2;
maxval = 0;
%% plot rasterplots and histograms
hh = gobjects(1, n_pos); % array of graphics objects 
for i = 1:n_pos
    ah = axes('position', rasterboxes(i,:));
    set(ah, 'color', 'w');
    ah.YAxis.Visible = 'off';
    %stim onset line
    plot([0,0], [0, length(sp_loc{i})], 'k')
    for tr = 1:length(sp_loc{i})
        ypos = [length(sp_loc{i})-tr+1; length(sp_loc{i})-tr];
        xpos = sp_loc{i}{tr};
        line([xpos; xpos], repmat(ypos,1, length(xpos)), 'color', 'k',...
            'linewidth', 1);
    end
    xlim([-500 1500])
    ylim([0 length(sp_loc{i})])
    set(ah, 'yticklabel', []); 
    xline(1000, 'color', [.7 .7 .7]);
    ah.YAxis.Visible = 'off';
        box off
        set(ah, 'xticklabel', [])
    hh(i) = axes('position', histboxes(i,:));
    histdata = histcounts(cell2mat(sp_loc{i}), edges)/...
        length(sp_loc{i})*1000/mean(diff(edges));
    plot([0 0], [0 1000], 'color', [.3 .3 .3])
    hold on
    bar(centers, histdata, 1, 'facecolor', [.3 .3 .3])
    xline(1000, 'color', [.7 .7 .7]);
    maxval = max(maxval, max(histdata));
    if i==1
        set(hh(i), 'fontsize', fs)
        xticks([0 1000]);
        box off
    else
        set(gca,'XColor','none') 
    end
    set(gca,'YColor', 'none')
    xlim([-500 1500])
end
for i =1:9
    ylim(hh(i), [0 maxval+1])
end

%% create density plot
aDense = axes('position',  [rasterboxes(9,1), s_s_ver*4, s_s_hor, s_s_ver*.8]);
plot_spikes = spiketable.spikeshapes{id};
density_plot(plot_spikes)
yl = double(ylim); xl = xlim;
yt = [ceil(yl(1)/10)*10, 0, floor(yl(2)/10)*10];
yticks(yt);
yspan = diff(yl);
xlabel('')
text(xl(2), double(yl(1)-yspan*.25), 'ms', 'fontsize', fs)
aDense.FontSize = fs;
%% text: brain region
axes('position', [rasterboxes(7,1), s_pad_ver*6 + s_s_ver*3,...
    s_s_hor*2, s_s_ver])
xlim([-1 1]); ylim([-1 1]);
text(0,0, reg_str, 'Fontsize', fs, 'horizontalalignment', 'center');
axis off

fname = sprintf('%s/figures/FigureS1_loc_%i.png', data_dir, id);
print(gcf, '-dpng', fname, '-r600');
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);
