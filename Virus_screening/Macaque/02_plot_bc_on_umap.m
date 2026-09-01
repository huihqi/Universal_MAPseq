% Clear workspace and load necessary data
clearvars
clc
load filt_neurons-clust.mat

% Setup output directory
output_dir = fullfile('figs', 'bc_umap');
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end
c=table2cell(readtable('analysis/glu/umap.csv','NumHeaderLines',0));
umap_coor=cell2mat(c);
    qc=~(contains(filt_neurons.clustid,'filtered')|...
        contains(filt_neurons.clustid,'gaba')|...
        contains(filt_neurons.clustid,'Glu')|...
        contains(filt_neurons.clustid,'subCTX')|...
        contains(filt_neurons.clustid,'SubCTX'));

is_barcoded = filt_neurons.is_barcoded(qc) == 1;

%indataset = filt_neurons.clustid == qc;
uniq_types=unique(filt_neurons.subclass(qc));
[~,clustidx]=ismember(filt_neurons.subclass,uniq_types);
filt_neurons.clustidx = clustidx;
uniq_clustidx=unique(clustidx);

% Define a custom colormap with darker shades of blue
dark_blue = [0 0 0.5]; % Dark blue color
map_dark_blue = repmat(dark_blue, numel(uniq_types), 1);

% Initialize figure
figure('Position', [50 50 800 400]);
hold on;

% Plot all cells in grey
scatter(umap_coor(:, 1), umap_coor(:, 2), 5, [0.6 0.6 0.6], 'filled', ...
    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05);
xl = xlim;
yl = ylim;
pbaspect([range(xl), range(yl), 1]);
axis off

% Plot barcoded cells in dark blue
hold on;
for i = 1:numel(uniq_clustidx)
    scatter(umap_coor(is_barcoded, 1), umap_coor(is_barcoded, 2), 5, dark_blue, 'filled', ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
end

% Add legend
legend({'all cells', 'barcoded'}, 'Location', 'eastoutside', 'NumColumns', 2, 'Interpreter', 'none')
axis off

% Save the figure
exportgraphics(gcf, fullfile(output_dir, ['cluster_glu_umap', '.tif']), 'ContentType', 'Image', 'Resolution', 300);

% Close the figure
close all;


% Define colormap
%cmap = distinguishable_colors(numel(uniq_types));
%rng(10)
%map1 = colormap('lines');

% Initialize figure
%figure('Position', [50 50 800 400]);
%hold on;

% Plot all cells in grey
%scatter(umap_coor(:, 1), umap_coor(:, 2), 5, [0.6 0.6 0.6], 'filled', ...
%    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05);
%xl = xlim;
%yl = ylim;
%pbaspect([range(xl), range(yl), 1]);
%axis off

% Plot barcoded cells in colors
%hold on;
%for i = 1:numel(uniq_clustidx)
%    scatter(umap_coor(is_barcoded, 1), umap_coor(is_barcoded, 2), 5, map1(i, :), 'filled', ...
%        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
%end
%legend({'all cells', 'barcoded'}, 'Location', 'eastoutside', 'NumColumns', 2, 'Interpreter', 'none')
%axis off

% Save the figure
%exportgraphics(gcf, fullfile(output_dir, ['cluster_whole_umap', '.tif']), 'ContentType', 'Image', 'Resolution', 300);

% Close the figure
%close all;
