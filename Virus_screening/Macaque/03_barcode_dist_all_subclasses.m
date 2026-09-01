clc
clearvars

% Load the data
filename = 'filt_neurons-clust-manual_anno.mat'; % Change this if needed
load(filename)

% Define possible barcodes
possible_bc = [4,2,1;...
               3,3,2;...
               1,3,2;...
               4,1,3;...
               4,3,4];

% Identify indices of barcoded neurons
is_barcoded = filt_neurons.is_barcoded;
bc_idx = find(is_barcoded);

% Get unique subclasses from the data
subclasses_to_plot = unique(filt_neurons.subclass);

% Initialize arrays to hold the statistics
num_cells_per_subclass = zeros(size(possible_bc, 1), length(subclasses_to_plot));

% Loop over each possible barcode and count the number of cells per subclass
for i = 1:size(possible_bc, 1)
    % Find neurons with the current barcode
    matched_bc{i} = find(ismember(filt_neurons.soma_bc(is_barcoded,:), possible_bc(i,:), 'rows'));
    matched_idx{i} = bc_idx(matched_bc{i});
    
    % Count the number of cells per unique subclass
    for j = 1:length(subclasses_to_plot)
        num_cells_per_subclass(i, j) = sum(strcmp(filt_neurons.subclass(matched_idx{i}), subclasses_to_plot{j}));
    end
end

% Plot the number of cells per unique subclass for each barcode
figure;
bar_handle = bar(num_cells_per_subclass, 'stacked');

% Create a colormap with unique colors for each subclass and shuffle the colors
num_colors = length(subclasses_to_plot);
cmap = jet(num_colors); % Change 'jet' to any other colormap if desired
shuffled_cmap = cmap(randperm(num_colors), :);

% Apply the shuffled colormap to the bar handles
for k = 1:num_colors
    bar_handle(k).FaceColor = 'flat';
    bar_handle(k).CData = shuffled_cmap(k, :);
end

set(gca, 'XTickLabel', {'EEEV', 'VEEV', 'WEEV', 'SFV', 'SINV'});
ylabel('Number of Neurons');
title('Subclass Distribution');
legend(subclasses_to_plot, 'Location', 'northeastoutside');
