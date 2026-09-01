clc
clearvars

% Load the data
filename = 'filt_neurons-clust'; % Change this if needed
load(filename)

% Define possible barcodes
possible_bc = [4,2,1;...
               3,3,2;...
               1,3,2;...
               4,1,3;...
               4,3,4];

% Specify the subclasses to plot
subclasses_to_plot = {'L2_3_IT','L4_5_IT','L4_IT','L5_IT','L5_6_NP','L6_IT','L6_IT_Car3','L6_CT','L6b'}; % Change these to your specific subclasses
%subclasses_to_plot = unique(filt_neurons.subclass);
% Identify indices of barcoded neurons
is_barcoded = filt_neurons.is_barcoded;
bc_idx = find(is_barcoded);

% Initialize arrays to hold the statistics
num_cells_per_subclass = zeros(size(possible_bc, 1), length(subclasses_to_plot));

% Loop over each possible barcode and count the number of cells per subclass
for i = 1:size(possible_bc, 1)
    % Find neurons with the current barcode
    matched_bc{i} = find(ismember(filt_neurons.soma_bc(is_barcoded,:), possible_bc(i,:), 'rows'));
    matched_idx{i} = bc_idx(matched_bc{i});
    
    % Count the number of cells per specified subclass
    for j = 1:length(subclasses_to_plot)
        num_cells_per_subclass(i, j) = sum(strcmp(filt_neurons.subclass(matched_idx{i}), subclasses_to_plot{j}));
    end
end

% Plot the number of cells per specified subclass for each barcode
figure;
bar(num_cells_per_subclass, 'stacked');
set(gca, 'XTickLabel', {'EEEV', 'VEEV', 'WEEV', 'SFV', 'SINV'});
ylabel('Number of Exc Neurons');
title('Subclass Distribution');
legend('L2/3 IT','L4/5 IT','L4 IT','L5 IT','L5/6 NP','L6 IT','Car3','L6 CT','L6b', 'Location', 'northeastoutside');
