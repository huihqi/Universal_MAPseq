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
%% 

% Identify indices of barcoded and non-barcoded neurons
is_barcoded = filt_neurons.is_barcoded;
bc_idx = find(is_barcoded);
non_bc_idx = find(~is_barcoded);

% Initialize arrays to hold the statistics
avg_reads_per_barcoded_neuron = zeros(length(subclasses_to_plot), size(possible_bc, 1));
avg_reads_per_non_barcoded_neuron = zeros(length(subclasses_to_plot), 1);

% Loop over each subclass and calculate the average reads per neuron for each barcode
for j = 1:length(subclasses_to_plot)
    % Find indices of barcoded and non-barcoded neurons of the current subclass
    subclass_bc_idx = find(strcmp(filt_neurons.subclass, subclasses_to_plot{j}) & is_barcoded);
    subclass_non_bc_idx = find(strcmp(filt_neurons.subclass, subclasses_to_plot{j}) & ~is_barcoded);

    % Loop over each barcode
    for i = 1:size(possible_bc, 1)
        % Find neurons with the current barcode
        matched_bc_idx = subclass_bc_idx(ismember(filt_neurons.soma_bc(subclass_bc_idx,:), possible_bc(i,:), 'rows'));

        % Calculate the average reads per barcoded neuron
        if ~isempty(matched_bc_idx)
            avg_reads_per_barcoded_neuron(j, i) = mean(filt_neurons.expmat(matched_bc_idx));
        end
    end

    % Select a random sample of non-barcoded neurons with the same size as barcoded neurons
    num_barcoded_neurons = length(subclass_bc_idx);
    if num_barcoded_neurons > 0
        if length(subclass_non_bc_idx) >= num_barcoded_neurons
            random_non_bc_idx = randsample(subclass_non_bc_idx, num_barcoded_neurons);
        else
            random_non_bc_idx = subclass_non_bc_idx; % If not enough non-barcoded neurons, use all of them
        end
        
        % Calculate the average reads per non-barcoded neuron
        avg_reads_per_non_barcoded_neuron(j) = mean(filt_neurons.expmat(random_non_bc_idx));
    end
end

% Plot the comparison
figure;
bar(avg_reads_per_barcoded_neuron, 'stacked');
hold on;
plot(avg_reads_per_non_barcoded_neuron, 'k--o', 'LineWidth', 2, 'MarkerSize', 6);
set(gca, 'XTickLabel', subclasses_to_plot);
ylabel('Average Number of Reads per Neuron');
title('Average Reads per Neuron: Barcoded vs Non-Barcoded');
legend(['EEEV', 'VEEV', 'WEEV', 'SFV', 'SINV', 'Non-Barcoded'], 'Location', 'northeastoutside');
hold off;
