clc
clearvars

% Load the dataset
filename = 'filt_neurons-clust'; % Change this to your actual filename
load(filename)

% Define possible barcodes to be matched
possible_bc = [4,2,1;...
               3,3,2;...
               1,3,2;...
               4,1,3;...
               4,3,4];

% Specify the subclasses to plot
subclasses_to_plot = {'L2_3_IT','L4_5_IT','L4_IT','L6_CT','L6b'}; % Change these to your specific subclasses

% Identify indices of barcoded and non-barcoded neurons
is_barcoded = filt_neurons.is_barcoded;
bc_idx = find(is_barcoded);
non_bc_idx = find(~is_barcoded);

% Loop over each barcode type
for i = 1:size(possible_bc, 1)
    % Initialize arrays to hold the statistics for the current barcode
    med_read_counts_barcoded = zeros(length(subclasses_to_plot), 1);
    med_read_counts_non_barcoded = zeros(length(subclasses_to_plot), 1);

    % Loop over each subclass and calculate the median reads per neuron for each barcode
    for j = 1:length(subclasses_to_plot)
        % Find indices of barcoded and non-barcoded neurons of the current subclass
        subclass_bc_idx = find(strcmp(filt_neurons.subclass, subclasses_to_plot{j}) & is_barcoded);
        subclass_non_bc_idx = find(strcmp(filt_neurons.subclass, subclasses_to_plot{j}) & ~is_barcoded);

        % Find neurons with the current barcode
        matched_bc_idx = subclass_bc_idx(ismember(filt_neurons.soma_bc(subclass_bc_idx,:), possible_bc(i,:), 'rows'));

        % Calculate the median reads per barcoded neuron
        if ~isempty(matched_bc_idx)
            med_read_counts_barcoded(j) = median(sum(filt_neurons.expmat(matched_bc_idx, :), 2));
        end

        % Calculate the median reads per non-barcoded neuron
        if ~isempty(subclass_non_bc_idx)
            med_read_counts_non_barcoded(j) = median(sum(filt_neurons.expmat(subclass_non_bc_idx, :), 2));
        end
    end

    % Plot the comparison for the current barcode
    figure;
    bar([med_read_counts_barcoded, med_read_counts_non_barcoded]);
    set(gca, 'XTickLabel', subclasses_to_plot);
    ylabel('Median Number of Reads per Neuron');
    title(['Median Reads per Neuron: Barcode ', num2str(possible_bc(i,:))]);
    legend('Barcoded', 'Non-Barcoded', 'Location', 'northeastoutside');
end
