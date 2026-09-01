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
num_barcodes_per_subclass = zeros(length(subclasses_to_plot), size(possible_bc, 1));

% Loop over each subclass and count the number of cells per barcode
for j = 1:length(subclasses_to_plot)
    % Find neurons of the current subclass
    subclass_idx = find(strcmp(filt_neurons.subclass, subclasses_to_plot{j}) & is_barcoded);
    
    % Count the number of cells per barcode for the current subclass
    for i = 1:size(possible_bc, 1)
        num_barcodes_per_subclass(j, i) = sum(ismember(filt_neurons.soma_bc(subclass_idx,:), possible_bc(i,:), 'rows'));
    end
end

% Plot the number of cells per barcode for each subclass
figure;
bar(num_barcodes_per_subclass, 'stacked');
set(gca, 'XTickLabel', {'L2/3 IT','L4/5 IT','L4 IT','L5 IT','L5/6 NP','L6 IT','Car3','L6 CT','L6b'});
ylabel('Number of Barcoded Neurons');
title('Barcode Distribution per Subclass');
legend('EEEV', 'VEEV', 'WEEV', 'SFV', 'SINV', 'Location', 'northeastoutside');
