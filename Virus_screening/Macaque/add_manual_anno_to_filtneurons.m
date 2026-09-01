clc
clearvars
load('filt_neurons-clust-viralidx.mat'); 

annotation = readtable('annotation.csv');
%%
% Extract uids (3rd column in annotation.csv) 
uids = string(annotation{:, 3}); 

% Find indices of matching neurons
[is_match, idx_filt] = ismember(filt_neurons.uid, uids);

% Get column names from annotation (excluding UID column)
%annotation_fields = annotation.Properties.VariableNames;
%annotation_fields(3) = []; % Remove the UID column

% add annotation (select columns) to filt_neurons
is_cell = annotation.is_cell;  
filt_neurons.is_cell = nan(size(filt_neurons.id));
filt_neurons.is_cell(is_match) = is_cell(idx_filt(is_match));

%%
dominant_barcode_count = annotation.dominant_barcode_count;  
filt_neurons.dominant_barcode_count = nan(size(filt_neurons.id)); 
filt_neurons.dominant_barcode_count(is_match) = dominant_barcode_count(idx_filt(is_match));

%%
virus = annotation.virus;
filt_neurons.virus = cell(size(filt_neurons.id));
filt_neurons.virus(is_match) = virus(idx_filt(is_match));


%%
save('filt_neurons-clust-manual_anno.mat', 'filt_neurons');

%%
load('filt_neurons-clust-manual_anno.mat');
%%
% Filter for cells with virus 'HQ065L_8' and is_cell = 1
vir = strcmp(filt_neurons.virus, 'HQ065L_5_1');
is_cell = filt_neurons.is_cell == 1;
valid_cells = vir & is_cell;  

subclass_data = filt_neurons.subclass(valid_cells);

subclass_categorical = categorical(subclass_data);
subclass_counts = countcats(subclass_categorical);
subclass_labels = categories(subclass_categorical);

%figure;
%pie(subclass_counts, subclass_labels);
%title('Distribution of Subclasses for Cells with Virus HQ065L_12 and is_cell = 1');


%%
viruses_of_interest = {'HQ065L_8', 'HQ065L_5', 'HQ065L_5_1','HQ065L_12'};

% Initialize a cell array to store subclass counts for each virus
subclass_counts = cell(length(viruses_of_interest), 1);
for i = 1:length(viruses_of_interest)
    % Filter data for the current virus
    is_current_virus = strcmp(filt_neurons.virus, viruses_of_interest{i});
    current_subclasses = filt_neurons.subclass(is_current_virus);
    
    % Convert to categorical to count unique subclasses
    current_subclasses_categorical = categorical(current_subclasses);
    subclass_counts{i} = countcats(current_subclasses_categorical);
end
% Create a figure
figure;

% Hold the plot to overlay bars
hold on;

% Define colors for each virus
colors = lines(length(viruses_of_interest));

% Loop through each virus and plot its subclass distribution
for i = 1:length(viruses_of_interest)
    % Create a vector of zeros for all subclasses
    subclass_vector = zeros(1, length(unique_subclasses));
    
    % Map current subclass counts to the vector
    [~, idx] = ismember(filt_neurons.subclass(is_current_virus), unique_subclasses);
    subclass_vector(idx) = subclass_counts{i};
    
    % Plot the bar for the current virus
    bar(i, subclass_vector, 'FaceColor', colors(i, :), 'BarWidth', 0.8);
end

% Set x-tick labels to unique subclasses
set(gca, 'XTick', 1:length(unique_subclasses), 'XTickLabel', unique_subclasses);

% Rotate x-tick labels for better readability
xtickangle(45);

% Add labels and title
xlabel('Subclass');
ylabel('Count');
title('Subclass Distribution for Different Viruses');

% Add a legend
legend(viruses_of_interest, 'Location', 'Best');

% Release the hold on the plot
hold off;
