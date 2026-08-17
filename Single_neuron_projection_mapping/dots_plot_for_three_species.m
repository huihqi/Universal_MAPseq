
combined_cluster_mats=cluster_mats;
fn = fieldnames(combined_cluster_mats);


last2_fields = fn(end-2:end);


result_mean = struct();
result_count_gt0 = struct();


result_mean_table = struct();
result_count_table = struct();

for f = 1:numel(last2_fields)
    this_field = last2_fields{f};
    

    ncol = [];
    for i = 1:numel(combined_cluster_mats)
        x = combined_cluster_mats(i).(this_field);
        if ~isempty(x)
            ncol = size(x, 2);
            break;
        end
    end
    

    if isempty(ncol)
        warning('%s empty', this_field);
        continue;
    end
    

    mean_mat = zeros(numel(combined_cluster_mats), ncol);
    count_mat = zeros(numel(combined_cluster_mats), ncol);
    

    for i = 1:numel(combined_cluster_mats)
        x = combined_cluster_mats(i).(this_field);
        
        if isempty(x)
            mean_mat(i, :) = 0;
            count_mat(i, :) = 0;
        else
            mean_mat(i, :) = mean(x, 1);
            count_mat(i, :) = sum(x > 0, 1);
        end
    end

    result_mean.(this_field) = mean_mat;
    result_count_gt0.(this_field) = count_mat;
    
  
    mean_varnames  = strcat("col", string(1:ncol), "_mean");
    count_varnames = strcat("col", string(1:ncol), "_count_gt0");
    
 
    result_mean_table.(this_field)  = array2table(mean_mat,  'VariableNames', cellstr(mean_varnames));
    result_count_table.(this_field) = array2table(count_mat, 'VariableNames', cellstr(count_varnames));
    

    result_mean_table.(this_field)  = addvars(result_mean_table.(this_field),  (1:numel(combined_cluster_mats))', 'Before', 1, 'NewVariableNames', 'cluster_id');
    result_count_table.(this_field) = addvars(result_count_table.(this_field), (1:numel(combined_cluster_mats))', 'Before', 1, 'NewVariableNames', 'cluster_id');
end


tbl_mean_MAPseq1_datanorm   = result_mean_table.MAPseq1_datanorm;
tbl_count_MAPseq1_datanorm  = result_count_table.MAPseq1_datanorm;

tbl_mean_MAPseq2_datanorm   = result_mean_table.MAPseq2_datanorm;
tbl_count_MAPseq2_datanorm  = result_count_table.MAPseq2_datanorm;

tbl_mean_MAPseq3_datanorm   = result_mean_table.MAPseq3_datanorm;
tbl_count_MAPseq3_datanorm  = result_count_table.MAPseq3_datanorm;
tbl_mean_MAPseq3_datanorm{:,[2 8:13]} = 0;
tbl_count_MAPseq3_datanorm{:,[2 8:13]} = 0;

% new_order = [1 2 8 7 3 4 9 6 5];
%new_order = [1 2 3 4 5 6];
%new_order = [8 4 5 6 3 7 2 1 9];  %for 9 clusters
new_order = [15 11 13 12 10 14 6 7 8 3 2 1 9 4 5];  %for 15 clusters
% parrot
%tbl_mean_MAPseq1_datanorm1  = tbl_mean_MAPseq1_datanorm(new_order, :);
%tbl_count_MAPseq1_datanorm1 = tbl_count_MAPseq1_datanorm(new_order, :);

%mouse
%tbl_mean_MAPseq1_datanorm1  = tbl_mean_MAPseq2_datanorm(new_order, :);
%tbl_count_MAPseq1_datanorm1 = tbl_count_MAPseq2_datanorm(new_order, :);

%zebra finch
tbl_mean_MAPseq1_datanorm1  = tbl_mean_MAPseq3_datanorm(new_order, :);
tbl_count_MAPseq1_datanorm1 = tbl_count_MAPseq3_datanorm(new_order, :);



mean_data  = tbl_mean_MAPseq1_datanorm1{:, 2:end};
count_data = tbl_count_MAPseq1_datanorm1{:, 2:end};

% counts
row_sum = sum(count_data, 2);


count_norm = zeros(size(count_data));
for i = 1:size(count_data, 1)
    if row_sum(i) > 0
        count_norm(i, :) = count_data(i, :) / row_sum(i);
    end
end


[nrow, ncol] = size(mean_data);
[xg, yg] = meshgrid(1:ncol, 1:nrow);

x = xg(:);
y = yg(:);
c = mean_data(:);        
s = count_norm(:);       
count_raw = count_data(:);

mask = (c ~= 0) & (count_raw ~= 0);

x = x(mask);
y = y(mask);
c = c(mask);
s = s(mask);


size_min = 1;
size_max = 300;
s_plot = size_min + s * (size_max - size_min);


figure('Color','w');
ax1 = axes('Position',[0.10 0.12 0.62 0.78]);   
scatter(ax1, x, y, s_plot, c, 'filled');
colormap(ax1, parula);

cb = colorbar(ax1);
cb.Label.String = 'Mean value';

ylim(ax1, [0.5, nrow + 0.5]);
xlim(ax1, [0.5, ncol + 0.5]);

set(ax1, 'XTick', 1:ncol);
set(ax1, 'YTick', 1:nrow);

xlabels = tbl_mean_MAPseq1_datanorm.Properties.VariableNames(2:end);
set(ax1, 'XTickLabel', xlabels, 'XTickLabelRotation', 45);

ylabels = string(tbl_mean_MAPseq1_datanorm1.cluster_id);
set(ax1, 'YTickLabel', ylabels);

xlabel(ax1, 'Target');
ylabel(ax1, 'Cluster');
title(ax1, 'Dot plot: color = mean, size = normalized counts');
set(ax1, 'YDir', 'reverse');
box(ax1, 'on');


cb.Position = [0.74 0.18 0.02 0.66];


q = [0.25 0.50 0.75 1.00];
s_legend = size_min + q * (size_max - size_min);


ax2 = axes('Position',[0.80 0.28 0.15 0.45]);
axis(ax2, [0 1 0 1]);
ax2.Visible = 'off';
hold(ax2, 'on');

x0 = 0.25;
y0 = [0.85 0.62 0.39 0.16];

scatter(ax2, repmat(x0,1,4), y0, s_legend, 'k', 'filled');

for i = 1:4
    text(ax2, 0.48, y0(i), sprintf('%.2f', q(i)), ...
        'VerticalAlignment', 'middle', 'FontSize', 10);
end

text(ax2, 0.00, 0.98, 'Normalized count', ...
    'FontWeight', 'bold', 'FontSize', 10, ...
    'VerticalAlignment', 'top');