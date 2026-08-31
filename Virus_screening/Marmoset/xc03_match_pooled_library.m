clearvars
clc
load data/filt_neurons_combined_clust
%%
% output_dir=fullfile(pwd,'figs','pooled_lib_0mismatch');
output_dir=fullfile(pwd,'figs','pooled_lib_1mismatch');
mkdir(output_dir);
%%

lib=readtable('Virus_Barcode1.csv');
lib.Barcode=upper(lib.Barcode);


%% sanity check: where are the barcoded cells
% plot only the ones that are glu
subclass_to_include=find(~contains(filt_neurons.subclass_names,{'NA','gaba','Nonneuronal'}));

%
uniq_slice=unique(filt_neurons.slice);
for nn=1:numel(uniq_slice)
    in_slice=filt_neurons.slice==uniq_slice(nn);
    in_subclass=ismember(filt_neurons.subclass_idx,subclass_to_include);
    figure;
    hold on;
    %plot all cells
    scatter(filt_neurons.pos(in_slice&~filt_neurons.is_barcoded,2),...
        filt_neurons.pos(in_slice&~filt_neurons.is_barcoded,1),...
        5, ...
        [0.9 0.9 0.9], ...
        'filled');
    scatter(filt_neurons.pos(in_slice&filt_neurons.is_barcoded&in_subclass,2),...
        filt_neurons.pos(in_slice&filt_neurons.is_barcoded&in_subclass,1),...
        5, ...
        [0.2 0.2 0.2], ...
        'filled');
    pbaspect([range(xlim),range(ylim),1]);
    title("Slice"+nn)
    axis off
    exportgraphics(gcf,fullfile(output_dir,"GluBCcells_Slice"+uniq_slice(nn)+".pdf"), ...
        'ContentType','vector')
    
end

%% using the is_barcoded flag, find barcoded cells that match the library pool
%truncate library pool
lib_bc=cell2mat(cellfun(@(x) x(1:size(filt_neurons.soma_bc,2)),lib.Barcode,...
    'UniformOutput',false));


[~,lib_bc]=ismember(lib_bc,'GTAC');




d=pdist2(filt_neurons.soma_bc,lib_bc,'hamming')*size(filt_neurons.soma_bc,2);

[min_d,idx]=min(d,[],2);
% figure;histogram(min_d,0:10);
% matched=min_d<=0;
matched=min_d<=1;






% [~,idx]=ismember(filt_neurons.soma_bc,lib_bc,'rows');



filt_neurons.variant_idx=idx;
filt_neurons.variant_data=lib;

c=histcounts(idx(filt_neurons.is_barcoded&matched),0.5:1:size(lib,1)+0.5);
figure('Position',[50 50 300 250]);
bar(c, ...
    'k');
set(gca,'xtick',1:size(lib,1), ...
    'xticklabel',lib.Var1);
xlabel('Variant');
ylabel('# of cells');
exportgraphics(gcf,fullfile(output_dir,'variant_count.pdf'), ...
    'ContentType','vector');

%close all

%% for each variatn with label, what are the fraction of cell types that were labeled? use subclass-level labels
uniq_idx=unique(idx);
uniq_idx(uniq_idx==0)=[];

subclass_counts=zeros(numel(uniq_idx),numel(filt_neurons.subclass_names));



for nn=1:numel(uniq_idx)
    is_variant=filt_neurons.variant_idx==uniq_idx(nn);
    % figure('Position',[50 50 300 250]);
    subclass_counts(nn,:)=histcounts(filt_neurons.subclass_idx(is_variant&filt_neurons.is_barcoded&matched), ...
        0.5:1:double(max(filt_neurons.subclass_idx))+0.5);
    % bar(subclass_counts, ...
    % 'k');
    % set(gca,'xtick',1:numel(filt_neurons.subclass_names), ...
    %     'xticklabel',filt_neurons.subclass_names, ...
    %     "TickLabelInterpreter",'none');
    % xlabel('Labeled subclass');
    % ylabel('Cell count');
    % title(lib.Var1(uniq_idx(nn)));
    % exportgraphics(gcf,fullfile(output_dir,lib.Var1(uniq_idx(nn))+".pdf"));
    % close all


end


figure;
bar(subclass_counts,'stacked');
%
inh_types=(contains(filt_neurons.subclass_names,{'NA','gaba','Low','Cd'}));
exc_types=~inh_types;
rng(11)
cmap_exc=hclrainbow(sum(exc_types),[90 270]);
cmap_inh=[hclrainbow(floor(sum(inh_types)/2),[310 360]);hclrainbow(ceil(sum(inh_types)/2),[0 60])];

cmap_exc=cmap_exc(randperm(size(cmap_exc,1)),:);
cmap_inh=cmap_inh(randperm(size(cmap_inh,1)),:);

cmap=[cmap_inh(1:2,:);...
    cmap_exc(1:8,:);...
    cmap_inh(3:5,:); ...
    cmap_exc(9,:);...
    cmap_inh(6:end,:)];

% rng(11);
% cmap=cmap(randperm(size(cmap,1)),:);
colororder(cmap);
set(gca,'xtick',1:size(lib,1), ...
    'xticklabel',lib.Var1, ...
    'xticklabelrotation',90);
xlabel('Variant');
ylabel('# of cells');
legend(filt_neurons.subclass_names, ...
    'NumColumns',2, ...
    'Location','eastoutside', ...
    'Interpreter','none')
exportgraphics(gcf,fullfile(output_dir,'stacked_histogram.pdf'), ...
    'ContentType','vector');



%
subclass_to_flip=(contains(filt_neurons.subclass_names,{'NA','gaba','Low','Cd'}));
subclass_counts(:,subclass_to_flip)=-subclass_counts(:,subclass_to_flip);

figure;
bar(subclass_counts,'stacked');

colororder(cmap);
set(gca,'xtick',1:size(lib,1), ...
    'xticklabel',lib.Var1, ...
    'xticklabelrotation',90);
xlabel('Variant');
ylabel('# of cells');
legend(filt_neurons.subclass_names, ...
    'NumColumns',2, ...
    'Location','eastoutside', ...
    'Interpreter','none')
exportgraphics(gcf,fullfile(output_dir,'stacked_histogram_doublesided.pdf'), ...
    'ContentType','vector');

subclass_to_exclude=(contains(filt_neurons.subclass_names,{'NA','Low'}));
figure;
bar(subclass_counts(:,~subclass_to_exclude),'stacked');
colororder(cmap(~subclass_to_exclude,:));
set(gca,'xtick',1:size(lib,1), ...
    'xticklabel',lib.Var1, ...
    'xticklabelrotation',90);
xlabel('Variant');
ylabel('# of cells');
legend(filt_neurons.subclass_names(~subclass_to_exclude), ...
    'NumColumns',2, ...
    'Location','eastoutside', ...
    'Interpreter','none')
exportgraphics(gcf,fullfile(output_dir,'stacked_histogram_noNA.pdf'), ...
    'ContentType','vector');


%% Create CSV with variant counts by subclass
% Create a table with variant names and subclass counts
variant_names = lib.Var1(uniq_idx);

% Create table with proper column names
variant_subclass_table = array2table(subclass_counts, ...
    'VariableNames', filt_neurons.subclass_names);

% Add variant name as first column
variant_subclass_table = addvars(variant_subclass_table, variant_names, ...
    'Before', 1, 'NewVariableNames', 'Variant');

% Write to CSV
writetable(variant_subclass_table, fullfile(output_dir, 'variant_subclass_counts.csv'), ...
    'WriteRowNames', false);

disp(['Variant-subclass count table saved to: ' fullfile(output_dir, 'variant_subclass_counts.csv')]);














