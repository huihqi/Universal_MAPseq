clearvars
clc
load filt_neurons_combined_clust
%%
output_dir=fullfile(pwd,'figs','pooled_lib');
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

[~,idx]=ismember(filt_neurons.soma_bc,lib_bc,'rows');

filt_neurons.variant_idx=idx;
filt_neurons.variant_data=lib;

c=histcounts(idx(filt_neurons.is_barcoded),0.5:1:size(lib,1)+0.5);
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
for nn=1:numel(uniq_idx)
    is_variant=filt_neurons.variant_idx==uniq_idx(nn);
    figure('Position',[50 50 300 250]);
    c1=histcounts(filt_neurons.subclass_idx(is_variant&filt_neurons.is_barcoded), ...
        0.5:1:double(max(filt_neurons.subclass_idx))+0.5);
    bar(c1, ...
    'k');
    set(gca,'xtick',1:numel(filt_neurons.subclass_names), ...
        'xticklabel',filt_neurons.subclass_names, ...
        "TickLabelInterpreter",'none');
    xlabel('Labeled subclass');
    ylabel('Cell count');
    title(lib.Var1(uniq_idx(nn)));
    exportgraphics(gcf,fullfile(output_dir,lib.Var1(uniq_idx(nn))+".pdf"));
    close all


end







