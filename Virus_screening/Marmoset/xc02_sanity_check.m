clearvars
clc
load(fullfile("data","filt_neurons_combined_clust.mat"));

%%
output_dir=fullfile(pwd,"figs","slice_images");
mkdir(output_dir);

%% color code cells by cell types on slices
% remove NA
plot_cells_on_slices(filt_neurons,output_dir,"subclass");
plot_cells_on_slices(filt_neurons, output_dir, "clust");
%%
plot_bccells_on_slices(filt_neurons,output_dir,"subclass");

%% since barcoded, visualize where barcoded cells are, and how many of them 
function plot_bccells_on_slices(filt_neurons,output_dir,level)
clr_factor=5; % this determines how faded non-barcoded cells are



output_dir=fullfile(output_dir,level+"bc");
mkdir(output_dir);

%level="subclass";

group_names=filt_neurons.(level+"_names");
group_idx=filt_neurons.(level+"_idx");

idx_toremove=find(ismember(group_names,"NA")|...
    contains(lower(group_names),"low_quality")|...
    contains(lower(group_names),"low_qual"));

uniq_idx=unique(group_idx);
cmap=hclrainbow(numel(uniq_idx)+1,[0 300]);
toplot=~ismember(group_idx,idx_toremove);
cmap(ismember(uniq_idx,idx_toremove),:)=0.8;
c=cmap(group_idx,:);




uniq_brains=unique(filt_neurons.brain_name);

for bb=1:numel(uniq_brains)
    in_brain=ismember(filt_neurons.brain_name,uniq_brains(bb));
    uniq_slices=unique(filt_neurons.slice(in_brain));
    for nn=1:numel(uniq_slices)
        in_slice=filt_neurons.slice==uniq_slices(nn);

        figure('Position',[50 50 1200 800])
        scatter(filt_neurons.pos(in_slice&in_brain&toplot,1), ...
            filt_neurons.pos(in_slice&in_brain&toplot,2), ...
            1, ...
            (c(in_slice&in_brain&toplot,:)+clr_factor)./(clr_factor+1),...
            'filled');
        set(gca,'ydir','reverse');
        hold on;
        scatter(filt_neurons.pos(in_slice&in_brain&filt_neurons.is_barcoded,1), ...
            filt_neurons.pos(in_slice&in_brain&filt_neurons.is_barcoded,2), ...
            1, ...
            c(in_slice&in_brain&filt_neurons.is_barcoded,:),...
            'filled');
        set(gca,'ydir','reverse');
        axis off;
        pbaspect([range(xlim),range(ylim),1])
        title(uniq_brains(bb)+", Slice"+uniq_slices(nn), ...
            'Interpreter','none')
        fname=level+"bc"+"_"+uniq_brains(bb)+"_Slice"+uniq_slices(nn)+".jpg";
        exportgraphics(gcf,fullfile(output_dir,fname), ...
            'Resolution',1200);
        close
    end
end

figure('Position',[50 50 1200 800])
hold on
for mm=1:numel(uniq_idx)
    if ~ismember(uniq_idx(mm),idx_toremove)
        scatter([],[],5,cmap(uniq_idx(mm),:),'filled');
    end
end
axis off
legend(group_names(~ismember(1:numel(group_names),idx_toremove)), ...
    'interpreter','none');
exportgraphics(gcf,fullfile(output_dir,level+"bc"+"_legend.jpg"), ...
    'Resolution',1200);
close
end




function plot_cells_on_slices(filt_neurons,output_dir,level)
output_dir=fullfile(output_dir,level);
mkdir(output_dir);

%level="subclass";

group_names=filt_neurons.(level+"_names");
group_idx=filt_neurons.(level+"_idx");

idx_toremove=find(ismember(group_names,"NA")|...
    contains(lower(group_names),"low_quality")|...
    contains(lower(group_names),"low_qual"));

uniq_idx=unique(group_idx);
cmap=hclrainbow(numel(uniq_idx)+1,[0 300]);
toplot=~ismember(group_idx,idx_toremove);

c=cmap(group_idx,:);




uniq_brains=unique(filt_neurons.brain_name);

for bb=1:numel(uniq_brains)
    in_brain=ismember(filt_neurons.brain_name,uniq_brains(bb));
    uniq_slices=unique(filt_neurons.slice(in_brain));
    for nn=1:numel(uniq_slices)
        in_slice=filt_neurons.slice==uniq_slices(nn);

        figure('Position',[50 50 1200 800])
        scatter(filt_neurons.pos(in_slice&in_brain&toplot,1), ...
            filt_neurons.pos(in_slice&in_brain&toplot,2), ...
            1, ...
            c(in_slice&in_brain&toplot,:),...
            'filled');
        set(gca,'ydir','reverse');
        axis off;
        pbaspect([range(xlim),range(ylim),1])
        title(uniq_brains(bb)+", Slice"+uniq_slices(nn), ...
            'Interpreter','none')
        fname=level+"_"+uniq_brains(bb)+"_Slice"+uniq_slices(nn)+".jpg";
        exportgraphics(gcf,fullfile(output_dir,fname), ...
            'Resolution',300);
        close
    end
end

figure('Position',[50 50 1200 800])
hold on
for mm=1:numel(uniq_idx)
    if ~ismember(uniq_idx(mm),idx_toremove)
        scatter([],[],5,cmap(uniq_idx(mm),:),'filled');
    end
end
axis off
legend(group_names(~ismember(1:numel(group_names),idx_toremove)), ...
    'interpreter','none');
exportgraphics(gcf,fullfile(output_dir,level+"_legend.jpg"), ...
    'Resolution',300);
close
end