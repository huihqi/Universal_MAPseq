output_dir=fullfile(pwd,"figs","match_barcodes");
mkdir(output_dir);



%%
load('filt_neurons-clust.mat');
%% Translate filt_neurons dom_bc code to MAPseq code
nucleotides='GTAC';
dom_bc=char(filt_neurons.soma_bc);
dom_bc(:)='N';
for nn=1:numel(filt_neurons.dom_bc)
    if ~isempty(filt_neurons.dom_bc{nn})
        dom_bc(nn,:)=nucleotides(filt_neurons.dom_bc{nn});
    end
end

%% Match to viral codebook (no mismatch)
load viral_codebook2.mat;

codes=char(codebook(:,2));

[~,viral_idx]=ismember(dom_bc,codes(:,1:size(dom_bc,2)),'rows');

sum(viral_idx>0)
filt_neurons.viral_idx=viral_idx;
save('filt_neurons-clust-viralidx.mat','filt_neurons','-v7.3');

%% which ones are barcoded? Add threshold
count_thresh=10;
complexity_thresh=0.2;
[~,c]=bccomplexity(filt_neurons.dom_bc);
is_bc=filt_neurons.dom_bc_count>=count_thresh&c>complexity_thresh; % unlike xc01, is_bc here combines both count threshold and complexity threshold

fprintf(['%u barcodes matched to projections.\n ' ...
    '%u barcodes in barcoded cells matched to projections.\n '], ...
    sum(viral_idx>0), ...
    sum(viral_idx>0&is_bc))

%% plot counts histogram
c=histcounts(viral_idx(is_bc),0.5:1:max(viral_idx)+0.5);
figure('Position',[50 50 250 200]);
bar(c, ...
    'k');
set(gca, ...
    'xtick',1:size(codebook,1), ...
    'xticklabel',codebook(:,1), ...
    'xticklabelrotation',90);
ylabel('Number of cells');
exportgraphics(gcf,fullfile(output_dir,'total_counts.pdf'), ...
    'ContentType','vector');
% plot by slice
uniq_slices=unique(filt_neurons.slice);
figure('Position',[50 50 450 300]);
for nn=1:numel(uniq_slices)
    in_slice=filt_neurons.slice==uniq_slices(nn);

    subplot(2,4,nn)
    c=histcounts(viral_idx(is_bc&in_slice),0.5:1:max(viral_idx)+0.5);
    bar(c, ...
        0.7,...
        'k');
    set(gca, ...
        'xtick',1:size(codebook,1), ...
        'xticklabel',codebook(:,1), ...
        'xticklabelrotation',90);
    ylabel('Number of cells');
    title("Slice "+uniq_slices(nn));

end
exportgraphics(gcf,fullfile(output_dir,'total_counts_byslice.pdf'), ...
    'ContentType','vector');

%close all;

%% plot on slices
output_dir = fullfile(pwd, "figs", "slice_imgs");
mkdir(output_dir);

downsample_rate = 0.3;
num_colors = max(max(viral_idx), 10); % Ensure enough colors
cmap = distinguishable_colors(num_colors, [1 1 1; 0 0 0]); % Avoid white/black

uniq_slices = unique(filt_neurons.slice);
downsample_toplot = rand(numel(filt_neurons.id), 1) <= downsample_rate;

for nn = 1:numel(uniq_slices)
    in_slice = filt_neurons.slice == uniq_slices(nn);

    figure('Position', [50, 50, 800, 600])
    scatter(filt_neurons.pos(in_slice & downsample_toplot, 1), ...
        filt_neurons.pos(in_slice & downsample_toplot, 2), ...
        5, [0.9, 0.9, 0.9], 'filled');
    hold on;

    for mm = 1:max(viral_idx)
        color_idx = mod(mm-1, size(cmap,1)) + 1; % Ensure index is within cmap range
        scatter(filt_neurons.pos(in_slice & (viral_idx == mm) & is_bc, 1), ...
            filt_neurons.pos(in_slice & (viral_idx == mm) & is_bc, 2), ...
            3, cmap(color_idx, :), 'filled');
    end

    pbaspect([range(xlim), range(ylim), 1]);
    set(gca, 'YDir', 'reverse');
    axis off;
    legend(["All"; string(codebook(:, 1))], 'Location', 'eastoutside')
    title("Slice " + uniq_slices(nn))

    exportgraphics(gcf, fullfile(output_dir, "Slice" + uniq_slices(nn) + ".pdf"), 'ContentType', 'vector');

    % Zoomed-in version
    if any(in_slice & (viral_idx > 0) & is_bc)
        for mm = 1:max(viral_idx)
            color_idx = mod(mm-1, size(cmap,1)) + 1; % Ensure index is within cmap range
            scatter(filt_neurons.pos(in_slice & (viral_idx == mm) & is_bc, 1), ...
                filt_neurons.pos(in_slice & (viral_idx == mm) & is_bc, 2), ...
                20, cmap(color_idx, :), 'filled');
        end

        xlim([min(filt_neurons.pos(in_slice & (viral_idx > 0) & is_bc, 1)) - 100, ...
            max(filt_neurons.pos(in_slice & (viral_idx > 0) & is_bc, 1)) + 100]);
        ylim([min(filt_neurons.pos(in_slice & (viral_idx > 0) & is_bc, 2)) - 100, ...
            max(filt_neurons.pos(in_slice & (viral_idx > 0) & is_bc, 2)) + 100]);
        pbaspect([range(xlim), range(ylim), 1]);

        exportgraphics(gcf, fullfile(output_dir, "Slice" + uniq_slices(nn) + "zoomin.pdf"), 'ContentType', 'vector');
    end

    close;
end

%% count by subclass
uniq_subclasses=unique(filt_neurons.subclass);
% uniq_subclasses

[~,subclass_idx]=ismember(filt_neurons.subclass,uniq_subclasses);
subclass_viral_count=zeros(numel(uniq_subclasses),size(codebook,1));
for nn=1:numel(uniq_subclasses)
    in_subclass=subclass_idx==nn;
    subclass_viral_count(nn,:)=histcounts(viral_idx(in_subclass&is_bc),0.5:1:size(codebook,1)+0.5);
end

figure;
imagesc(subclass_viral_count./sum(subclass_viral_count,2));
set(gca,'ytick',1:numel(uniq_subclasses), ...
    'yticklabel',uniq_subclasses, ...
    "TickLabelInterpreter",'none', ...
    'xtick',1:size(codebook,1), ...
    'xticklabel',codebook(:,1));

figure;
imagesc(subclass_viral_count./sum(subclass_viral_count));
set(gca,'ytick',1:numel(uniq_subclasses), ...
    'yticklabel',uniq_subclasses, ...
    "TickLabelInterpreter",'none', ...
    'xtick',1:size(codebook,1), ...
    'xticklabel',codebook(:,1));

to_keep=~ismember(codebook(:,1),'VEEV'); % plot after removing VEEV

figure;
imagesc(subclass_viral_count(:,to_keep)./sum(subclass_viral_count(:,to_keep),2));
set(gca,'ytick',1:numel(uniq_subclasses), ...
    'yticklabel',uniq_subclasses, ...
    "TickLabelInterpreter",'none', ...
    'xtick',1:size(codebook(to_keep,:),1), ...
    'xticklabel',codebook(to_keep,1));

%%
rng(100);
cmap=hclrainbow(12);
cmap=cmap(randperm(size(cmap,1)),:);

figure('Position',[50 50 250 200]);
h=bar(subclass_viral_count(1:10,to_keep)','stack');
for nn=1:numel(h)
h(nn).FaceColor=cmap(nn,:);
end
legend(uniq_subclasses(1:10), ...
    'Interpreter','none', ...
    'Location','eastoutside');
set(gca,'xtick',1:sum(to_keep), ...
    'xticklabel',codebook(to_keep,1), ...
    'xticklabelrotation',90);
ylabel('# of neurons');
% colormap(cmap)
exportgraphics(gcf,fullfile(output_dir,'celltype_dist.pdf'), ...
    'ContentType','vector');




%%
subclass_viral_count(:,end+1)=histcounts(subclass_idx,0.5:1:max(subclass_idx)+0.5);
cortical_exc=subclass_viral_count(1:10,[to_keep;true]);

figure;h=bar((cortical_exc./sum(cortical_exc))','stack');
for nn=1:numel(h)
h(nn).FaceColor=cmap(nn,:);
end
legend(uniq_subclasses(1:10), ...
    'Interpreter','none', ...
    'Location','eastoutside');
set(gca,'xtick',1:sum(to_keep)+1, ...
    'xticklabel',[codebook(to_keep,1);{'All'}]);
ylabel('Fraction of neurons');
exportgraphics(gcf,fullfile(output_dir,'celltype_fraction.pdf'), ...
    'ContentType','vector');
% colormap(cmap)

%%


%% check a few cells for each virus
fov = "MAX_Pos1_008_007";
output_dir = fullfile(pwd, 'figs', fov);
mkdir(output_dir);

[~, fov_idx] = ismember(fov, filt_neurons.fov_names);
in_fov = filt_neurons.fov == fov_idx;

data_folder = "C:\Users\maryam.majeed\Allen Institute\Barcoded Connectomics - General\Maryam\Analysis\20250121_MacaqueBARseq_alder_BB_frontal_7_maxprojection";

% Load Barcode Image
files = dir(fullfile(data_folder, fov, "RGB", "*bcseq*.jpg"));
if isempty(files)
    error("No bcseq images found for FOV: " + fov);
end
im_bc = imread(fullfile(files(1).folder, files(1).name));

% Load Gene Image
files = dir(fullfile(data_folder, fov, "RGB", "*hyb*.jpg"));
if isempty(files)
    error("No hyb images found for FOV: " + fov);
end
im_gene = imread(fullfile(files(1).folder, files(1).name));

% Load Cell Mask
S = load(fullfile(data_folder, fov, "aligned", 'cellmask.mat'));
im_mask = label2rgb(S.maski, 'jet', 'w', 'shuffle');

% Make a mask with viral identities
cellid = rem(filt_neurons.id(in_fov), 10000);
viral_idx_infov = viral_idx(in_fov) .* is_bc(in_fov); 

im_maskbyvirus = S.maski(:);
[~, I] = ismember(im_maskbyvirus, cellid);
im_maskbyvirus(I > 0) = viral_idx_infov(I(I > 0));
im_maskbyvirus(I == 0) = 0;
im_maskbyvirus = reshape(im_maskbyvirus, size(S.maski));

% Use a better colormap
num_colors = max(max(viral_idx), 10);
cmap = distinguishable_colors(num_colors, [1 1 1; 0 0 0]); % Avoid black/white
im_maskbyvirus_colored = label2rgb(im_maskbyvirus, cmap, 'w');

% Display & Save Images
figure; imshow(im_bc); title("Barcode Image");
figure; imshow(im_gene); title("Gene Image");
figure; imshow(im_mask); title("Cell Mask");
figure; imshow(im_maskbyvirus_colored); title("Mask by Virus");

imwrite(im_gene, fullfile(output_dir, 'gene.jpg'));
imwrite(im_bc, fullfile(output_dir, 'bc.jpg'));
imwrite(im_mask, fullfile(output_dir, 'mask.jpg'));
imwrite(im_maskbyvirus_colored, fullfile(output_dir, 'maskbyvirus.jpg'));



%imagesc(im_maskbyvirus)

% 
% fov_idx=
% in_fov=
% celloi=