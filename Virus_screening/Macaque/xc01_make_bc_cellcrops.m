

load filt_neurons-viralidx.mat % needs to have uid, batch_num, viral_idx, is_bc
%%
crop_rad=60; % crop size is crop_rad*2+1
% lbound=1; % how big are the images?
% ubound=3200; % how big are the images?
spacer_width=5; % spacer bar width between crops

output_dir=fullfile(pwd,'figs','cell_crops');
%data_dir = fullfile(pwd, 'processed');
data_dir=fullfile('\\allen\programs\celltypes\workgroups\barcoded_connectomics\20250121_MacaqueBARseq_alder_BB_frontal_7_maxprojection\','processed');
temp_dir=fullfile(pwd,'temp');


mkdir(output_dir)

mkdir(temp_dir);

%% Get the list of cells to plot. Plot cells that are is_bc and has a valid viral_idx
filt_neurons.batch_num=ones(size(filt_neurons.id));
filt_neurons=add_uid(filt_neurons);
cell_idx=find(filt_neurons.is_bc&filt_neurons.viral_idx>0);
numel(cell_idx)



fov=filt_neurons.fov_names(filt_neurons.fov(cell_idx));
pixelregion=[filt_neurons.pos40x(cell_idx,1)+[0,crop_rad*2],...
        filt_neurons.pos40x(cell_idx,2)+[0, crop_rad*2]];
% pixelregion(pixelregion<=lbound)=lbound;
% pixelregion(pixelregion>=ubound)=ubound;
% pixelregion1=mat2cell(pixelregion, ones(size(pixelregion,1),1),[2 2]); %
% if tif files
%how many digits is needed?
dig=ceil(log10(numel(cell_idx)));
file_idx=cellfun(@(x) sprintf("%0"+dig+"d",x),num2cell(1:numel(cell_idx)), ...
    'UniformOutput',false)';

fname=fullfile(output_dir,file_idx+"_uid_"+filt_neurons.uid(cell_idx)+".jpg");


spacer=round(0.9*ones(crop_rad*2+1,spacer_width,3)*255);
%%        
currfov={};
for nn=1:numel(cell_idx)
    if ~ismember(fov(nn),currfov)% check if this is an fov that has worked on previously. If not, get images
        currfov=fov{nn};
        [im_gene,im_hyb,im_bc1,im_bc2,mask_rgb]=get_imgs(currfov,data_dir,temp_dir,crop_rad);
    end
    % assemble cropped image
    assembled_img=assemble_crop_imgs(im_gene,im_hyb,im_bc1,im_bc2,mask_rgb,spacer,pixelregion(nn,:));
    % writeimage
    imwrite(assembled_img,fname{nn},'quality',100);

end



%% write csv

[~,file_name]=fileparts(fname);
uid=filt_neurons.uid(cell_idx);
is_cell=cell(numel(uid),1);
index=(1:numel(cell_idx))';
dominant_barcode=cell2mat(filt_neurons.dom_bc(cell_idx));
ntcode='GTACN';
dominant_barcode=ntcode(dominant_barcode);
dominant_barcode_count=filt_neurons.dom_bc_count(cell_idx);
virus=filt_neurons.viral_codebook(filt_neurons.viral_idx(cell_idx),1);

T=table(index,file_name,uid,dominant_barcode,dominant_barcode_count,virus,is_cell);
% check if output file exists, in case it overwrites a previous annotation
csv_fname=fullfile(output_dir,'annotation.csv');
if ~isfile(csv_fname)
    writetable(T,csv_fname);
else
    X=input('Annotation file exists. Overwrite?(Y/N)','s');
    if ismember(upper(X),'N')
        warning('Annotation file exists and not overwritten, check existing file.')
    elseif ismember(upper(X),'Y')
        warning('Annotation file overwritten.');

        writetable(T,csv_fname);
    end
end


%%
function assembled_img=assemble_crop_imgs(im_gene,im_hyb,im_bc1,im_bc2,mask_rgb,spacer,thispixelregion)
thispixelregion=uint16(thispixelregion);
im_gene_crop=im_gene(thispixelregion(3):thispixelregion(4), ...
    thispixelregion(1):thispixelregion(2),:);
im_hyb_crop=im_hyb(thispixelregion(3):thispixelregion(4), ...
    thispixelregion(1):thispixelregion(2),:);
im_bc1_crop=im_bc1(thispixelregion(3):thispixelregion(4), ...
    thispixelregion(1):thispixelregion(2),:);
im_bc2_crop=im_bc2(thispixelregion(3):thispixelregion(4), ...
    thispixelregion(1):thispixelregion(2),:);
mask_rgb_crop=mask_rgb(thispixelregion(3):thispixelregion(4), ...
    thispixelregion(1):thispixelregion(2),:);
assembled_img=[im_gene_crop,spacer,...
    im_hyb_crop,spacer,...
    im_bc1_crop,spacer,...
    im_bc2_crop,spacer,...
    mask_rgb_crop];
end


 function [im_gene,im_hyb,im_bc1,im_bc2,mask_rgb]=get_imgs(currfov,data_dir,temp_dir, crop_rad)
     fov_rgbpath=fullfile(data_dir, currfov, 'RGB');
 
     % im_gene=imread(fullfile(fov_rgbpath,'RGBalignedfixedn2vgeneseq01.jpg'));
     % im_hyb=imread(fullfile(fov_rgbpath,'RGBalignedn2vhyb01.jpg'));
     % im_bc1=imread(fullfile(fov_rgbpath,'RGBalignedfixedregn2vbcseq01.jpg'));
     % im_bc2=imread(fullfile(fov_rgbpath,'RGBalignedfixedregn2vbcseq02.jpg'));
     im_gene=padarray(imread(fullfile(fov_rgbpath,'RGBalignedfixedn2vgeneseq01.tif')), ...
         [crop_rad,crop_rad,0]);
     im_hyb=padarray(imread(fullfile(fov_rgbpath,'RGBalignedn2vhyb01.tif')), ...
         [crop_rad,crop_rad,0]);
     im_bc1=padarray(imread(fullfile(fov_rgbpath,'RGBalignedfixedregn2vbcseq01.tif')), ...
         [crop_rad,crop_rad,0]);
     im_bc2=padarray(imread(fullfile(fov_rgbpath,'RGBalignedfixedregn2vbcseq02.tif')), ...
         [crop_rad,crop_rad,0]);
 
     %% figure;imshow(im_gene)
     % extract cellmask, look for aligned.zip
     aligned_folder = fullfile(data_dir, currfov, 'aligned');
     % mask_fidx = dir(fullfile(aligned_folder, 'cellmask.mat'));
     % if isempty(mask_fidx)
     %    error('cellmask.mat not found in the aligned folder: %s', aligned_folder);
     % end
     S = load(fullfile(aligned_folder,'cellmask.mat')); 
     %aligned_files=unzip(fullfile(data_dir, currfov, 'aligned.zip'),temp_dir);
    
     % find cellmask.mat
     %mask_fidx=find(contains(aligned_files,'cellmask.mat'));
     %S=load(aligned_files{mask_fidx});
     rng(10);
     mask_rgb=padarray(label2rgb(S.maski,'jet','k','shuffle'), ...
         [crop_rad,crop_rad,0]);
     %% cleanup
     % for mm=1:numel(aligned_files)
     %     delete(aligned_files{mm})
     % end
     % rmdir(fullfile(temp_dir,currfov),'s');
 end