% change these manually
clearvars
clc

fname={fullfile(pwd,'..','reprocessed data','20250328_marmoset_AltoCtx_BunnyLH_3_2_3-3_Angela1_maxprojection','filt_neurons_ctx_fixedvglut2.mat')};


brain_name={'AltoCtx_20250328'};

batches=[6];

%%
% load each filt_neurons, combine them
for nn=1:length(fname)

    S=load(fname{nn});
    S.filt_neurons=S.filt_neurons_ctx;
    % for old datasets, add uid first
    if ~isfield(S.filt_neurons,'uid')
        S.filt_neurons=add_uid(S.filt_neurons);
    end
    % change batch numbers
    S.filt_neurons.batch_num(:)=uint8(batches(nn));

    % %for 0317_BunnyLH, update slice numbers. %% These don't match, verify slice numbers and fix
    % %this later.
    % if nn==4
    %      old_slices=1:16;
    %      new_slices=[71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146];
    %      [~,I]=ismember(S.filt_neurons.slice, old_slices);
    %      S.filt_neurons.slice(I>0)=new_slices(I(I>0));
    % end


    % for some datasets, reset pos on each slice.
    uniq_slices=unique(S.filt_neurons.slice);
    % only process cells in split slices
    uniq_slices(isnan(uniq_slices))=[];
    uniq_slices(uniq_slices==0)=[];

    for mm=1:numel(uniq_slices)
        in_slice=S.filt_neurons.slice==uniq_slices(mm);
        S.filt_neurons.pos(in_slice,:)=S.filt_neurons.pos(in_slice,:)-min(S.filt_neurons.pos(in_slice,:),1)+1;
    end

    % % for some datasets that have previously been clustered, remove the
    % % clusters
    % if isfield(S.filt_neurons,'clust_idx')
    %     S.filt_neurons=rmfield(S.filt_neurons,{'clust_idx','clust_names','subclass_idx','subclass_names'});
    % end
    % 
    % % for this dataset some batches had more genes than others. Only keep
    % % the genes that are shared.
    % if nn>1
    %     common_genes=intersect(S.filt_neurons.genes(:,1),filt_neurons.genes(:,1));
    %     missing_genes_dataset1=~ismember(filt_neurons.genes(:,1),common_genes);
    %     filt_neurons.genes(missing_genes_dataset1,:)=[];
    %     filt_neurons.expmat(:,missing_genes_dataset1)=[];
    % 
    %     [~,sorted_idx_dataset2]=ismember(filt_neurons.genes(:,1),S.filt_neurons.genes(:,1));
    %     S.filt_neurons.genes=S.filt_neurons.genes(sorted_idx_dataset2(sorted_idx_dataset2>0),:);
    %     S.filt_neurons.expmat=S.filt_neurons.expmat(:,sorted_idx_dataset2(sorted_idx_dataset2>0));
    % 
    % end

    % add brain name and update uid
    S.filt_neurons.brain_name=repmat(string(brain_name{nn}),numel(S.filt_neurons.uid),1);
    % for 0317_Alto, switch slice176-180 to Alto_LH
    % slice_range=176:180;
    % if nn==3
    %     I=ismember(S.filt_neurons.slice, slice_range);
    %     S.filt_neurons.brain_name(I)=repmat("AltoLH_20250317",sum(I),1);
    % end

    S.filt_neurons.uid=S.filt_neurons.brain_name+"_"+S.filt_neurons.uid;
    if nn==1
        filt_neurons=S.filt_neurons;
    else
        filt_neurons=combine_neurons_multibatch(filt_neurons,S.filt_neurons);
    end
end
%filt_neurons.expmat=sparse(full(filt_neurons.expmat));
%%
% change to a single gene codebook
% filt_neurons.genes=filt_neurons.genes{1};
% filter out cells in nan slices
to_remove=isnan(filt_neurons.slice)|filt_neurons.slice==0;
filt_neurons_nonan=subset_neurons(filt_neurons,filt_neurons.uid(~to_remove));


%% save
filt_neurons=filt_neurons_nonan;
data_folder=fullfile(pwd,"data");

mkdir(data_folder);

save(fullfile(data_folder,'filt_neurons_combined.mat'), ...
    'filt_neurons', ...
    '-v7.3');
%% write uid and brain_name to csv
writematrix(filt_neurons.uid,fullfile(data_folder,'uid.csv'));
writematrix(filt_neurons.brain_name,fullfile(data_folder,'brain_name.csv'));


