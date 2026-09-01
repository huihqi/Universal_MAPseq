clearvars
clc

load filt_neurons-bc-somas.mat
output_dir=fullfile(pwd,'figs','cell_crops');
mkdir(output_dir)

data_dir=fullfile('\\allen\programs\celltypes\workgroups\barcoded_connectomics\20250121_MacaqueBARseq_alder_BB_frontal_7_maxprojection\','processed');
%% 

%% Translate filt_neurons dom_bc code to MAPseq code
nucleotides='GTAC';
dom_bc=char(filt_neurons.soma_bc);
dom_bc(:)='N';
for nn=1:numel(filt_neurons.dom_bc)
    if ~isempty(filt_neurons.dom_bc{nn})
        dom_bc(nn,:)=nucleotides(filt_neurons.dom_bc{nn});
    end
end
%%



%% Match to viral codebook (no mismatch)
load viral_codebook2_corr;

codebook(:,2)=upper(codebook(:,2));

codes=char(codebook(:,2));

[~,viral_idx]=ismember(dom_bc,codes(:,1:size(dom_bc,2)),'rows');

sum(viral_idx>0)
filt_neurons.viral_idx=viral_idx;
filt_neurons.viral_codebook=codebook;



count_thresh=10;
complexity_thresh=0.2;
[~,c]=bccomplexity(filt_neurons.dom_bc);
is_bc=filt_neurons.dom_bc_count>=count_thresh&c>complexity_thresh; % unlike xc01, is_bc here combines both count threshold and complexity threshold
filt_neurons.is_bc=is_bc;



save('filt_neurons-viralidx.mat','filt_neurons','-v7.3');
