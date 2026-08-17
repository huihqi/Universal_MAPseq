clear all;
close all;

%% load all data

load colors30.mat

load Parrot_normalized_targets.mat
parrotdata=data;
parrottargets=["NLC", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM","NLC", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
clear data

load Zebra_finch_normalized_targets.mat
zebratargets=["HVC","LMAN","X","DLM","DM","nXIIts","RAM","RA","c_HVC","c_LMAN","c_X","c_DLM","c_DM","c_nXIIts","c_RAM"];
zebradata=data;
clear data


load MOp_MAPseq2.mat
mousedata=MAPseq2;
mousetargets=targets;

clear MAPseq2 data_mapseq2 targets

%% deal with target regions

%parrot remove NLC
parrotdata=parrotdata(sum(parrotdata(:,[2:8,10:16]),2)>0,[2:8,10:16]);
parrottargets=parrottargets([2:8,10:16]);

parrotdata2=[sum(parrotdata(:,[1,2]),2),parrotdata(:,3:7),sum(parrotdata(:,[8,9]),2),parrotdata(:,10:14)];
parrottargets2=["ctxI","strI","thalI",'midI','nXIItsI','PAM/RAMI','ctxC','strC','thalC','midC','nXIItsC','PAM/RAMC'];

%zebrafinch remove HVC, RA
zebradata=zebradata(sum(zebradata(:,[2:7,10:15]),2)>0,[2:7,10:15]);
zebratargets=zebratargets([2:7,10:15]);
zebratargets2=parrottargets2;

%deal with mouse
ctxI=[1,2];
ctxC=[3:5];
strI=[6:8];
strC=[9:11];
thalI=[12:15];
thalC=[16:19];
midI=[20:22];
myI=[23:24];
myC=[25:26];

%original code
%{ 
mousedata2=[sum(mousedata(:,ctxI),2),sum(mousedata(:,strI),2),sum(mousedata(:,thalI),2),...
    sum(mousedata(:,midI),2),sum(mousedata(:,myI),2),...
    sum(mousedata(:,ctxC),2),sum(mousedata(:,strC),2),sum(mousedata(:,thalC),2),...
   sum(mousedata(:,midI),2),zeros(size(mousedata,1),1)];
%}
%new matrix for bilateral mouse and parrot
mousedata2=[sum(mousedata(:,ctxI),2),sum(mousedata(:,strI),2),sum(mousedata(:,thalI),2),...
    sum(mousedata(:,midI),2),zeros(size(mousedata,1),1),sum(mousedata(:,myI),2),...
    sum(mousedata(:,ctxC),2),sum(mousedata(:,strC),2),sum(mousedata(:,thalC),2),...
   zeros(size(mousedata,1),1),zeros(size(mousedata,1),1),sum(mousedata(:,myC),2)];

%collapse my in the birds
parrotdata3=[parrotdata2(:,[1:4]),sum(parrotdata2(:,5:6),2),...
    parrotdata2(:,[7:10]),sum(parrotdata2(:,11:12),2)];


zebradata3=[zebradata(:,[1:4]),sum(zebradata(:,5:6),2),...
    zebradata(:,[7:10]),sum(zebradata(:,11:12),2)];

%targets=["ctxI","strI","thalI",'midI','myI','ctxC','strC','thalC','midC','myC'];

%targets=["ctxI","strI","thalI",'midI','nXIIts','myI','ctxC','strC','thalC','midC','nXIItsc','myC'];
targets=["ctxI","strI","thalI",'midI','nXIIts','myI'];

%max normalize everyone;
%A=parrotdata3./repmat(max(parrotdata3,[],2),1,10);
%B=zebradata3./repmat(max(zebradata3,[],2),1,10);
%C=mousedata2./repmat(max(mousedata2,[],2),1,10);


%try out different sets

%parrot vs zebra ipsi, splitting nXIIts and RAM
%MAPseq1=parrotdata2(sum(parrotdata2(:,1:6),2)>0,1:6);
%MAPseq2=zebradata(sum(zebradata(:,1:6),2)>0,1:6);

%parrot vs mouse ipsi only
%MAPseq1=parrotdata3(sum(parrotdata3(:,1:5),2)>0,1:5); 
%MAPseq2=mousedata2(sum(mousedata2(:,1:5),2)>0,1:5);

%parrot vs mouse ipsi only and seperate nXIIts
%mousedata_nXIIts_0 = [mousedata2(:,1:6)];
%MAPseq1=parrotdata2(sum(parrotdata2(:,1:6),2)>0,1:6); %seperate nXIIts
%MAPseq2=mousedata_nXIIts_0(sum(mousedata_nXIIts_0(:,1:6),2)>0,1:6);

%parrot vs mouse vs zebra finch bilateral
MAPseq1=parrotdata2;
MAPseq2=mousedata2;
MAPseq3=zebradata(sum(zebradata(:,1:6),2)>0,1:12);
%parrot vs mouse contralateral
%MAPseq1=parrotdata3(sum(parrotdata3(:,6:10),2)>0,6:10);
%MAPseq2=mousedata2(sum(mousedata2(:,6:10),2)>0,6:10);

%parrot vs mouse bilateral split nxIIts
%MAPseq1=parrotdata2;
%MAPseq2=mousedata2;




%csvwrite('parrot.csv',parrotdata3);
%csvwrite('mouse.csv',mousedata2);
%csvwrite('zebra.csv',zebradata3);



%%

% log normalization
%matching=randperm(size(MAPseq2,1), size(MAPseq1,1));% to take same number of neurons, used matching is provided
matching=randperm(size(MAPseq2,1), 1905); %take same number of neuron to parrot without forebrain 1058 for zebra finch, 1539 for mouse,  1905 for bilateral mouse
matching2=randperm(size(MAPseq3,1), 1058);
pre_inputdata=[MAPseq1;MAPseq2(matching,:);MAPseq3(matching2,:)];
randorder=randperm(size(pre_inputdata,1));% used randorder is provided
inputdata=pre_inputdata(randorder,:);
datanorm=inputdata./repmat(sum(inputdata,2),1,size(inputdata,2));
logdata=log(1+datanorm*100);
scaled_logdata=logdata/max(max(logdata));
data=scaled_logdata;

%cellular ID
ID=(1:size(inputdata,1))';
%ident
%pre_ident=[repmat(1,size(MAPseq1,1),1);repmat(2,size(MAPseq2,1),1)]; %old version
pre_ident = [repmat(1,size(MAPseq1,1),1); repmat(2,size(MAPseq2(matching,:),1),1);repmat(3,size(MAPseq3(matching2,:),1),1)];
ident=pre_ident(randorder,:);
uniqident=unique(ident);

Z = linkage(data,'ward','euclidean');

% get dendrogram leaf order without showing an extra figure
tmpFig = figure('Visible','off');
[~,~,perm] = dendrogram(Z,0);
close(tmpFig);

num_clusters = 15;

% original MATLAB cluster IDs
Traw = cluster(Z,'maxclust',num_clusters);

% reorder according to dendrogram
sorteddata = data(perm,:);
identsorted = ident(perm);
IDsorted = ID(perm);
Tsorted_raw = Traw(perm);


%% Relabel clusters according to heatmap position
% Cluster 1 = top cluster in heatmap
% Cluster 2 = second cluster
% ...
% Cluster 9 = bottom cluster

clusterOrder = unique(Tsorted_raw,'stable');

T = zeros(size(Traw));

for k = 1:num_clusters
    T(Traw == clusterOrder(k)) = k;
end

Tsorted = T(perm);


%% color settings

clusterColors = colors30_(1:num_clusters,:);

% Three species:
% 1 = parrot
% 2 = mouse
% 3 = zebra finch
identColors = colors30_([25 27 29],:);

speciesNames = ["Parrot","Mouse","Zebra finch"];


%% make species identity color bar

identColorBar = zeros(length(ident),1,3);

for i = 1:length(uniqident)

    idx = find(identsorted == uniqident(i));

    identColorBar(idx,1,:) = ...
        repmat(reshape(identColors(i,:),1,1,3), ...
        length(idx),1,1);

end


%% make cluster color bar

clusterColorBar = zeros(length(Tsorted),1,3);

for k = 1:num_clusters

    idx = find(Tsorted == k);

    clusterColorBar(idx,1,:) = ...
        repmat(reshape(clusterColors(k,:),1,1,3), ...
        length(idx),1,1);

end


%% plotting single-neuron projectome
% dendrogram | cluster bar | species bar | heatmap

figure;
set(gcf,'Units','normalized',...
    'Position',[0,0,0.45,1]);


%% Heatmap
subplot(1,100,18:100);

imagesc(sorteddata);
colormap('parula');

set(gca,...
    'XTick',1:numel(targets),...
    'XTickLabel',targets,...
    'XTickLabelRotation',45,...
    'YAxisLocation','right',...
    'YTick',[]);

ax = gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 12;

ylabel('Single neurons');


%% species color bar
subplot(1,100,14:16);

image(identColorBar);

set(gca,...
    'YTick',[],...
    'XTick',[]);

axis off;


%% cluster color bar
subplot(1,100,11:12);

image(clusterColorBar);

set(gca,...
    'YTick',[],...
    'XTick',[]);

axis off;


%% dendrogram
subplot(1,100,1:9);

H = dendrogram(Z,0,...
    'Orientation','left',...
    'Reorder',perm(end:-1:1));

set(gca,...
    'YTick',[],...
    'XTick',[]);

axis off;


%% first set all dendrogram branches to black

for j = 1:length(H)
    H(j).Color = [0 0 0];
    H(j).LineWidth = 2.0;
end


%% find row range for every cluster

H_mn = zeros(1,num_clusters);
H_mx = zeros(1,num_clusters);

for k = 1:num_clusters

    rows = find(Tsorted == k);

    H_mn(k) = min(rows);
    H_mx(k) = max(rows);

end


%% color dendrogram branches according to cluster

N = length(perm);

for k = 1:num_clusters

    for j = 1:length(H)

        branchRows = N - H(j).YData + 1;

        if min(branchRows) >= H_mn(k) && ...
           max(branchRows) <= H_mx(k)

            H(j).Color = clusterColors(k,:);

        end

    end

end


%% print cluster positions in heatmap

fprintf('\n===== Cluster positions in heatmap =====\n');

for k = 1:num_clusters

    rows = find(Tsorted == k);

    fprintf('Cluster %d: heatmap rows %d to %d (n=%d)\n',...
        k,...
        min(rows),...
        max(rows),...
        numel(rows));

end


%% cluster composition summary

cluster_summary = table(...
    'Size',[num_clusters 8],...
    'VariableTypes',...
    {'double','double','double','double','double',...
     'double','double','double'},...
    'VariableNames',...
    {'Cluster',...
     'N',...
     'N_Parrot',...
     'N_Mouse',...
     'N_Zebra',...
     'Frac_Parrot',...
     'Frac_Mouse',...
     'Frac_Zebra'});


%% mean projection profile

mean_profile = zeros(num_clusters,size(datanorm,2));


for k = 1:num_clusters

    % neurons belonging to this cluster
    idx = find(T == k);

    nParrot = sum(ident(idx) == 1);
    nMouse  = sum(ident(idx) == 2);
    nZebra  = sum(ident(idx) == 3);

    nTotal = numel(idx);


    cluster_summary.Cluster(k) = k;

    cluster_summary.N(k) = nTotal;

    cluster_summary.N_Parrot(k) = nParrot;
    cluster_summary.N_Mouse(k)  = nMouse;
    cluster_summary.N_Zebra(k)  = nZebra;

    cluster_summary.Frac_Parrot(k) = nParrot/nTotal;
    cluster_summary.Frac_Mouse(k)  = nMouse/nTotal;
    cluster_summary.Frac_Zebra(k)  = nZebra/nTotal;


    % mean projection using row-normalized data
    mean_profile(k,:) = mean(datanorm(idx,:),1);

end


disp(cluster_summary);


%% mean projection heatmap

figure;

imagesc(mean_profile);

colormap('parula');
colorbar;

set(gca,...
    'XTick',1:numel(targets),...
    'XTickLabel',targets,...
    'XTickLabelRotation',45,...
    'YTick',1:num_clusters,...
    'YTickLabel',1:num_clusters);

xlabel('Projection targets');
ylabel('Cluster');

title('Mean projection profile');

ax = gca;
ax.FontSize = 12;


cluster_mats = struct();

for k = 1:num_clusters
    idx = find(T == k);              
    
    idx1 = idx(ident(idx) == 1);     
    idx2 = idx(ident(idx) == 2);     
    idx3 = idx(ident(idx) == 3); 
    % keep index
    cluster_mats(k).idx_all = idx;
    cluster_mats(k).idx_MAPseq1 = idx1;
    cluster_mats(k).idx_MAPseq2 = idx2;
    cluster_mats(k).idx_MAPseq3 = idx3;
   
    cluster_mats(k).MAPseq1_datanorm = datanorm(idx1, :);
    cluster_mats(k).MAPseq2_datanorm = datanorm(idx2, :);
    cluster_mats(k).MAPseq3_datanorm = datanorm(idx3, :);
    
    % cluster_mats(k).MAPseq1_inputdata = inputdata(idx1, :);
    % cluster_mats(k).MAPseq2_inputdata = inputdata(idx2, :);
    
    fprintf('Cluster %d: MAPseq1=%d rows, MAPseq2=%d rows\n', k, numel(idx1), numel(idx2));
end


save('cluster_mats_parrot_mouse_08082026_three_species.mat','cluster_mats','targets');

outdir = 'cluster_exports_08082026';
if ~exist(outdir,'dir'); mkdir(outdir); end

for k = 1:num_clusters
    writematrix(cluster_mats(k).MAPseq1_datanorm, fullfile(outdir, sprintf('cluster%d_MAPseq1_datanorm.csv',k)));
    writematrix(cluster_mats(k).MAPseq2_datanorm, fullfile(outdir, sprintf('cluster%d_MAPseq2_datanorm.csv',k)));
    writematrix(cluster_mats(k).MAPseq3_datanorm, fullfile(outdir, sprintf('cluster%d_MAPseq3_datanorm.csv',k)));
end
clear
%Extract the matrix of each cluster
%C2_M1 = cluster_mats(2).MAPseq1_datanorm;
%C2_M2 = cluster_mats(2).MAPseq2_datanorm;