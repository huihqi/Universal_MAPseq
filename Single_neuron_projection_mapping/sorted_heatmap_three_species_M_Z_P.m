
S = combined_cluster_mats;   % combined_cluster_mats

num_clusters = numel(S);


nT = [];
for k = 1:num_clusters
    if ~isempty(S(k).MAPseq1_datanorm), nT = size(S(k).MAPseq1_datanorm,2); break; end
    if ~isempty(S(k).MAPseq2_datanorm), nT = size(S(k).MAPseq2_datanorm,2); break; end
    if ~isempty(S(k).MAPseq3_datanorm), nT = size(S(k).MAPseq3_datanorm,2); break; end
end
if isempty(nT), error('All clusters are empty.'); end


nRows_each = arrayfun(@(x) numel(x.idx_all), S(:));
totalRows = sum(nRows_each);
endRow = cumsum(nRows_each);


Xall   = nan(totalRows, nT);
srcAll = zeros(totalRows, 1);   
% 1 = MAPseq1 (Parrot)
% 2 = MAPseq2 (Mouse)
% 3 = MAPseq3 (Zebra finch)
% 0 = unknown

rowPtr = 1;

for k = 1:num_clusters
    idx_all = S(k).idx_all(:);
    idx1    = S(k).idx_MAPseq1(:);
    idx2    = S(k).idx_MAPseq2(:);
    idx3    = S(k).idx_MAPseq3(:);
    A1      = S(k).MAPseq1_datanorm;
    A2      = S(k).MAPseq2_datanorm;
    A3      = S(k).MAPseq3_datanorm;

    if isempty(idx_all)
        continue
    end


    map1 = containers.Map('KeyType','double','ValueType','double');
    for i = 1:numel(idx1), map1(idx1(i)) = i; end

    map2 = containers.Map('KeyType','double','ValueType','double');
    for i = 1:numel(idx2), map2(idx2(i)) = i; end

    map3 = containers.Map('KeyType','double','ValueType','double');
    for i = 1:numel(idx3), map3(idx3(i)) = i; end


    Xmix = nan(numel(idx_all), nT);
    src  = zeros(numel(idx_all), 1);

    for r = 1:numel(idx_all)
        g = idx_all(r);
        if ~isempty(A1) && isKey(map1, g)
            Xmix(r,:) = A1(map1(g), :);
            src(r) = 1;
        elseif ~isempty(A2) && isKey(map2, g)
            Xmix(r,:) = A2(map2(g), :);
            src(r) = 2;
        elseif ~isempty(A3) && isKey(map3, g)
            Xmix(r,:) = A3(map3(g), :);
            src(r) = 3;
        else
            src(r) = 0;
        end
    end


    Xclu = Xmix;
    Xclu(isnan(Xclu)) = 0;

    if size(Xclu,1) >= 2
 
        D = pdist(Xclu, 'euclidean');

        Z = linkage(D, 'ward');


        [~,~,perm] = dendrogram(Z, 0);


        Xmix = Xmix(perm, :);
        src  = src(perm);
    end


    rr = rowPtr:(rowPtr + numel(idx_all) - 1);
    Xall(rr,:)   = Xmix;
    srcAll(rr,:) = src;

    rowPtr = rowPtr + numel(idx_all);
end


rowBar = srcAll;
rowBar(rowBar==0) = NaN;


fig = figure('Color','w','Units','pixels','Position',[50 50 2600 900]);

leftMargin   = 0.05;
rightMargin  = 0.02;
bottomMargin = 0.08;
topMargin    = 0.06;

gap  = 0.01;
barW = 0.010;          ）
mainW = 1 - leftMargin - rightMargin - barW - gap;
H = 1 - bottomMargin - topMargin;


% 1 = MAPseq1 (Parrot)
% 2 = MAPseq2 (Mouse)
% 3 = MAPseq3 (Zebra finch)

axL = axes('Parent',fig,...
    'Position',[leftMargin, bottomMargin, barW, H]);

hBar = imagesc(axL, rowBar);

axis(axL,'tight');
axis(axL,'off');



% speciesColors = colors30_([25 27 29],:);

% or
speciesColors = [

    0.35  0.35  0.65;    % blue-purple

    0.10  0.65  0.58;    % teal

    0.90  0.45  0.15     % orange

];

colormap(axL, speciesColors);

% categorical values 1, 2, 3
caxis(axL,[0.5 3.5]);

% unknown = NaN -> transparent/white
set(hBar,'AlphaData',~isnan(rowBar));
set(axL,'Color','w');


axH = axes('Parent',fig,'Position',[leftMargin+barW+gap, bottomMargin, mainW, H]);
imagesc(axH, Xall);
axis(axH,'tight');
colormap(axH, parula);
caxis(axH,[0 1]);
cb = colorbar(axH);
cb.Label.String = 'datanorm';

set(axH, 'XTick', [], 'YTick', [], 'TickLength', [0 0]);


% box(axH, 'off');


box(axH, 'on');  

title(axH, 'All clusters in order; within each cluster: Euclidean clustering (Ward)');


linkaxes([axL, axH],'y');

%{

hold(axH,'on');
for k = 1:num_clusters-1
    yline(axH, endRow(k)+0.5, '-', 'Color', [0 0 0], 'LineWidth', 0.75);
end
hold(axH,'off');
%}


tickLen = 0.035;        
x0 = 0.5;                

hold(axH,'on');
for k = 1:num_clusters-1
    yb = endRow(k) + 0.5; 
    plot(axH, [x0, x0+tickLen], [yb, yb], 'k-', 'LineWidth', 1.0);
end
hold(axH,'off');


% set(axH,'XTick',1:nT,'XTickLabel',targets(1:nT),'XTickLabelRotation',45);