load MOp_MAPseq2.mat
mousedata=MAPseq2;
mousetargets=targets;

clear MAPseq2 data_mapseq2 targets

load Parrot_normalized_targets.mat
parrotdata=data;
parrottargets=["NLC", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM","NLC", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
clear data

%parrot remove NLC
%parrotdata=parrotdata(sum(parrotdata(:,[2:8,10:16]),2)>0,[2:8,10:16]); %bilateral
%parrottargets=parrottargets([2:8,10:16]);
%parrotdata=parrotdata(sum(parrotdata(:,[2:8]),2)>0,[2:8]); %ipsilateral

MOI=[1];
SSI=[2];
MOC=[3:4];
SSC=[5];
STRI=[6:8];
STRC=[9:11];
THI=[12:15];
THC=[16:19];
MIDI=[20:22];
MYI=[23:24];
MYC=[25:26];


mousedata_PBA = [sum(mousedata(:,MOI),2), ...
    sum(mousedata(:,SSI),2), ...
    sum(mousedata(:,STRI),2), ...
    sum(mousedata(:,THI),2), ...
    sum(mousedata(:,MIDI),2), ...
    zeros(size(mousedata,1),1), ...
    sum(mousedata(:,MYI),2), ...
    sum(mousedata(:,MOC),2), ...
    sum(mousedata(:,SSC),2), ...
    sum(mousedata(:,STRC),2), ...
    sum(mousedata(:,THC),2), ...
    zeros(size(mousedata,1),1), ...
    zeros(size(mousedata,1),1), ...
    sum(mousedata(:,MYC),2)];

zebradata=mousedata_PBA;
clear data
%zebradata=zebradata(sum(zebradata(:,1:12),2)>0,1:12);
%idx = (zebradata(:,1) == 0) & (zebradata(:,2) == 0);  
idx = (zebradata(:,6) == 0) & (zebradata(:,13) == 0); %for mouse
X = zebradata(idx, :);                      
Xb=X>0; %for zebra finch

% Xb = mousedata_PBA > 0;                    % n×k logical for mouse
%Xb=parrotdata > 0;  %for parrot
nA = sum(Xb, 1);               % 1×k
co = double(Xb') * double(Xb); % k×k: #(Ai & Bj)

P = co ./ nA';                 
P(nA==0, :) = NaN;             
%P(nA==0, :) = 0;

% ----- 4) draw heatmap
figure('Color','w');
imagesc(P);
axis tight;


% rgb = readmatrix('/Users/huihqi/Downloads/Vocalization/parula_256.csv'); % 256×3
% colormap(rgb);


colormap(parula);


set(gca, 'Color', [1 1 1]);
h = findobj(gca, 'Type', 'Image');
set(h, 'AlphaData', ~isnan(P));

caxis([0 1]);
cb = colorbar;
cb.Label.String = "P(B|A)";

names = {'MOI','SSI','STRI','THI','MIDI','Amb','MYI','MOC','SSC','STRC','THC','MIDC','Amb','MYC'};
%names = {"NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"};
%names = {"NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"};
%names={"LMAN","X","DLM","DM","nXIIts","RAM"};
%names={"DLM","DM","nXIIts","RAM"};


set(gca, 'XTick', 1:numel(names), 'XTickLabel', names, 'XTickLabelRotation', 90);
set(gca, 'YTick', 1:numel(names), 'YTickLabel', names);
xlabel("B (columns)");
ylabel("A (rows)");

% square

fmt = '%.2f';

[k, k2] = size(P);
assert(k == k2, 'P should be square (k×k).');
%{
for i = 1:k
    for j = 1:k
        if ~isnan(P(i,j))
      
            if P(i,j) > 0.5
                txtColor = 'w';
            else
                txtColor = 'k';
            end

            text(j, i, sprintf(fmt, P(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 9, ...
                'Color', txtColor);
        end
    end
end
%}

title("P(B|A) heatmap");