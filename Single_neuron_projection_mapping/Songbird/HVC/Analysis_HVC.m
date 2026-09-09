%filter barcodes by thresholds
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end

sourcesite=[1];
projsite=[2:5];
projthresh=5;
sourcethresh=75;
%filter barcodes
B1=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);
Bseq1=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);

B2=B1(B1(:,2) ./ B1(:,1)<0.2,:);
Bseq2=Bseq1(B1(:,2) ./ B1(:,1)<0.2,:);

projsite1=[2:5];
B=B2(max(B2(:,projsite1),[],2) ./ B2(:,1) <0.2,:);
Bseq=Bseq2(max(B2(:,projsite1),[],2) ./ B2(:,1) <0.2,:);
B(B==1)=0;
B(B==2)=0;
%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);


%B_filt=B(max(B,[],2)<cellbodyThresh,:);
%Bnorm_filt=Bnorm(max(B,[],2)<cellbodyThresh,:);l
%Bseq_filt=Bseq(max(B,[],2)<cellbodyThresh,:);

%targets only
%B_filt_tar=B_filt(:,projsite);
%Bnorm_filt_tar=Bnorm_filt(:,projsite);
B_tar=B(:,projsite1);
Bnorm2_tar=Bnorm(:,projsite1);
%save('filtBCmat_75.mat','B','Bnorm','B_filt','Bseq','Bnorm_filt','Bseq_filt','B_filt_tar','Bnorm_filt_tar');
save('filtBCmat_bird_66_R.mat','Bnorm','B_tar','Bseq','Bnorm2_tar','B','-v7.3');

data1=Bnorm2_tar;
Bnorm_tar1=data1(:,[1:2 4]);

Bnorm_tar1 = Bnorm_tar1(any(Bnorm_tar1, 2), :);
%targets=1:4;
% targets=["RA", "X", "Av", "CB"];
targets=["RA", "X", "CB"];
data_tar=Bnorm_tar1(:,1:3);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;        % Element-wise division

clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')

data1=data(:, 1:7);

data1=Bnorm_tar1;
thr = 0;              % change if needed (e.g., thr = 1)
B = data1 > thr;          % binarize: 1 = projects to that area

k = sum(B, 2);        % how many areas each neuron projects to

n1 = sum(k == 1);
n2 = sum(k == 2);
n3 = sum(k == 3);
n4 = sum(k == 4);
n5 = sum(k == 5);
n6 = sum(k == 6);
n7 = sum(k == 7);

fprintf('Exactly 1 area: %d\n', n1);
fprintf('Exactly 2 areas: %d\n', n2);
fprintf('Exactly 3 areas: %d\n', n3);
fprintf('Exactly 4 areas: %d\n', n4);
fprintf('Exactly 5 areas: %d\n', n5);
fprintf('Exactly 6 areas: %d\n', n6);
fprintf('Exactly 7 areas: %d\n', n7);
