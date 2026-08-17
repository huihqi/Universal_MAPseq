%filter barcodes by thresholds
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end

sourcesite=[1];
projsite=[2:22];
projthresh=5;
sourcethresh=100;
%filter barcodes
B1=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);
Bseq1=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2) & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);

%B2=B1(B1(:,2) ./ B1(:,1)<0.1,:);
%Bseq2=Bseq1(B1(:,2) ./ B1(:,1)<0.1,:);

%projsite1=[1 3:11 14:22];
B=B1(max(B1(:,projsite),[],2) ./ B1(:,1) <0.2,:);
Bseq=Bseq1(max(B1(:,projsite),[],2) ./ B1(:,1) <0.2,:);

B(B==1)=0;
B(B==2)=0;
%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);


rows=all(Bnorm(:, [2 3 12 13 14])==0,2);
Bnorm1=Bnorm(rows,:);
Bseq2=Bseq(rows,:);
B2=B(rows,:);
%projsite1=[4:11 15:22];
%B=Bnorm1(max(Bnorm1(:,projsite),[],2) ./ Bnorm1(:,12) <0.1,:);
%Bseq=Bseq1(max(Bnorm1(:,projsite),[],2) ./ Bnorm1(:,12) <0.1,:);

%B_filt=B(max(B,[],2)<cellbodyThresh,:);
%Bnorm_filt=Bnorm(max(B,[],2)<cellbodyThresh,:);
%Bseq_filt=Bseq(max(B,[],2)<cellbodyThresh,:);

%targets only
%B_filt_tar=B_filt(:,projsite);
%Bnorm_filt_tar=Bnorm_filt(:,projsite);
B_tar=B2(:,projsite);
Bnorm_tar=Bnorm1(:,projsite);
%save('filtBCmat_75.mat','B','Bnorm','B_filt','Bseq','Bnorm_filt','Bseq_filt','B_filt_tar','Bnorm_filt_tar');
%save('Parrot2_R_5times.mat','B2','Bnorm1','B_tar','Bseq2','Bnorm_tar','-v7.3');

Bnorm_tar1=Bnorm_tar(:,[4:10 15:21]);
Bnorm_tar1 = Bnorm_tar1(any(Bnorm_tar1, 2), :);
%targets=1:16;
targets=["NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
% targets=["NLC", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
data_tar=Bnorm_tar1(:,1:14);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;        % Element-wise division

clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')

