%filter barcodes by thresholds
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end

sourcesite=[1];
projsite=[2:17];
projthresh=5;
sourcethresh=75;
%filter barcodes
B1=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);
Bseq1=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh & barcodematrix(:,sourcesite)> max(barcodematrix(:,projsite),[],2),:);

B2=B1(B1(:,2) ./ B1(:,1)<0.1,:);
Bseq2=Bseq1(B1(:,2) ./ B1(:,1)<0.1,:);

projsite1=[3:17];
B=B2(max(B2(:,projsite1),[],2) ./ B2(:,1) <0.1,:);
Bseq=Bseq2(max(B2(:,projsite1),[],2) ./ B2(:,1) <0.1,:);
B(B==1)=0;
B(B==2)=0;
%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);

rows=all(Bnorm(:, 2)==0,2);
Bnorm1=Bnorm(rows,:);
Bseq3=Bseq(rows,:);l
B3=B(rows,:);
%B_filt=B(max(B,[],2)<cellbodyThresh,:);
%Bnorm_filt=Bnorm(max(B,[],2)<cellbodyThresh,:);l
%Bseq_filt=Bseq(max(B,[],2)<cellbodyThresh,:);

%targets only
%B_filt_tar=B_filt(:,projsite);
%Bnorm_filt_tar=Bnorm_filt(:,projsite);
B_tar=B3(:,projsite1);
Bnorm2_tar=Bnorm1(:,projsite1);
%save('filtBCmat_75.mat','B','Bnorm','B_filt','Bseq','Bnorm_filt','Bseq_filt','B_filt_tar','Bnorm_filt_tar');
save('filtBCmat_bird_GP_R.mat','B3','Bnorm1','B_tar','Bseq3','Bnorm2_tar','-v7.3');

