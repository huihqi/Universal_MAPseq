function rankplot(bcn,prefix,donotplot)
%plot rank plot of counts
%% read in sequences and counts of data split into the different libraries

bcn=1:65;
prefix='Mouse_CP_Screening';
for i=1:length(bcn)
data(i).rawcounts=dlmread([prefix,'_BC',int2str(bcn(i)),'_BU_counts.txt']);
end
save('rawdata.mat','data','-v7.3');

%use this to look at the sequence rank plot of every libary and choose a threshold for preprocessing.sh 
load rawdata
if ~exist('donotplot')
    figure;for i=1:length(bcn);semilogy(data(i).rawcounts);title(int2str(i));xlabel('Sequence rank');
        ylabel('Read count');
        pause();end

        saveas([fig,'RankPlot_Fig/',int2str(i),'.tiff'])
elseif donotplot~=1
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(i));
        xlabel('Sequence rank');
        ylabel('Read count');
        pause();end

        saveas([fig,'RankPlot_Fig/',int2str(i),'.tiff'])
end