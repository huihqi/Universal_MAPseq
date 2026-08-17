function rankplot(bcn,prefix,donotplot)
%plot rank plot of counts
%% read in sequences and counts of data split into the different libraries

bcn=1:60;
prefix='Songbird_RA';
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
        saveas([gcf,'RankPlot_Fig/',int2str(i),'.tiff']);
        
elseif donotplot~=1
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(i));
        xlabel('Sequence rank');
        ylabel('Read count');
        
        pause();end
        saveas([gcf,'RankPlot_Fig/',int2str(i),'.tiff']);
        
end




Macaque_CP_Screening_V2_BU_counts.txt



% store the image

mkdir('RankPlot_Fig'); % Ensure the folder exists

for i = 1:length(bcn)
    % Create a new figure without displaying it
    fig = figure('Visible', 'off');
    
    % Plot the data using semilogy
    semilogy(data(i).rawcounts);
    title(int2str(bcn(i)));
    xlabel('Sequence rank');
    ylabel('Read count');
    
    % Generate the filename for saving
    filename = fullfile('RankPlot_Fig', [int2str(bcn(i)), '.tiff']);
    
    % Save the figure using print (faster than saveas)
    print(fig, filename, '-dtiff');
    
    % Close the figure to free memory
    close(fig);
end