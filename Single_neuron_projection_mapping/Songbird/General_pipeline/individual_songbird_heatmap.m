
% draw each birds matrix
% bird O_L

newCols = [4 1 2 7 3 6 5 9 14 11 12 17 13 16 15];

O_L_new = Bnorm2_tar(:, newCols);

targets=1:15;
data_tar=O_L_new(:,1:15);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;        % Element-wise division

clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')
% DM percentage out of decending projection

cols = 4:7; %DLM DM nXIIts,PAM/RAM
nRows_DLM_DM_nXIIts_PAM = sum( any(O_L_new(:, cols) > 0, 2) );
nRows_DM = sum( any(O_L_new(:, 5) > 0, 2) );
percentage=nRows_DM / nRows_DLM_DM_nXIIts_PAM; % 32.9%

nRows_nXIIts = sum( any(O_L_new(:, 6) > 0, 2) );
percentage_nXIIts=nRows_nXIIts / nRows_DLM_DM_nXIIts_PAM; % 58.4%

% bird_B_R

newCols = [4 1 2 7 3 6 5 9 14 11 12 17 13 16 15];

B_R_new = Bnorm2_tar(:, newCols);

targets=1:15;
data_tar=B_R_new(:,1:15);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;        % Element-wise division

clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')
% DM percentage out of decending projection

cols = 4:7; %DLM DM nXIIts,PAM/RAM
nRows_DLM_DM_nXIIts_PAM = sum( any(B_R_new(:, cols) > 0, 2) );
nRows_DM = sum( any(B_R_new(:, 5) > 0, 2) );
percentage=nRows_DM / nRows_DLM_DM_nXIIts_PAM; % 48.1%

nRows_nXIIts = sum( any(B_R_new(:, 6) > 0, 2) );
percentage_nXIIts=nRows_nXIIts / nRows_DLM_DM_nXIIts_PAM; % 54.7%

% bird_GP_R

newCols = [4 1 2 7 3 6 5 9 14 11 12 17 13 16 15];

GP_R_new = Bnorm2_tar(:, newCols);

targets=1:15;
data_tar=GP_R_new(:,1:15);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;        % Element-wise division

clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')
