% Parrot 5
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end
Bnorm=barcodematrix./repmat(x,size(barcodematrix,1),1);

barcodematrix1 = barcodematrix(:, [13:24, 1:12, 25:end]);
Bnorm1=Bnorm(:, [13:24, 1:12, 25:end]);

filtered_data = barcodematrix1(any(barcodematrix1(:, 5:22) > 0, 2), :);
filtered_data_Bnorm=Bnorm1(any(barcodematrix1(:, 5:22) > 0, 2), :);

filtered_data_75 = filtered_data(filtered_data(:, 1) > 75, :);
filtered_data_75_Bnorm=filtered_data_Bnorm(filtered_data(:, 1) > 75, :);

denominator = max(filtered_data_75(:, [1, 2, 4]), [], 2);
idx1 = filtered_data_75(:, 3) ./ denominator < 0.1;
filtered_data_non_10 = filtered_data_75(idx1, :);
filtered_data_non_10_Bnorm = filtered_data_75_Bnorm(idx1, :);
%filtered_data_non_10 = filtered_data_75(filtered_data_75(:, 3) ./ denominator < 0.1, :);
%filtered_data_non_10(filtered_data_non_10 <= 2) = 0;

filtered_data_non_10(filtered_data_non_10 <= 2) = 0;
filtered_data_non_10_Bnorm(filtered_data_non_10 <= 2) = 0;
idx2 = any(filtered_data_non_10(:, 5:25) ~= 0, 2);

filtered_data_non_10_1 = filtered_data_non_10(idx2, :);
filtered_data_non_10_1_Bnorm = filtered_data_non_10_Bnorm(idx2, :);

%filtered_data_non_10_1 = filtered_data_non_10(any(filtered_data_non_10(:, 5:25) ~= 0, 2), :);
max_1_4  = max(filtered_data_non_10_1(:, 1:4), [], 2);
max_5_25 = max(filtered_data_non_10_1(:, 5:25), [], 2);
idx2 = max_1_4 >= 5 .* max_5_25;
filtered_data_final = filtered_data_non_10_1(idx2, :);
filtered_data_final_Bnorm = filtered_data_non_10_1_Bnorm(idx2, :);

%filtered_data_final_binary=filtered_data_final;
%filtered_data_final_binary(:, 5:25) = filtered_data_final_binary(:, 5:25) > 0;
idx4 = any(filtered_data_final_Bnorm(:, [6:12 18:24]) > 0, 2);
data_tar = filtered_data_final_Bnorm(idx4, [6:12 18:24]);
%data_tar=filtered_data_final_Bnorm(:,[6:12 18:24]);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;  
targets=[ "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
'Colormap','parula')
save('data_5_norm.mat','data','-v7.3');

% Parrot 4
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end
Bnorm=barcodematrix./repmat(x,size(barcodematrix,1),1);

load('barcodematrixParrot_AACd_parrot_4_07162026.mat')
filtered_data = barcodematrix(any(barcodematrix(:, 5:22) > 0, 2), :);
filtered_data_Bnorm=Bnorm(any(barcodematrix(:, 5:22) > 0, 2), :);

filtered_data_75 = filtered_data(filtered_data(:, 1) > 75, :);
filtered_data_75_Bnorm=filtered_data_Bnorm(filtered_data(:, 1) > 75, :);

denominator = max(filtered_data_75(:, [1, 2, 4]), [], 2);
idx1 = filtered_data_75(:, 3) ./ denominator < 0.1;
filtered_data_non_10 = filtered_data_75(idx1, :);
filtered_data_non_10_Bnorm = filtered_data_75_Bnorm(idx1, :);
%filtered_data_non_10 = filtered_data_75(filtered_data_75(:, 3) ./ denominator < 0.1, :);
%filtered_data_non_10(filtered_data_non_10 <= 2) = 0;

filtered_data_non_10(filtered_data_non_10 <= 2) = 0;
filtered_data_non_10_Bnorm(filtered_data_non_10 <= 2) = 0;
idx2 = any(filtered_data_non_10(:, 5:25) ~= 0, 2);

filtered_data_non_10_1 = filtered_data_non_10(idx2, :);
filtered_data_non_10_1_Bnorm = filtered_data_non_10_Bnorm(idx2, :);

%filtered_data_non_10_1 = filtered_data_non_10(any(filtered_data_non_10(:, 5:25) ~= 0, 2), :);
max_1_4  = max(filtered_data_non_10_1(:, 1:4), [], 2);
max_5_25 = max(filtered_data_non_10_1(:, 5:25), [], 2);
idx2 = max_1_4 >= 5 .* max_5_25;
filtered_data_final = filtered_data_non_10_1(idx2, :);
filtered_data_final_Bnorm = filtered_data_non_10_1_Bnorm(idx2, :);

%filtered_data_final_binary=filtered_data_final;
%filtered_data_final_binary(:, 5:25) = filtered_data_final_binary(:, 5:25) > 0;
idx4 = any(filtered_data_final_Bnorm(:, [6:12 18:24]) > 0, 2);
data_tar = filtered_data_final_Bnorm(idx4, [6:12 18:24]);
%data_tar=filtered_data_final_Bnorm(:,[6:12 18:24]);
row_sums = sum(data_tar, 2);                     % Sum each row (along columns)
data = data_tar ./ row_sums;  
targets=[ "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM", "NAO", "MO", "LPO", "DMm", "DM", "nXIIts", "PAM/RAM"];
clustergram(data,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
'Colormap','parula')

save('data_4_norm.mat','data','-v7.3');
