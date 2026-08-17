
function conditionalP4(matrix,targetlabels,cutoff,min_region,min_cts,pscale,dicescale,rarityscale,dscale,dolabel,doN_t)

% This function computes the following projection metrics between brain regions:
% - Conditional probability P(B|A)
% - Dice score (observed)
% - Rarity score
% - Dice score (estimated)
% - Observed/Estimated ratio with over/under-representation plotting
% data,targets,0,0,5,1,0.5,60,[1/16 1*16],1,1
% INPUTS:
% matrix      : Binary or continuous input data matrix (e.g., projection data)
% cutoff      : Threshold to binarize the matrix (e.g., 0 means any non-zero is treated as projection)
% min_region  : Minimum number of neurons projecting to a region to include it (e.g., 10)
% min_cts     : Minimum count of observed or estimated co-projections.
%               If both are below this value (e.g., 5), the region pair is marked as 'Few' in the plot.
% pscale      : Scaling range for conditional probability P(B|A) (e.g., [0 1])
% dicescale   : Scaling range for both observed and estimated Dice scores (e.g., [0 1])
% rarityscale : Scaling range for Rarity score (e.g., [1 10]); values range from 1 to ∞
% dscale      : Log-scale range for Observed/Estimated ratio plot (e.g., [1/10 10]); 
%               a value of 1 indicates non-random expectation.
% dolabel     : If set to 1, significance labels will be added to the Observed/Estimated plot.
% doN_t       : (Optional) If provided, calculates total N values following 
%               Han & Kebschull et al., 2018 (Nature). 
%               Note: Not recommended when the number of regions is very large due to computational load.

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);

matrix1=matrix.*matrix_binary_tmp;
matrix2=matrix1(:,sum(matrix_binary_tmp,1)>=min_region);
matrix3=matrix2(sum(matrix_binary_tmp2,2)~=0,:);

Plabels=targetlabels(sum(matrix_binary_tmp,1)>=min_region);


%find N_t counts
if exist('doN_t','var')
    [N_obs, k] = size(matrix_binary);
    N_a = sum(matrix_binary, 1);  % count 1s in each dimension

% Initialize polynomial coefficients array
    coeffs = zeros(1, k+1);
   
% First coefficient (highest degree)
    coeffs(1) = N_obs - sum(N_a);
   
% Calculate other coefficients using nchoosek
    for i = 2:k
   % (-1)^(i-1) gives alternating signs
        coeffs(i) = -(-1)^(i-1) * sum(prod(nchoosek(N_a,i),2));
    end
   
% Last coefficient (constant term)
    coeffs(end) = -(-1)^k * prod(N_a);
   
% Find roots
    r = roots(coeffs);
   
% Take the largest positive real root as N_total
    N_t = round(max(r));
else
    N_t=size(matrix_binary,1);
end


coP=zeros(length(Plabels));
coO=zeros(length(Plabels));
comb=nchoosek(1:length(Plabels),2);
P=sum(matrix_binary,1)/N_t;
for i=1:size(comb,1)
    r=comb(i,1);c=comb(i,2);
    coP(r,c)=P(r)*P(c);
    coP(c,r)=P(r)*P(c);
    coO(r,c)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
    coO(c,r)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
end
coE=coP*N_t;



% Calculate the original p-value (two-tailed) for the observed outcome
[X, Y] = size(coO);
pvalueP = nan(X, Y);  % 결과 저장할 행렬

for r = 1:X
    for c = 1:Y
        k_obs = coO(r, c);        % 관측된 co-occurrence 수
        p_null = coP(r, c);       % null 모델 하의 발생 확률
        n = N_t;                  % 전체 뉴런 수 (trial 횟수)

        % Binomial 분포 계산
        x = 0:n;
        prob = binopdf(x, n, p_null);

        % 관측값의 확률
        p_k = binopdf(k_obs, n, p_null);

        % Two-tailed p-value: 관측값보다 "덜 가능성 있는" 값들의 확률 합
        pvalueP(r, c) = sum(prob(prob <= p_k));
    end
end



% Apply Bonferroni correction
p_adjP = pvalueP * nchoosek(length(Plabels),2);
signifP=strings(size(p_adjP));
signifP(p_adjP>=0.05)="NS";
signifP(p_adjP<0.05 & p_adjP>=0.01)="*";
signifP(p_adjP<0.01 & p_adjP>=0.001)="**";
signifP(p_adjP<0.001 & p_adjP>=0.0001)="***";
signifP(p_adjP<0.0001)="****";
signifP(p_adjP==1)="";


for i=1:length(Plabels)
    for j=1:length(Plabels)
        if max(coE(j,i), coO(j,i))<min_cts
            signifP(j,i)="Few";
        end
    end
    signifP(i,i)="";
end


[p_adjP_idxr, p_adjP_idxc]=find((p_adjP<0.05)==1);
p_adjP_idx=[p_adjP_idxr p_adjP_idxc];


j=1;
conditionalp=zeros(length(Plabels));
numerator=zeros(length(Plabels));
for i=1:length(Plabels)
    conditionalp(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1)./sum(matrix_binary(:,i)==1);
    numerator(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1);
    j=j+1;
end

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(conditionalp');clim([0 pscale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(conditionalp')
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

conditionalp_dice=numerator;
for i=1:length(Plabels)
    for j=1:length(Plabels)
        conditionalp_dice(i,j)=2*conditionalp_dice(i,j)/(sum(matrix_binary(:,i)==1)+sum(matrix_binary(:,j)==1));
    end
end
    
pre_sizefactor=zeros(size(numerator));
for i=1:length(Plabels)
    for j=1:length(Plabels)
        pre_sizefactor(i,j)=((N_t/sum(matrix_binary(:,i)==1))+(N_t/sum(matrix_binary(:,j)==1)))/2;
    end
end

sizefactor=1./pre_sizefactor;

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(conditionalp_dice);clim([0 dicescale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(conditionalp_dice)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

rarity=pre_sizefactor;
figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(rarity);clim([1 rarityscale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(sizefactor)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;




% let's divide P(B|A) with P(B) to check over or under representation
divbyPb=N_t*conditionalp./sum(matrix_binary);

figure;imagesc(divbyPb);clim(dscale);colormap(redblue);
set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
set(gca,'ColorScale','log');
ax = gca;
ax.YTick=[];
ax.XTick=[];
[rows, cols] = size(divbyPb);
for i = 1:rows
   for j = 1:cols
        if dolabel==1
        text(j, i, num2str(divbyPb(i,j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',10);
        text(j, i-0.3, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        else
        text(j, i, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',13);
        end
    end
end
hold on;
for i = 1:length(divbyPb)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;


findfigs;

save('contionalP.mat','conditionalp','conditionalp_dice','sizefactor','matrix_binary','N_t','signifP','coE','coO','p_adjP','Plabels','numerator');