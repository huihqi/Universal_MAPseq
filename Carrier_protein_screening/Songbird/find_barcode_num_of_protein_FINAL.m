Final_Matrix=[];
sample_info=[];
for sample=1:12
        name=['sample',num2str(sample)];
        index=sample;
        %matrix1=['Total_Barcode','Target_Barcode','Thalamus','SNr','Medullla'];
        %matrix1=[0,0,0,0,0];
        %writematrix(matrix1, [name,'_inj_Barcode_gt200_','ALL_GENE_BARCODE_10_brief.csv'],'WriteMode', 'append');
        %writematrix(matrix1, [name,'ALL_GENE_BARCODE_10_brief.csv'],'WriteMode', 'append');

        Carrier_Protein_list = {'AP2B1', 'AP2M1', 'CACNG4', 'CADM1', 'DNM1', 'ERC2', 'GAP43', 'PTPRA', 'RAB3A', 'RAB6B', 'SNAP25', 'STX1B', 'STX6', 'SYNAPSINI', 'VAMP2_SYNAPTOPHYSIN', 'SYNAPTOPHYSIN', 'SYNGR3', 'SYT1', 'VAMP1', 'VAMP2', 'VAMP4', 'STX1B_4nl', 'STX1B_PP7', 'SYNAPTOPHYSIN_4nl', 'SYNAPTOPHYSIN_PP7', 'VAMP2_4nl', 'VAMP2_PP7', 'DNM1_4nl', 'DNM1_PP7', 'mVAMP2_HK16', 'MAPP_JK100'}; 
        for CP=1:length(Carrier_Protein_list) 
                Carrier_Protein=Carrier_Protein_list{CP};
                load([name,'_','banfilt_',Carrier_Protein,'.mat'])
                load([name,'_barcodematrixBird_CP_Screening_',name,'__',Carrier_Protein,'.mat'])
                dimensions = size(barcodematrix);
                numRows = dimensions(1);
                numCols = dimensions(2);
                
                sample_num = numCols;
                %injection_site_started = sample_num - 1;
                injection_site_sample_info=sample_num; 
                proj=sample_num-1;
                
                
                projthresh = 5;
                sourcethresh = 10;
                %sourcesite=49:61;
                %projsite=1:48;
                sourcesite=injection_site_sample_info;
                projsite=1:proj;
                B=barcodematrix(max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
                if ~isempty(B)
                    if length(bcnfilt)>1
                        
                        B=B(max(B(:,projsite),[],2)<5000,:);
                        %B=barcodematrix(max(barcodematrix(:,sourcesite),[],2)>sourcethresh & max(barcodematrix(:,sourcesite),[],2)> max(barcodematrix(:,projsite),[],2),:);
                        %Bseq=refbarcodes(max(barcodematrix(:,sourcesite),[],2)>sourcethresh & max(barcodematrix(:,sourcesite),[],2)> max(barcodematrix(:,projsite),[],2),:);
                        %B=barcodematrix(max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
                        
                        %save('filtbarcode_source_AP2M1.mat','B');
                        %save('filtbarcode_source_AP2M1_seq.mat','Bseq');
                        
                        %barcode_target= B(max(B(:,projsite),[],2)>projthresh,:);
                        
                        %save('filtbarcode_target_AP2M1.mat','B')
                        
                        %sample_num = 60;
                        %injection_site_started = sample_num - 12;
                        %injection_site_sample_info=injection_site_started:sample_num;
                        %Carrier_Protein='VAMP4';
                        %writematrix(Carrier_Protein, [name,'_inj_Barcode_gt200_','ALL_GENE_BARCODE_10.csv'],'WriteMode', 'append');
                        %writematrix(Carrier_Protein, [name,'ALL_GENE_BARCODE_10_brief.csv'],'WriteMode', 'append');
                        
                        %for j = 1:13
                        j=1;
                        k=injection_site_sample_info(j);
                        disp(k);
                        B=B(max(B(:,k),[],2)>sourcethresh,:);
                        dimensions_B = size(B);
                        numRows_B = dimensions_B(1);
                        disp(numRows_B);
                        
                        barcode_target= B(max(B(:,projsite),[],2)>projthresh,:);
                        dimensions_barcode_target = size(barcode_target);
                        numRows_barcode_target = dimensions_barcode_target(1);
                        disp(numRows_barcode_target);
                        
                        
                        %numCols = dimensions(2);
                        RA = [];
                        X = [];
                        Cerebellum = [];
                        
                        
                        final = [];        
                        %For total mouse
                        RA_ori = [1,4,7,10,13,16,19,22,25,28,31,34];
                        
                        
                        X_ori = [2,5,8,11,14,17,20,23,26,29,32,35];
                        
                        %Thalamus_ori = [3,8,11,15,19,21,25,29,33,37,41,45,49];
                        
                        Cerebellum_ori = [3,6,9,12,15,18,21,24,27,30,33,36];
        
        
                        %For total mouse
                        %Medulla_ori = [1,5,9,13,17,23,27,31,35,39,43,47,51];
                        
                        
                        %SNr_ori = [2,6,10,14,18,22,26,30,34,38,42,46,50];
                        
                        %Thalamus_ori = [3,8,11,15,19,21,25,29,33,37,41,45,49];
                        
                        %Cerebellum_ori = [4,7,12,16,20,24,28,32,36,40,44,48,52];
                        
                        for i = 1:length(bcnfilt)  % Corrected loop range
                             region = bcnfilt(i);  % Get the current sample
                            
                             if ismember(region, RA_ori(index))
                                 RA = [RA, i];  % Append index to Medulla
                             elseif ismember(region, X_ori(index))
                                 X = [X, i];          % Append index to SNr
                             %elseif ismember(region, Thalamus_ori(index))
                                 %Thalamus = [Thalamus, i]; % Append index to Thalamus
                             elseif ismember(region, Cerebellum_ori(index))
                                 Cerebellum = [Cerebellum, i]; % Append index to Cerebellum
                             end
                         end 
                         
                         if ~isempty(RA)
                            RA_list = B(:,RA);
                         else
                            zeroColumn = zeros(size(B, 1), size(B, 1));
                            RA_list = zeroColumn(:,1);
                         end
                         
                         if ~isempty(X)
                            X_list = B(:,X);
                         else
                            zeroColumn = zeros(size(B, 1), size(B, 1));
                            X_list = zeroColumn(:,1);
                         end
                         
                         
                    else
                        zeroColumn = zeros(size(B, 1), size(B, 1));
                        RA_list = zeroColumn(:,1);
                        %Thalamus_list = zeroColumn(:,1);
                        X_list = zeroColumn(:,1);
                    end
                     All_Barcode=B(:,end);
                     Carrier_Protein_Column = repmat({Carrier_Protein}, size(B, 1), 1);
                     Carrier_Protein_Column = cellstr(Carrier_Protein_Column);
                     Sample_Column=repmat({name}, size(B, 1), 1);
                     Sample_Column = cellstr(Sample_Column);
                     final=[RA_list,X_list,All_Barcode];
                     information=[Carrier_Protein_Column,Sample_Column];
                     sample_info=[sample_info;information];
                     Final_Matrix=[Final_Matrix;final];
                 
                end
                
        
        end
end
% Save the matrix to a .mat file
save('Final_Matrix_5000.mat', 'Final_Matrix','sample_info');
sample_info_table = cell2table(sample_info);
Final_Matrix_table = array2table(Final_Matrix);
final_table = [Final_Matrix_table, sample_info_table];
writetable(final_table, 'Final_matrix_With_Protein_sample_info_5000.csv');
