
load ../rawdata
prefix='Mouse_CP_Screening'
brain='sample1_'
%gene_list={'AP2B1', 'AP2M1', 'CACNG4', 'CADM1', 'DNM1', 'ERC2', 'GAP43', 'PTPRA', 'RAB3A', 'RAB6B', 'SNAP25', 'STX1B', 'STX6', 'SYNAPSINI', 'VAMP2_SYNAPTOPHYSIN', 'SYNAPTOPHYSIN', 'SYNGR3', 'SYT1', 'VAMP1', 'VAMP2', 'VAMP4', 'STX1B_4nl', 'STX1B_PP7', 'SYNAPTOPHYSIN_4nl', 'SYNAPTOPHYSIN_PP7', 'VAMP2_4nl', 'VAMP2_PP7', 'DNM1_4nl', 'DNM1_PP7', 'mVAMP2_HK16', 'MAPP_JK100', 'mVAMP2_HZ120', 'Spike_in_RNA'};

%process MAPseq data to produce barcode matrix
gene_list={'AP2B1', 'AP2M1', 'CACNG4', 'CADM1', 'DNM1', 'ERC2', 'GAP43', 'PTPRA', 'RAB3A', 'RAB6B', 'SNAP25', 'STX1B', 'STX6', 'SYNAPSINI', 'VAMP2_SYNAPTOPHYSIN', 'SYNAPTOPHYSIN', 'SYNGR3', 'SYT1', 'VAMP1', 'VAMP2', 'VAMP4', 'STX1B_4nl', 'STX1B_PP7', 'SYNAPTOPHYSIN_4nl', 'SYNAPTOPHYSIN_PP7', 'VAMP2_4nl', 'VAMP2_PP7', 'DNM1_4nl', 'DNM1_PP7', 'mVAMP2_HK16', 'MAPP_JK100'};
for k = 1:length(gene_list)
    j=gene_list{k};
    j
    
    bcn=[1,2,3,4,53];
    %bcn=1:length(bcn);

    %filter out SSI with no barcodes
    filesize=[];
    for i=1:length(bcn)
        file=dir([brain,prefix,'_',int2str(bcn(i)),'_',j,'_counts.txt']);
        %if ~isempty(file)
            %filesize(i) = file.bytes;  % If file exists, get its size
        %else
            % If file does not exist, set size to 0 or skip depending on your needs
            %filesize(i) = 0;  % Alternatively, you can skip this iteration
            %warning(['File not found: thresholds/', prefix, '_', int2str(bcn(i)), '_', j, '_counts.txt']);
        %end
        filesize(i)=file.bytes;
    end

    bcnfilt=bcn(filesize>0);
    save([brain,'banfilt','_',j,'.mat'],'bcnfilt')
    data=struct();
    for i=1:length(bcnfilt)
    data(i).counts=dlmread([brain,prefix,'_',int2str(bcnfilt(i)),'_',j,'_counts.txt']);
    data(i).reads=int8(char(textread([brain,prefix,'_',int2str(bcnfilt(i)),'_',j,'_seq.txt'],'%s')));
    end
    
    %save(['data_',brain,'_',j,'.mat'],'data','-v7.3');
 
 
 
 
    %% Finish error correction by reading in bowtie alignments and finding connected graph components
    
    %load alignments
    positions1=[];
    for i=1:length(bcnfilt)
    positions1(i).x=dlmread([brain,'bowtie',int2str(bcnfilt(i)),'_',j,'_2u_1.txt']);
    positions1(i).y=dlmread([brain,'bowtie',int2str(bcnfilt(i)),'_',j,'_2u_3.txt']);
    clustermatrix1(i).C=sparse(positions1(i).x,positions1(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
    end
    %save(['clustermatrix1','_',j,'.mat'],'clustermatrix1','-v7.3');
    % 
    %load clustermatrix1
    
    %find connected components
    graph=[];
    for i=1:length(bcnfilt)
    i
        [graph(i).S,graph(i).G]=graphconncomp(clustermatrix1(i).C,'Directed',false); %find the connected graph components
    end
    % save(['graph1','_',j,'.mat'],'graph');
    % 
    % load graph1
    
    %collapse barcodes to most abundant member of the connected graph component
    
    for i=1:length(bcnfilt)
    x=1:graph(i).S;
    [tf,loc]=ismember(x,graph(i).G,'R2012a');
    collapsedreads=data(i).reads(loc,:);
    collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
    [corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
    corrected(i).reads2u=collapsedreads(ix,:);
    end
    
    
    %% remove reads containing homopolymers
    minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
    for i=1:length(bcnfilt)
        a=findhomopolymers(corrected(i).reads2u,minrunlength);
        corrected(i).freads=corrected(i).reads2u(~a,:);
        corrected(i).fcounts=corrected(i).counts2u(~a,:);
    end
    
    for i=1:length(bcnfilt)
    data(i).BCseqff=corrected(i).freads;
    data(i).BCcountsff=corrected(i).fcounts;
    end
    %save('data.mat','data','-v7.3');
    
 
 
 
 
      %{ 
        %% do collapse for spike-in molecules
        
                                
        %read in raw spike in sequences
        %bcnfilt=[10:16,18:23,25];
        for i=1:length(bcnfilt)
        spikes(i).counts=dlmread(['thresholds/',prefix,'_','spikes',int2str(bcnfilt(i)),'_counts.txt']);
        spikes(i).reads=int8(char(textread(['thresholds/',prefix,'_','spikes',int2str(bcnfilt(i)),'_seq.txt'],'%s')));
        end
                           
                           
        %read in bowtie alignments of spikes
        positions2=[];
        
        for i=1:length(bcnfilt)
        positions2(i).x=dlmread(['thresholds/','bowtiespikes',int2str(bcnfilt(i)),'_2u_1.txt']);
        positions2(i).y=dlmread(['thresholds/','bowtiespikes',int2str(bcnfilt(i)),'_2u_3.txt']);
        clustermatrix2(i).C=sparse(positions2(i).x,positions2(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
        end
        %save('clustermatrix2.mat','clustermatrix2','-v7.3');
        
        %load clustermatrix2
        graph2=[];
        for i=1:length(bcnfilt)
        i
            [graph2(i).S,graph2(i).G]=graphconncomp(clustermatrix2(i).C,'Directed',false); %find the connected graph components
        end
        %save('graph2.mat','graph2');
        
        
        %load graph2
        for i=1:length(bcnfilt)
        x=1:graph2(i).S;
        [tf,loc]=ismember(x,graph2(i).G,'R2012a');
        collapsedreads=spikes(i).reads(loc,:);
        collapsedcounts=accumarray(graph2(i).G',spikes(i).counts);%'
        [spikes(i).counts2u,ix]=sort(collapsedcounts,'descend');
        spikes(i).reads2u=collapsedreads(ix,:);
        end
        
        
        save(['spikes',prefix,'_',brain,'.mat'],'spikes')
      %}
        %% load data.m and match the barcodes from the different areas against each other to produce a barcode matrix
        % this is code modified to in situ MAPseq, so it will run without an
        % injection site.
        
        
        %% find overlap of barcodes
        %% set up reference set of barcodes
        %load data
        %load data.mat
        
        %collect all suquences detected in the target sites
    refbarcodes_tmp=[];
    for i=1:length(data)
    refbarcodes_tmp=[refbarcodes_tmp;data(i).BCseqff];  
    end
    refbarcodes=unique(refbarcodes_tmp,'rows'); %unique barocodes to get reference set

    %% construct barcodematrix by matching barcodes in target sites to reference barcodes

    barcodematrix=zeros(size(refbarcodes,1),length(data));%initiate the barcode matrix

    for i=1:length(data)
        %pull out reads and counts into new variables for ease of use
        BCreads=data(i).BCseqff;
        BCmolcounts=data(i).BCcountsff;
        [ind,loc]=ismember(BCreads,refbarcodes,'rows');
        barcodematrix(loc(loc~=0),i)=BCmolcounts(ind);   
    end


    %set thresholds of how many molecule counts a barcode has to have to be
    %considered real. by default we leave this at 0
    threshold=0; %lowest molecule count considered trustworthy
    barcodematrix=barcodematrix(sum(barcodematrix>threshold,2)~=0,:);
    refbarcodes=refbarcodes(sum(barcodematrix>threshold,2)~=0,:);

    save([brain,'barcodematrix',prefix,'_',brain,'_',j,'.mat'],'barcodematrix','refbarcodes','-v7.3')

end						   
