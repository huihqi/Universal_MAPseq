library(tidyverse)
source("dataset.R")


main = function() {
    #annotate_whole()
    annotate_glu()
    annotate_gaba()
    #annotate_subglu()
}

annotate_whole = function(output_dir = "analysis/whole") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "GLU_1",
        2, "GABA_2",
        3, "GABA_3",
        4, "GABA_4",
        5, "GABA_5",
        6, "GLU_6",
        7, "GLU_7",
        8, "GLU_8",
        9, "GLU_9",
        10, "GABA_10",
        11, "GLU_11",
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_glu = function(output_dir = "analysis/glu") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "L2_IT",
        2, "L3_IT",
        3, "L4_IT",
        4, "L6CT_L6b",
        5, "Car3",
        6, "glu_6",
        7, "L6_IT",
        8, "L5_IT",
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_gaba = function(output_dir = "analysis/gaba") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,
        1, "gaba_1","gaba_1",
        2, "gaba_2","gaba_2",
        3, "gaba_3","gaba_3",
        4, "gaba_4","gaba_4",
        5, "gaba_5","gaba_5",
        6, "gaba_6","gaba_6",
        7, "gaba_7","gaba_7",
        8, "gaba_8","gaba_8",
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}
#99 clusters total, including 4 low quality clusters
annotate_subglu = function(output_dir = "analysis/subglu") {
    subclass_annotation = read_csv("analysis/glu/cluster_annotation.csv")
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,~notes,
        "Car3_1","Car3 1","Car3"," ",
        "Car3_2","Car3 2","Car3"," ",
        "Car3_3","Car3 3","Car3"," ",
        "Car3_4","Car3 4","Car3"," ",
        
        "glu_1_1","glu_1 1","glu_1"," ",
        "glu_1_2","glu_1 2","glu_1"," ",
        "glu_1_3","glu_1 3","glu_1"," ",
        "glu_1_4","glu_1 4","glu_1"," ",
        "glu_1_5","glu_1 5","glu_1"," ",
        
        "glu_2_1","glu_2 1","glu_2"," ",
        "glu_2_2","glu_2 2","glu_2"," ",
        "glu_2_3","glu_2 3","glu_2"," ",
        "glu_2_4","glu_2 4","glu_2"," ",
        
        "glu_4_1","glu_4 1","glu_4"," ",
        "glu_4_2","glu_4 2","glu_4"," ",
        "glu_4_3","glu_4 3","glu_4"," ",
        "glu_4_4","glu_4 4","glu_4"," ",
        "glu_4_5","glu_4 5","glu_4"," ",
        "glu_4_6","glu_4 6","glu_4"," ",
        
        "glu_5_1","glu_5 1","glu_5"," ",
        "glu_5_2","glu_5 2","glu_5"," ",
        "glu_5_3","glu_5 3","glu_5"," ",
        "glu_5_4","glu_5 4","glu_5"," ",
        "glu_5_5","glu_5 5","glu_5"," ",
        "glu_5_6","glu_5 6","glu_5"," ",
        
        "glu_16_1","glu_16 1","glu_16"," ",
        "glu_16_2","glu_16 2","glu_16"," ",
        
        "glu_20_1","glu_20 1","glu_20"," ",
        "glu_20_2","glu_20 2","glu_20"," ",
        "glu_20_3","glu_20 3","glu_20"," ",
        "glu_20_4","glu_20 4","glu_20"," ",
        "glu_20_5","glu_20 5","glu_20"," ",
        "glu_20_6","glu_20 6","glu_20"," ",
        
        "HIP_1_1","HIP_1 1","HIP_1"," ",
        "HIP_1_2","HIP_1 2","HIP_1"," ",
        "HIP_1_3","HIP_1 3","HIP_1"," ",
        "HIP_1_4","HIP_1 4","HIP_1"," ",
        "HIP_1_5","HIP_1 5","HIP_1"," ",
        "HIP_1_6","HIP_1 6","HIP_1"," ",
        "HIP_1_7","HIP_1 7","HIP_1"," ",
        "HIP_1_8","HIP_1 8","HIP_1"," ",
        
        "HIP_2_1","HIP_2 1","HIP_2"," ",
        "HIP_2_2","HIP_2 2","HIP_2"," ",
        "HIP_2_3","HIP_2 3","HIP_2"," ",
        "HIP_2_4","HIP_2 4","HIP_2"," ",
        
        "HIP_3_1","HIP_3 1","HIP_3"," ",
        "HIP_3_2","HIP_3 2","HIP_3"," ",
        "HIP_3_3","HIP_3 3","HIP_3"," ",
        "HIP_3_4","HIP_3 4","HIP_3"," ",
        
        "L2_3_IT_1","L2_3_IT 1","L2_3_IT"," ",
        "L2_3_IT_2","L2_3_IT 2","L2_3_IT"," ",
        "L2_3_IT_3","L2_3_IT 3","L2_3_IT"," ",
        "L2_3_IT_4","L2_3_IT 4","L2_3_IT"," ",
        "L2_3_IT_5","L2_3_IT 5","L2_3_IT"," ",
        "L2_3_IT_6","L2_3_IT 6","L2_3_IT"," ",
        
        "L4_5_IT_1","L4_5_IT 1","L4_5_IT"," ",
        "L4_5_IT_2","L4_5_IT 2","L4_5_IT"," ",
        "L4_5_IT_3","L4_5_IT 3","L4_5_IT"," ",
        "L4_5_IT_4","L4_5_IT 4","L4_5_IT"," ",
        "L4_5_IT_5","L4_5_IT 5","L4_5_IT"," ",
        "L4_5_IT_6","L4_5_IT 6","L4_5_IT"," ",
        
        "L5_ET_1","L5_ET 1","L5_ET"," ",
        "L5_ET_2","L5_ET 2","L5_ET"," ",
        "L5_ET_3","L5_ET 3","L5_ET"," ",
        "L5_ET_4","L5_ET 4","L5_ET"," ",
        "L5_ET_5","L5_ET 5","L5_ET"," ",
        
        "L5_IT_1","L5_IT 1","L5_IT"," ",
        "L5_IT_2","L5_IT 2","L5_IT"," ",
        "L5_IT_3","L5_IT 3","L5_IT"," ",
        "L5_IT_4","L5_IT 4","L5_IT"," ",
        "L5_IT_5","L5_IT 5","L5_IT"," ",
        
        "L6_CT_1","L6_CT 1","L6_CT"," ",
        "L6_CT_2","L6_CT 2","L6_CT"," ",
        "L6_CT_3","L6_CT 3","L6_CT"," ",
        "L6_CT_4","L6_CT 4","L6_CT"," ",
        "L6_CT_5","L6_CT 5","L6_CT"," ",
        "L6_CT_6","L6_CT 6","L6_CT"," ",
        
        "L6_IT_1","L6_IT 1","L6_IT"," ",
        "L6_IT_2","L6_IT 2","L6_IT"," ",
        "L6_IT_3","L6_IT 3","L6_IT"," ",
        "L6_IT_4","L6_IT 4","L6_IT"," ",
        "L6_IT_5","L6_IT 5","L6_IT"," ",
        
        "L6b_1","L6b 1","L6b"," ",
        "L6b_2","L6b 2","L6b"," ",
        "L6b_3","L6b 3","L6b"," ",
        "L6b_4","L6b 4","L6b"," ",
        "L6b_5","L6b 5","L6b"," ",
        "L6b_6","L6b 6","L6b"," ",
        
        "NP_1","NP 1","NP"," ",
        "NP_2","NP 2","NP"," ",
        "NP_3","NP 3","NP"," ",
        "NP_4","NP 4","NP"," ",
        "NP_5","NP 5","NP"," ",
        "NP_6","NP 6","NP"," ",
        
        "RSP_UL_1","RSP_UL 1","RSP_UL"," ",
        "RSP_UL_2","RSP_UL 2","RSP_UL"," ",
        "RSP_UL_3","RSP_UL 3","RSP_UL"," ",
        "RSP_UL_4","RSP_UL 4","RSP_UL"," ",
        "RSP_UL_5","RSP_UL 5","RSP_UL"," ",
        
        "SUB_1","SUB 1","SUB"," ",
        "SUB_2","SUB 2","SUB"," ",
        "SUB_3","SUB 3","SUB"," ",
        "SUB_4","SUB 4","SUB"," ",
        "SUB_5","SUB 5","SUB"," "
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

if (sys.nframe() == 0) {
    main()
}

