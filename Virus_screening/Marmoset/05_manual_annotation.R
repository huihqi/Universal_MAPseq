library(tidyverse)
source("dataset.R")
#source("common.R")


main = function() {
    #annotate_whole()
    annotate_glu()
    annotate_gaba()
    #annotate_subglu()
}

annotate_whole = function(output_dir = "analysis/whole") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "GABA_1",
        2, "GLU_1",
        3, "Low_qual",
        4, "GLU_2",
        5, "GABA_2",
        6, "GLU_3",
        7, "GLU_4",
        8, "GLU_5",
        9, "GABA_3",
        10, "GLU_6",
        11, "GLU_7",
        12, "GLU_8",
        13, "GLU_9",
        14, "GLU_10",
        15, "Cd_1",
        16, "Cd_2",
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_glu = function(output_dir = "analysis/glu") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "TH_1",
        2, "Cl",
        3, "Pul",
        4, "CM/VP",
        5, "Ctx",
        6, "MD_PulA_senosry",
        7, "medial_ventral",
        8, "LGN", 
        9, "HIP",
        10, "Hb",


    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_gaba = function(output_dir = "analysis/gaba") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,
        1, "gaba_Cd_Pu_1","gaba_Cd_Pu",
        2, "gaba_Cd_Pu_2","gaba_Cd_Pu",
        3, "gaba_med_vent","gaba_med_vent",
        4, "gaba_Ctx","gaba_Ctx",
        5, "gaba_TH_1","gaba_TH",
        6, "gaba_TH_2","gaba_TH",
        7, "gaba_TH_3","gaba_TH",
        8, "gaba_Unknown","gaba_Unknown",
        9, "gaba_Rt","gaba_Rt",
     
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}
#99 clusters total, including 4 low quality clusters
annotate_subglu = function(output_dir = "analysis/subglu") {
    subclass_annotation = read_csv("analysis/glu/cluster_annotation.csv")
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,~notes,
        "Ctx_1","Ctx_1","Ctx"," ",
        "Ctx_2","Ctx_2","Ctx"," ",
        "Ctx_3","Ctx_3","Ctx"," ",
        "Ctx_4","Ctx_4","Ctx"," ",
        "Ctx_5","Ctx_5","Ctx"," ",
        "Ctx_6","Ctx_6","Ctx"," ",
        "Ctx_7","Ctx_7","Ctx"," ",
        "Ctx_8","Ctx_8","Ctx"," ",
        
        "Hb_1","Hb_1","Hb"," ",
        "Hb_2","Hb_2","Hb"," ",
        "Hb_3","Hb_3","Hb"," ",
        "Hb_4","Hb_4","Hb"," ",
       
        "HIP_1","HIP_1","HIP"," ",
        "HIP_2","HIP_2","HIP"," ",
        "HIP_3","HIP_3","HIP"," ",
        "HIP_4","HIP_4","HIP"," ",
        "HIP_5","HIP_5","HIP"," ",
        
        "LD_1","LD_1","LD"," ",
        "LD_2","LD_2","LD"," ",
        "LD_3","LD_3","LD"," ",
        "LD_4","LD_4","LD"," ",
        "LD_5","LD_5","LD"," ",
        
        
        "LGN_1","LGN_1","LGN"," ",
        "LGN_2","LGN_2","LGN"," ",
        "LGN_3","LGN_3","LGN"," ",
        "LGN_4","LGN_4","LGN"," ",
        "LGN_5","LGN_5","LGN"," ",
        "LGN_6","LGN_6","LGN"," ",
        "LGN_7","LGN_7","LGN"," ",
        
        "MD_PulA_senosry_1","MD_PulA_sensory_1","MD_PulA_sensory"," ",
        "MD_PulA_senosry_2","MD_PulA_sensory_2","MD_PulA_sensory"," ",
        "MD_PulA_senosry_3","MD_PulA_sensory_3","MD_PulA_sensory"," ",
        "MD_PulA_senosry_4","MD_PulA_sensory_4","MD_PulA_sensory"," ",
        "MD_PulA_senosry_5","MD_PulA_sensory_5","MD_PulA_sensory"," ",
        
        "Medioventral_1","Medioventral_1","Medioventral"," ",       
        "Medioventral_2","Medioventral_2","Medioventral"," ",       
        "Medioventral_3","Medioventral_3","Medioventral"," ",       
        "Medioventral_4","Medioventral_4","Medioventral"," ",       
        "Medioventral_5","Medioventral_5","Medioventral"," ",       
        "Medioventral_6","Medioventral_6","Medioventral"," ",    
        
        "mRT_ZI_1","mRT_ZI_1","mRT_ZI"," ",         
        "mRT_ZI_2","mRT_ZI_2","mRT_ZI"," ",         
        "mRT_ZI_3","mRT_ZI_3","mRT_ZI"," ",         
        "mRT_ZI_4","mRT_ZI_4","mRT_ZI"," ",         
        "mRT_ZI_5","mRT_ZI_5","mRT_ZI"," ",         
        "mRT_ZI_6","mRT_ZI_6","mRT_ZI"," ",         
        "mRT_ZI_7","mRT_ZI_7","mRT_ZI"," ",       
        
        
        "Pul_1","Pul_1","Pul"," ",
        "Pul_2","Pul_2","Pul"," ",
        "Pul_3","Pul_3","Pul"," ",
        "Pul_4","Pul_4","Pul"," ",
        "Pul_5","Pul_5","Pul"," ",
        "Pul_6","Pul_6","Pul"," ",
        
        "mRT_ZI_1","mRT_ZI_1","mRT_ZI"," ",         
        "mRT_ZI_2","mRT_ZI_2","mRT_ZI"," ",         
        "mRT_ZI_3","mRT_ZI_3","mRT_ZI"," ",         
        "mRT_ZI_4","mRT_ZI_4","mRT_ZI"," ",         
        "mRT_ZI_5","mRT_ZI_5","mRT_ZI"," ",         
        "mRT_ZI_6","mRT_ZI_6","mRT_ZI"," ",         
        "mRT_ZI_7","mRT_ZI_7","mRT_ZI"," ",         
        "mRT_ZI_8","mRT_ZI_8","mRT_ZI"," ",         
        
        "TH_1_1","TH_Low_qual_1","TH_Low_qual"," ",
        "TH_1_2","TH_Low_qual_2","TH_Low_qual"," ",
        "TH_1_3","TH_Low_qual_3","TH_Low_qual"," ",
        "TH_1_4","TH_Low_qual_4","TH_Low_qual"," ",
        "TH_1_5","TH_Low_qual_5","TH_Low_qual"," ",
        "TH_1_6","TH_Low_qual_6","TH_Low_qual"," ",
        "TH_1_7","TH_Low_qual_7","TH_Low_qual"," ",
        "TH_1_8","TH_Low_qual_8","TH_Low_qual"," ",
      

    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

if (sys.nframe() == 0) {
    main()
}

