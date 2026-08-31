
library(tidyverse)
library(SingleCellExperiment)
library(rasterpdf)
source("dataset.R")


main = function() {
    #plot_analysis("analysis/whole", "figs/whole",min_counts = 20)
    plot_analysis("analysis/glu", "figs/glu",min_counts = 20)
    plot_analysis("analysis/gaba", "figs/gaba",min_counts = 20)
    #plot_subglu(min_counts = 20)
    #plot_analysis("analysis/IT45","figs/IT45")
}

plot_analysis = function(input_dir, output_dir, max_cells = 2000000,min_genes=5,min_counts=20) {
    #input_dir="analysis/whole"
    #output_dir="figs/whole"
    #max_cells = 2000000
    #min_genes=5
    #min_counts=10
    
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)    
    barseq = load_dataset(input_dir, max_cells,min_genes=min_genes,min_counts=min_counts)
    if (!("aligned_x" %in% colnames(colData(barseq)))) {
        barseq$aligned_x<-barseq$pos_x
        barseq$aligned_y<-barseq$pos_y
    }

    plot_umap(barseq, output_dir,filename="preliminary_analysis_4.pdf")
    plot_gradients(barseq, output_dir)
    plot_marker_scores(barseq, output_dir)
    plot_slice_summary(barseq, output_dir)
    plot_slices(barseq, output_dir)
    plot_cell_types(barseq, output_dir)
}

plot_subglu = function(input_dir = "analysis/subglu", output_dir = "figs/subglu",min_genes=5,min_counts=20) {
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)
    barseq = load_dataset(input_dir,min_genes=min_genes,min_counts=min_counts)
    if (!("aligned_x" %in% colnames(colData(barseq)))) {
      barseq$aligned_x<-barseq$pos_x
      barseq$aligned_y<-barseq$pos_y
    }
    
    
    
    plot_umap(barseq, output_dir,1)
    plot_gradients(barseq, output_dir)

    ct_dir = file.path(output_dir, "cell_types")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)    
    reducedDim(barseq, "UMAP") = as.data.frame(read_csv(file.path(input_dir, "subumap.csv")))
    subclass_id = read_csv("analysis/glu/cluster.csv")
    subclass_id = subclass_id[match(colnames(barseq), subclass_id$sample),]$label
    subclass_annotation = deframe(read_csv("analysis/glu/cluster_annotation.csv"))
    subclass_name = subclass_annotation[subclass_id]
    for (my_subclass in sort(unique(subclass_name))) {
        filename = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', my_subclass)))
        subdata = barseq
        subdata$label[subclass_name!=my_subclass] = NA
        plot_umap(subdata[, !is.na(subdata$label)], ct_dir, filename=paste0(filename, "_umap.pdf"))
        plot_gradients(subdata[, !is.na(subdata$label)], ct_dir, filename=paste0(filename, "_gradients.pdf"))
        plot_slice_summary(subdata, ct_dir, filename=paste0(filename, "_slice.pdf"))
    }
}

 load_dataset = function(input_dir, max_cells=2000000,min_genes=5,min_counts=20) {
    barseq = load_barseq(normalization_factor = 10,min_genes=min_genes,min_counts=min_counts)
    barseq$slice = as.factor(barseq$slice)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    barseq$n_genes = colSums(counts(barseq)>0)
    
    clusters = read_csv(file.path(input_dir, "cluster.csv"))
    cluster_annotation = deframe(read_csv(file.path(input_dir, "cluster_annotation.csv")))
    barseq = barseq[, as.character(clusters$sample)]
    barseq$cluster_id = as.factor(clusters$label)
    colLabels(barseq) = cluster_annotation[as.vector(barseq$cluster_id)]
    reducedDim(barseq, "UMAP") = as.data.frame(read_csv(file.path(input_dir, "umap.csv")))
    
    if (ncol(barseq) > max_cells) {
        set.seed(17)
        cells_to_keep = sample.int(ncol(barseq), max_cells)
        barseq = barseq[, cells_to_keep]
    }
    return(barseq)
}

 plot_umap = function(barseq, output_dir, text_size=3, filename="preliminary_analysis.pdf") {
   point_size = 70 / sqrt(ncol(barseq))
   pdf(file.path(output_dir, filename))
   # try alpha = 0.05
   text_name = ifelse("subclass" %in% names(colData(barseq)), "subclass", "label")
   p = plot_raster_umap(barseq, color_by="label", text_by=text_name, text_size=text_size, point_size=point_size) +
     ggtitle("Preliminary clustering") +
     guides(col="none") +
     coord_fixed()
   print(p)
   # p = plot_raster_umap(barseq, color_by="brain_name", text_by=text_name, text_size=text_size, point_size=point_size) +
   #   ggtitle("Brains") +
   #   guides(col="none") +
   #   coord_fixed() + 
   #   theme(legend.position = "right",
   #         legend.text = element_text(size=10)) +
   #   labs(color = "brain_name")
   # print(p)
   n_cells <- ncol(barseq)
   sample_size <- ceiling(n_cells * 0.05)  # 5% of cells
   sampled_cells <- sample(n_cells, sample_size)
   
   p = plot_raster_umap(barseq[,sampled_cells], color_by="brain_name", text_by=text_name, text_size=text_size, point_size = 0.1, point_alpha = 0.3) +
     ggtitle("Brains") +
     guides(colour = guide_legend(title = "Brain Names")) +  # explicitly add and title the legend
     theme(legend.position = "right",  # position the legend
           legend.text = element_text(size=8)) +  # adjust text size if needed
     guides(color = guide_legend(override.aes = list(size=5))) +
     coord_fixed()
   print(p)
   
   p = plot_raster_umap(barseq, color_by = "n_genes", point_size = 0.1, point_alpha = 0.3) +
     ggtitle("Number of detected genes") +
     coord_fixed()
   print(p)
   p = plot_raster_umap(barseq, color_by = "SLC17A7", point_size = 0.1, point_alpha = 0.3) +
     ggtitle("Excitatory marker") +
     coord_fixed()
   print(p)
   p = plot_raster_umap(barseq, color_by = "GAD2", point_size = 0.1, point_alpha = 0.3) +
     ggtitle("Inhibitory marker") +
     coord_fixed()
   print(p)
   stats = tibble(cell_type = barseq$label, n_genes=barseq$n_genes, gaba = logcounts(barseq)["GAD2",], glu = logcounts(barseq)["SLC17A7",])
   p = ggplot(stats, aes(cell_type, n_genes)) +
     geom_boxplot() +
     theme_classic() +
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x=NULL, y="# genes detected")
   print(p)
   p = ggplot(stats, aes(cell_type, glu)) +
     geom_boxplot() +
     theme_classic() +
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x=NULL, y="SLC17A7 expression")
   print(p)
   p = stats %>%
     group_by(cell_type) %>%
     summarize(glu = mean(glu>0)) %>%
     ggplot(aes(cell_type, glu)) +
     geom_col() +
     theme_classic() +
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x=NULL, y="Fraction cells expressing SLC17A7") +
     geom_hline(yintercept = 0.95, linetype = "dashed")
   print(p)
   p = ggplot(stats, aes(cell_type, gaba)) +
     geom_boxplot() +
     theme_classic() +
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x=NULL, y="GAD2 expression")
   print(p)
   p = stats %>%
     group_by(cell_type) %>%
     summarize(gaba = mean(gaba>0)) %>%
     ggplot(aes(cell_type, gaba)) +
     geom_col() +
     theme_classic() +
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x=NULL, y="Fraction cells expressing GAD2") +
     geom_hline(yintercept = 0.95, linetype = "dashed")
   print(p)
   dev.off()
 }
 

plot_raster_umap = function(barseq, color_by, text_by=NULL, point_size=0.1, point_alpha = 0.3, text_size=3, raster.dpi=150) {
    to_plot = as_tibble(reducedDim(barseq, "UMAP"))
    colnames(to_plot) = c("UMAP_1", "UMAP_2")
    if (color_by %in% rownames(barseq)) {
        to_plot[[color_by]] = logcounts(barseq)[color_by,]
    } else {
        to_plot[[color_by]] = barseq[[color_by]]
    }
    result=to_plot %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2, col=.data[[color_by]])) +
        ggrastr::geom_point_rast(size=point_size, alpha = point_alpha, shape=16, stroke=0, raster.dpi = raster.dpi) +
        theme_classic()
    if (is.numeric(to_plot[[color_by]])) {
        result = result + scale_color_viridis_c()
    }
    if (!is.null(text_by)) {
        to_plot[[text_by]] = barseq[[text_by]]
        to_plot = to_plot %>%
            group_by(.data[[text_by]]) %>%
            summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        result = result +
            geom_text(data=to_plot, aes(label=.data[[text_by]]), size=text_size, color="black")
    }
    return(result)
}

plot_gradients = function(barseq, output_dir, filename="gradients.pdf") {
    set.seed(17)
    point_size = 70 / sqrt(ncol(barseq))
    barseq$depth = sqrt(barseq$depth_x^2 + barseq$depth_y^2)
    barseq$slice = as.numeric(barseq$slice)
    barseq = barseq[, sample.int(ncol(barseq))]
    pdf(file.path(output_dir, filename))
    p = plot_raster_umap(barseq, color_by = "slice", point_size=point_size) +
        ggtitle("Anterior-Posterior") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Slice") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "aligned_x", point_size=point_size) +
        ggtitle("Medial-Lateral") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Medial->Lateral") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "aligned_y", point_size=point_size) +
        ggtitle("Dorso-ventral") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Ventral->Dorsal") +
        coord_fixed()
    print(p)
    dev.off()
}

plot_marker_scores = function(barseq, output_dir) {
    marker_scores = load_marker_scores(barseq)
    MetaMarkers::plot_marker_scores(
        marker_scores, reducedDim(barseq, "UMAP"),
        normalize_scores=FALSE, rasterize_umap=TRUE, point_size=0.1
    ) +
        theme_classic() +
        ggtitle("Marker score")
    ggsave(file.path(output_dir, "marker_scores.pdf"), width = 15, height=15)
}

load_marker_scores = function(barseq) {
  marker_scores = read_csv("analysis/marker_scores.csv") %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
  marker_scores = t(marker_scores) / colMaxs(marker_scores,useNames = FALSE)
  marker_scores = marker_scores[,as.character(colnames(barseq))]
  return(marker_scores)
}

plot_slice_summary = function(barseq, output_dir, filename="all_slices.pdf") {    
    # n_slices = max(as.numeric(unfactor(barseq$slice)),na.rm=TRUE)
    slice_numbers = unique(as.numeric(unfactor(barseq$slice)))
    slice_numbers=sort(slice_numbers[!is.na(slice_numbers)])
    brain_names=unique(as.character(barseq$brain_name))
    pdf(file.path(output_dir, filename))
    
    for (brain in brain_names) {
      # Loop through unique slices for this brain
      for (slice_id in slice_numbers) {
        # Only plot if there are cells in this slice for this brain
        if (any(barseq$slice == slice_id & barseq$brain_name == brain)) {
          print(plot_slice_brain(barseq, 
                                 slice_id = slice_id, 
                                 brain = brain))
        }
      }
    }
    
    # for (i in 1:n_slices) {
    #     print(plot_slice(barseq, i))
    # }
    dev.off()
}



plot_slices = function(barseq, output_dir) {
    slice_dir = file.path(output_dir, "slices")
    dir.create(slice_dir, showWarnings = FALSE, recursive=TRUE)
    slice_numbers = unique(as.numeric(unfactor(barseq$slice)))
    slice_numbers=sort(slice_numbers[!is.na(slice_numbers)])
    brain_names=unique(as.character(barseq$brain_name))
    
    # n_slices = max(as.numeric(unfactor(barseq$slice)))
    unique_labels = levels(as.factor(barseq$label))
    for (brain in brain_names) {
      # Loop through unique slices for this brain
      for (slice_id in slice_numbers) {
        # Only plot if there are cells in this slice for this brain
        if (any(barseq$slice == slice_id & barseq$brain_name == brain)) {
          pdf(file.path(slice_dir, paste0(brain,"slice_", slice_id, ".pdf")))
          print(plot_slice_brain(barseq, 
                                 slice_id = slice_id, 
                                 brain = brain))
          for (l in unique_labels) {
            print(plot_slice_ct(barseq, slice_id, l,brain))
          }
          dev.off()
        }
      }
    }
    
    
    # for (i in 1:n_slices) {
    #     pdf(file.path(slice_dir, paste0("slice_", i, ".pdf")))
    #     print(plot_slice(barseq, i))
    #     for (l in unique_labels) {
    #         print(plot_slice_ct(barseq, i, l))
    #     }
    #     dev.off()
    # }
}

plot_cell_types = function(barseq, output_dir) {
    ct_dir = file.path(output_dir, "cell_types")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)
    slice_numbers = unique(as.numeric(unfactor(barseq$slice)))
    slice_numbers=sort(slice_numbers[!is.na(slice_numbers)])
    brain_names=unique(as.character(barseq$brain_name))
    
    # n_slices = max(as.numeric(unfactor(barseq$slice)))
    for (l in unique(barseq$label)) {
        file_name = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', l)))
        pdf(file.path(ct_dir, paste0(file_name, ".pdf")))
        for (brain in brain_names) {
          # Loop through unique slices for this brain
          for (slice_id in slice_numbers) {
            # Only plot if there are cells in this slice for this brain
            if (any(barseq$slice == slice_id & barseq$brain_name == brain)) {
              print(plot_slice_ct(barseq, slice_id, l,brain))
            }
          }
        }
        dev.off()
    }
}

plot_slice = function(barseq, slice_id, raster.dpi=150) {
    as.data.frame(colData(barseq)) %>%
        mutate(label = as.factor(label)) %>%
        filter(slice == slice_id) %>%
        ggplot(aes(x = aligned_x, y = aligned_y, col = label)) +
        ggrastr::geom_point_rast(size=0.5, shape=16, stroke=0, raster.dpi = raster.dpi) +
        theme_classic() +
        labs(x=NULL, y=NULL, col=NULL, title=paste0("Slice ", slice_id)) +
        coord_equal() +
        lims(x = range(barseq$aligned_x), y = range(barseq$aligned_y)) +
        theme(legend.text=element_text(size=6), legend.position = "bottom") +
        scale_color_discrete(na.value="gray90", drop=FALSE)
}

plot_slice_brain = function(barseq, slice_id, brain, raster.dpi=300) {
  # as.data.frame(colData(barseq)) %>%
  #   mutate(label = as.factor(label)) %>%
  #   filter(slice == slice_id,
  #          brain_name == brain) %>%
  #   ggplot(aes(x = aligned_x, y = aligned_y, col = label)) +
  #   ggrastr::geom_point_rast(size=0.5, shape=16, stroke=0, raster.dpi = raster.dpi) +
  #   theme_classic() +
  #   labs(x=NULL, y=NULL, col=NULL, 
  #        title=paste0("Brain: ", brain, " - Slice ", slice_id)) +
  #   coord_equal() +
  #   lims(x = range(barseq$aligned_x), y = range(barseq$aligned_y)) +
  #   theme(legend.text=element_text(size=6), legend.position = "bottom") +
  #   scale_color_discrete(na.value="gray90", drop=FALSE)
  filtered_data <- as.data.frame(colData(barseq)) %>%
    filter(slice == slice_id,
           brain_name == brain)
  
  filtered_data %>%
    mutate(label = as.factor(label)) %>%
    ggplot(aes(x = aligned_x, y = aligned_y, col = label)) +
    ggrastr::geom_point_rast(size=0.5, shape=16, stroke=0, raster.dpi = raster.dpi) +
    theme_classic() +
    labs(x=NULL, y=NULL, col=NULL, 
         title=paste0("Brain: ", brain, " - Slice ", slice_id)) +
    coord_equal() +
    lims(x = range(filtered_data$aligned_x), 
         y = range(filtered_data$aligned_y)) +
    theme(legend.text=element_text(size=6), legend.position = "bottom") +
    scale_color_discrete(na.value="gray90", drop=FALSE)
  
  
  
}




plot_slice_ct = function(barseq, slice_id, ct,brain) {
    all_labels = as.factor(barseq$label)
    my_cols = make_cluster_cols(all_labels)
    full_plot = plot_slice_brain(barseq, slice_id,brain)
    my_cols[names(my_cols) != ct] = "gray90"
    my_plot = full_plot +
        ggtitle(paste0(brain,"_Slice ", slice_id, ", Cluster ", ct)) +
        scale_color_manual(values = my_cols, drop=FALSE)
    return(my_plot)
}

make_cluster_cols = function(cluster_names) {
    cluster_names = levels(cluster_names)
    result = scales::hue_pal()(length(cluster_names))
    names(result) = cluster_names
    return(result)
}

if (sys.nframe() == 0) {
    main()
}

