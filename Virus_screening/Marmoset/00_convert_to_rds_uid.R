library(tidyverse)
library(hdf5r)
library(SingleCellExperiment)


convert_currentdata = function() {
  sce = convert_v7_filtneurons(file.path("data","filt_neurons_combined.mat"),
                               file.path("data","uid.csv"),
                               file.path("data","brain_name.csv"))
  saveRDS(sce, file.path("data","converted_barseq_data.rds"))    
}

convert_v7_filtneurons = function(mat_file,uid_file,brain_file) {
  # mat_file="filt_neurons_combined.mat"
  #uid_file="uid.csv"
  #brain_file="brain_name.csv"
  h5_file = H5File$new(mat_file)
  #h5_file$ls()
  
  n_genes = length(readDataSet(h5_file[["filt_neurons/expmat/jc"]]))-1
  gene_refs = readDataSet(h5_file[["filt_neurons/genes"]])
  gene_ids = gene_refs$dereference()
  genes = lapply(gene_ids, readDataSet)
  genes = lapply(genes, intToUtf8)
  genes = unlist(genes[1:n_genes])
  
  # updated to read sample_name and brain_name from separate csv files. This 
  # allows string type uid and brain_name to be read in efficiently
  sample_name=read.csv(uid_file,
                       header=FALSE)
  sample_name=unlist(sample_name)
  brain_name=read.csv(brain_file,
                       header=FALSE)
  names(brain_name)=c('brain_name')
  
  position = readDataSet(h5_file[["filt_neurons/pos"]])
  depth = readDataSet(h5_file[["filt_neurons/depth"]])
  angle = readDataSet(h5_file[["filt_neurons/angle"]])
  slice = readDataSet(h5_file[["filt_neurons/slice"]])
  fov_position = readDataSet(h5_file[["filt_neurons/pos40x"]])
  
  expr = Matrix::sparseMatrix(
    i=readDataSet(h5_file[["filt_neurons/expmat/ir"]])+1,
    p=readDataSet(h5_file[["filt_neurons/expmat/jc"]]),
    x=as.vector(readDataSet(h5_file[["filt_neurons/expmat/data"]])),
    dims=c(length(sample_name), length(genes)),
    dimnames = list(sample_name, genes)
  )
  expr = t(expr)
  
  metadata = data.frame(
    slice = as.vector(slice),
    pos_x = position[,1],
    pos_y = position[,2],
    fov_x = fov_position[,1],
    fov_y = fov_position[,2],
    angle = as.vector(angle),
    depth_x = as.vector(angle),
    depth_y = as.vector(angle),
    brain = (brain_name)
  )
  sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
  return(sce)
}

add_ccf = function(filename) {
    barseq = readRDS(filename)
    ccf = read_csv("neurons_registered.csv") %>%
        select(sample_id=id, ccf_name=name, ccf_red=red, ccf_green=green, ccf_blue=blue,
               ccf_x = x_CCF, ccf_y = y_CCF, ccf_z = z_CCF)
    ccf$sample_id = as.character(ccf$sample_id)
    barseq$sample_id = colnames(barseq)
    new_coldata = left_join(as.data.frame(colData(barseq)), ccf, by="sample_id")
    colData(barseq) = DataFrame(new_coldata)
    colnames(barseq) = barseq$sample_id
    saveRDS(barseq, filename)
}

if (sys.nframe() == 0) {
  convert_currentdata()
}

