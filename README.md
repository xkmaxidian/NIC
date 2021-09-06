# NIC
Overview:

This is code to do joint learning dimension reduction and clustering of single-cell multiomics data given in the "experiment" section of the paper: 

Wenming Wu, Xiaoke Ma*. "Network-based integrative analysis of single-cell transcriptomic and epigenomic data for cell types".

The coding here is a generalization of the algorithm given in the paper. NIC is written in the MATLAB programming language. To use, please download the NIC folder and follow the instructions provided in the ¡°README.doc¡±.

Files:

clustering_NIC.m - The main function.

main_NIC.m - A script with single cell multiomics data to show how to run the code.

dataset1.mat, dataset2.mat - Two simulated single cell multiomics datasets for the experiment. 

mESC.mat - A real single cell multiomics data utilized in the cell type clustering example. mESC for mouse embryonic stem cells (mESCs), including 13 cells cultured in 2i media and 64 serum-grown cells, which were profiled by parallel single-cell methylation and transcriptome sequencing technique scM&T-seq

DR_nmf.m - NMF is utilized to extract features to preprocess the single cell multiomics data.

constructW.m - Compute adjacent matrix W.

PMI.m - PMI matrix construction.

ClusteringMeasure_new.m - Cell type clustering performance evaluation.

rand_index.m - Program for calculating the Adjusted Rand Index ( Hubert & Arabie).

Contact:

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn 
