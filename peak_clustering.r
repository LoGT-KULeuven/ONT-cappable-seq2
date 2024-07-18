#!/usr/bin/env Rscript

## packages
if (!require("dplyr")) {
    install.packages("dplyr",repos = "http://cran.us.r-project.org")
}

library(dplyr)
library(tools)

## read in the .counts.withTR 
mod_tss_peaks <- function(input, strand_s, cov_min = 1, merge_w = 5){
    # add column names. 
    input %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, 
                  prominence = V5, strand_peak = V6, width = V10,
                  start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16) %>%
    # remove irrelevant columns
    dplyr::select(-V4, -V7, -V8, -V11) %>%
    ## cluster by peak
    dplyr::group_by(start_peak, end_peak) %>%
    ## get full coverage (summed over all peak positions)
    dplyr::mutate(full_cov_peak = sum(cov))%>%
    ## only keep the peak summit position(s)
    dplyr::filter(cov == max(cov)) %>%
    ## In case plateau of max(cov) present
    dplyr::mutate(decision_v = ifelse(strand_s == "+", 
                                    min(end_cov), max(end_cov))) %>%
    dplyr::filter(end_cov == decision_v) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(chr, end_cov) %>%
    
    # group by nucleotide sequence (chr)
    dplyr::group_by(chr) %>%
    # index is the upperbound for clustering with the respective peak; index1 are the clusters
    dplyr::mutate(index = lag(end_cov, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>%
    # cluster by group
    dplyr::group_by(chr, index1) %>%
    # new variable full_cov_clust is the full coverage of the grouped peak.
    dplyr::mutate(full_cov_clust = sum(full_cov_peak))%>%
    # Keep as representative of each cluster the peak with maximal coverage
    dplyr::filter(cov == max(cov),
                  cov >= cov_min)%>%
    dplyr::mutate(RPM = 1000000*full_cov_clust/mapped_reads)

}


cluster_peaks <- function(inputFile,strand,clusterw,covm,out){
        counts <- read.table(inputFile)
        peaks <- mod_tss_peaks(counts,strand,merge_w=clusterw,cov_min=covm)
        write.csv(peaks, out)
}

inputFile=commandArgs(TRUE)[1]
strand=commandArgs(TRUE)[2]
clusterw=as.numeric(commandArgs(TRUE)[3])
covm=as.numeric(commandArgs(TRUE)[4])
out=commandArgs(TRUE)[5]
cluster_peaks(inputFile,strand,clusterw,covm,out)
