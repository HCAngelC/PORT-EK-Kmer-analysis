#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2024.12.25
##############################################################
# Input: A full list of k-mers associated with corresponding loci using bootstrapping method
# Object: A matrix of the enrichment of k-mer in corresponding SARS-CoV-2 genes
##############################################################

# R functions

## The deers exclusive group
bootstrapping_kmer_over_genes_100 <- function(df) {
  gene_btsp <- sample_n(df, 100, replace = T)
  
  gene_btsp.A <- dplyr::filter(gene_btsp, group == "A")
  df_out_A <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_enriched_bd(gene_btsp.A) %>% dplyr::mutate(group = "A")
  
  gene_btsp.B <- dplyr::filter(gene_btsp, group == "B")
  df_out_B <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_favourable_bd(gene_btsp.B) %>% dplyr::mutate(group = "B")
  
  gene_btsp.C <- dplyr::filter(gene_btsp, group == "C")
  df_out_C <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched_bd(gene_btsp.B) %>% dplyr::mutate(group = "C")
  
  gene_btsp.D <- dplyr::filter(gene_btsp, group == "D")
  df_out_D <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable_bd(gene_btsp.C) %>% dplyr::mutate(group = "D")
  
  df_out <- dplyr::bind_rows(df_out_A, df_out_B, df_out_C, df_out_D) 
  df_out <- replace(df_out, is.na(df_out), 0)
  
  return(df_out)
}

bootstrapping_kmer_over_genes_500 <- function(df) {
  gene_btsp <- sample_n(df, 500, replace = T)
  
  gene_btsp.A <- dplyr::filter(gene_btsp, group == "A")
  df_out_A <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_enriched_bd(gene_btsp.A) %>% dplyr::mutate(group = "A")
  
  gene_btsp.B <- dplyr::filter(gene_btsp, group == "B")
  df_out_B <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_favourable_bd(gene_btsp.B) %>% dplyr::mutate(group = "B")
  
  gene_btsp.C <- dplyr::filter(gene_btsp, group == "C")
  df_out_C <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched_bd(gene_btsp.B) %>% dplyr::mutate(group = "C")
  
  gene_btsp.D <- dplyr::filter(gene_btsp, group == "D")
  df_out_D <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable_bd(gene_btsp.C) %>% dplyr::mutate(group = "D")
  
  df_out <- dplyr::bind_rows(df_out_A, df_out_B, df_out_C, df_out_D) 
  df_out <- replace(df_out, is.na(df_out), 0)
  
  return(df_out)
}

bootstrapping_kmer_over_genes_1000 <- function(df) {
  gene_btsp <- sample_n(df, 1000, replace = T)
  
  gene_btsp.A <- dplyr::filter(gene_btsp, group == "A")
  df_out_A <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_enriched_bd(gene_btsp.A) %>% dplyr::mutate(group = "A")
  
  gene_btsp.B <- dplyr::filter(gene_btsp, group == "B")
  df_out_B <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_favourable_bd(gene_btsp.B) %>% dplyr::mutate(group = "B")
  
  gene_btsp.C <- dplyr::filter(gene_btsp, group == "C")
  df_out_C <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched_bd(gene_btsp.B) %>% dplyr::mutate(group = "C")
  
  gene_btsp.D <- dplyr::filter(gene_btsp, group == "D")
  df_out_D <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable_bd(gene_btsp.C) %>% dplyr::mutate(group = "D")
  
  df_out <- dplyr::bind_rows(df_out_A, df_out_B, df_out_C, df_out_D) 
  df_out <- replace(df_out, is.na(df_out), 0)
  
  return(df_out)
}

bootstrapping_kmer_over_genes_5000 <- function(df) {
  gene_btsp <- sample_n(df, 5000, replace = T)
  
  gene_btsp.A <- dplyr::filter(gene_btsp, group == "A")
  df_out_A <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_enriched_bd(gene_btsp.A) %>% dplyr::mutate(group = "A")
  
  gene_btsp.B <- dplyr::filter(gene_btsp, group == "B")
  df_out_B <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_favourable_bd(gene_btsp.B) %>% dplyr::mutate(group = "B")
  
  gene_btsp.C <- dplyr::filter(gene_btsp, group == "C")
  df_out_C <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched_bd(gene_btsp.B) %>% dplyr::mutate(group = "C")
  
  gene_btsp.D <- dplyr::filter(gene_btsp, group == "D")
  df_out_D <- calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable_bd(gene_btsp.C) %>% dplyr::mutate(group = "D")
  
  df_out <- dplyr::bind_rows(df_out_A, df_out_B, df_out_C, df_out_D) 
  df_out <- replace(df_out, is.na(df_out), 0)
  
  return(df_out)
}

bootstrap_deer_breakdown_100 <- do.call(rbind, replicate(5000, bootstrapping_kmer_over_genes_100(df_matrix_wszystko_deer_breakdown), simplify = F))
bootstrap_deer_breakdown_500 <- do.call(rbind, replicate(5000, bootstrapping_kmer_over_genes_500(df_matrix_wszystko_deer_breakdown), simplify = F))
bootstrap_bat_breakdown_1000 <- do.call(rbind, replicate(5000, bootstrapping_kmer_over_genes_1000(df_matrix_wszystko_bat_breakdown), simplify = F))
bootstrap_bat_breakdown_5000 <- do.call(rbind, replicate(5000, bootstrapping_kmer_over_genes_5000(df_matrix_wszystko_bat_breakdown), simplify = F))
