#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2024.02.28
##############################################################
# Input: A full list of k-mers associated with corresponding loci (complete ORF1ab)
# Object: A matrix of the enrichment of k-mer in corresponding SARS-CoV-2 genes
##############################################################

# R functions

## The deers exclusive group
calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_deer_enriched <- function(df) {

ORF1ab <- df %>% dplyr::filter(gene == "ORF1ab")
S <- df %>% dplyr::filter(gene == "S")
ORF3a <- df %>% dplyr::filter(gene == "ORF3a")
E <- df %>% dplyr::filter(gene == "E")
M <- df %>% dplyr::filter(gene == "M")
ORF6 <- df %>% dplyr::filter(gene == "ORF6")
ORF7a <- df %>% dplyr::filter(gene == "ORF7a")
ORF7b <- df %>% dplyr::filter(gene == "ORF7b")
ORF8 <- df %>% dplyr::filter(gene == "ORF8")
N <- df %>% dplyr::filter(gene == "N")
ORF10 <- df %>% dplyr::filter(gene == "ORF10")

zamowienie_corona_gene <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum(ORF1ab$n_enriched)[1]/dim(ORF1ab)[1])/(21555-266+1)), ((sum(S$n_enriched)[1]/dim(S)[1])/(25384-21563+1)), ((sum(ORF3a$n_enriched)[1]/dim(ORF3a)[1])/(26220-25393+1)), ((sum(E$n_enriched)[1]/dim(E)[1])/(26472-26245+1)), ((sum(M$n_enriched)[1]/dim(M)[1])/(27191-26523+1)), ((sum(ORF6$n_enriched)[1]/dim(ORF6)[1])/(27387-27202+1)),
((sum(ORF7a$n_enriched)[1]/dim(ORF7a)[1])/(27759-27394+1)), ((sum(ORF7b$n_enriched)[1]/dim(ORF7b)[1])/(27887-27756+1)), ((sum(ORF8$n_enriched)[1]/dim(ORF8)[1])/(28259-27894+1)),
((sum(N$n_enriched)[1]/dim(N)[1])/(29533-28274+1)),
((sum(ORF10$n_enriched)[1]/dim(ORF10)[1])/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

## The deers favorable group
calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_deer_favourable <- function(df) {

ORF1ab <- df %>% dplyr::filter(gene == "ORF1ab")
S <- df %>% dplyr::filter(gene == "S")
ORF3a <- df %>% dplyr::filter(gene == "ORF3a")
E <- df %>% dplyr::filter(gene == "E")
M <- df %>% dplyr::filter(gene == "M")
ORF6 <- df %>% dplyr::filter(gene == "ORF6")
ORF7a <- df %>% dplyr::filter(gene == "ORF7a")
ORF7b <- df %>% dplyr::filter(gene == "ORF7b")
ORF8 <- df %>% dplyr::filter(gene == "ORF8")
N <- df %>% dplyr::filter(gene == "N")
ORF10 <- df %>% dplyr::filter(gene == "ORF10")

zamowienie_corona_gene <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum(ORF1ab$n_enriched)[1]/dim(ORF1ab)[1])/(21555-266+1)-((sum(ORF1ab$n_depleted)[1]/dim(ORF1ab)[1])/(21555-266+1))), ((sum(S$n_enriched)[1]/dim(S)[1])/(25384-21563+1))-((sum(S$n_depleted)[1]/dim(S)[1])/(25384-21563+1)), ((sum(ORF3a$n_enriched)[1]/dim(ORF3a)[1])/(26220-25393+1))-((sum(ORF3a$n_depleted)[1]/dim(ORF3a)[1])/(26220-25393+1)), ((sum(E$n_enriched)[1]/dim(E)[1])/(26472-26245+1))-((sum(E$n_depleted)[1]/dim(E)[1])/(26472-26245+1)), ((sum(M$n_enriched)[1]/dim(M)[1])/(27191-26523+1))-((sum(M$n_depleted)[1]/dim(M)[1])/(27191-26523+1)), ((sum(ORF6$n_enriched)[1]/dim(ORF6)[1])/(27387-27202+1))-((sum(ORF6$n_depleted)[1]/dim(ORF6)[1])/(27387-27202+1)),
((sum(ORF7a$n_enriched)[1]/dim(ORF7a)[1])/(27759-27394+1))-((sum(ORF7a$n_depleted)[1]/dim(ORF7a)[1])/(27759-27394+1)), ((sum(ORF7b$n_enriched)[1]/dim(ORF7b)[1])/(27887-27756+1))-((sum(ORF7b$n_depleted)[1]/dim(ORF7b)[1])/(27887-27756+1)), ((sum(ORF8$n_enriched)[1]/dim(ORF8)[1])/(28259-27894+1))-((sum(ORF8$n_depleted)[1]/dim(ORF8)[1])/(28259-27894+1)),
((sum(N$n_enriched)[1]/dim(N)[1])/(29533-28274+1))-((sum(N$n_depleted)[1]/dim(N)[1])/(29533-28274+1)),
((sum(ORF10$n_enriched)[1]/dim(ORF10)[1])/(29674-29558+1))-((sum(ORF10$n_depleted)[1]/dim(ORF10)[1])/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

## The humans exclusive group
calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched <- function(df) {

ORF1ab <- df %>% dplyr::filter(gene == "ORF1ab")
S <- df %>% dplyr::filter(gene == "S")
ORF3a <- df %>% dplyr::filter(gene == "ORF3a")
E <- df %>% dplyr::filter(gene == "E")
M <- df %>% dplyr::filter(gene == "M")
ORF6 <- df %>% dplyr::filter(gene == "ORF6")
ORF7a <- df %>% dplyr::filter(gene == "ORF7a")
ORF7b <- df %>% dplyr::filter(gene == "ORF7b")
ORF8 <- df %>% dplyr::filter(gene == "ORF8")
N <- df %>% dplyr::filter(gene == "N")
ORF10 <- df %>% dplyr::filter(gene == "ORF10")

zamowienie_corona_gene <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum(ORF1ab$n_depleted)[1]/dim(ORF1ab)[1])/(21555-266+1)), ((sum(S$n_depleted)[1]/dim(S)[1])/(25384-21563+1)), ((sum(ORF3a$n_depleted)[1]/dim(ORF3a)[1])/(26220-25393+1)), ((sum(E$n_depleted)[1]/dim(E)[1])/(26472-26245+1)), ((sum(M$n_depleted)[1]/dim(M)[1])/(27191-26523+1)), ((sum(ORF6$n_depleted)[1]/dim(ORF6)[1])/(27387-27202+1)),
((sum(ORF7a$n_depleted)[1]/dim(ORF7a)[1])/(27759-27394+1)), ((sum(ORF7b$n_depleted)[1]/dim(ORF7b)[1])/(27887-27756+1)), ((sum(ORF8$n_depleted)[1]/dim(ORF8)[1])/(28259-27894+1)),
((sum(N$n_depleted)[1]/dim(N)[1])/(29533-28274+1)),
((sum(ORF10$n_depleted)[1]/dim(ORF10)[1])/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

## The humans favourable group
calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable <- function(df) {

ORF1ab <- df %>% dplyr::filter(gene == "ORF1ab")
S <- df %>% dplyr::filter(gene == "S")
ORF3a <- df %>% dplyr::filter(gene == "ORF3a")
E <- df %>% dplyr::filter(gene == "E")
M <- df %>% dplyr::filter(gene == "M")
ORF6 <- df %>% dplyr::filter(gene == "ORF6")
ORF7a <- df %>% dplyr::filter(gene == "ORF7a")
ORF7b <- df %>% dplyr::filter(gene == "ORF7b")
ORF8 <- df %>% dplyr::filter(gene == "ORF8")
N <- df %>% dplyr::filter(gene == "N")
ORF10 <- df %>% dplyr::filter(gene == "ORF10")

zamowienie_corona_gene <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum(ORF1ab$n_depleted)[1]/dim(ORF1ab)[1])/(21555-266+1)-((sum(ORF1ab$n_enriched)[1]/dim(ORF1ab)[1])/(21555-266+1))), ((sum(S$n_depletedd)[1]/dim(S)[1])/(25384-21563+1))-((sum(S$n_enriched)[1]/dim(S)[1])/(25384-21563+1)), ((sum(ORF3a$n_depleted)[1]/dim(ORF3a)[1])/(26220-25393+1))-((sum(ORF3a$n_enriched)[1]/dim(ORF3a)[1])/(26220-25393+1)), ((sum(E$n_depleted)[1]/dim(E)[1])/(26472-26245+1))-((sum(E$n_enriched)[1]/dim(E)[1])/(26472-26245+1)), ((sum(M$n_depleted)[1]/dim(M)[1])/(27191-26523+1))-((sum(M$n_enriched)[1]/dim(M)[1])/(27191-26523+1)), ((sum(ORF6$n_depleted)[1]/dim(ORF6)[1])/(27387-27202+1))-((sum(ORF6$n_enriched)[1]/dim(ORF6)[1])/(27387-27202+1)),
((sum(ORF7a$n_depleted)[1]/dim(ORF7a)[1])/(27759-27394+1))-((sum(ORF7a$n_enriched)[1]/dim(ORF7a)[1])/(27759-27394+1)), ((sum(ORF7b$n_depleted)[1]/dim(ORF7b)[1])/(27887-27756+1))-((sum(ORF7b$n_enriched)[1]/dim(ORF7b)[1])/(27887-27756+1)), ((sum(ORF8$n_depleted)[1]/dim(ORF8)[1])/(28259-27894+1))-((sum(ORF8$n_enriched)[1]/dim(ORF8)[1])/(28259-27894+1)),
((sum(N$n_depleted)[1]/dim(N)[1])/(29533-28274+1))-((sum(N$n_enriched)[1]/dim(N)[1])/(29533-28274+1)),
((sum(ORF10$n_depleted)[1]/dim(ORF10)[1])/(29674-29558+1))-((sum(ORF10$n_enriched)[1]/dim(ORF10)[1])/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}
