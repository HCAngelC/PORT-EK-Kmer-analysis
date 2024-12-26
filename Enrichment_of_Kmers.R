#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2024.12.25
##############################################################
# Input: A full list of k-mers associated with corresponding loci
# Object: A matrix of the enrichment of k-mer in corresponding SARS-CoV-2 genes
##############################################################

# R functions

## The deers exclusive group
calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_enriched_bd <- function(df) {

nsp1 <- df %>% dplyr::filter(gene == "nsp1")
nsp2 <- df %>% dplyr::filter(gene == "nsp2")
nsp3 <- df %>% dplyr::filter(gene == "nsp3")
nsp4 <- df %>% dplyr::filter(gene == "nsp4")
nsp5 <- df %>% dplyr::filter(gene == "nsp5")
nsp6 <- df %>% dplyr::filter(gene == "nsp6")
nsp7 <- df %>% dplyr::filter(gene == "nsp7")
nsp8 <- df %>% dplyr::filter(gene == "nsp8")
nsp9 <- df %>% dplyr::filter(gene == "nsp9")
nsp10 <- df %>% dplyr::filter(gene == "nsp10")
nsp11 <- df %>% dplyr::filter(gene == "nsp11")
nsp12 <- df %>% dplyr::filter(gene == "nsp12")
nsp13 <- df %>% dplyr::filter(gene == "nsp13")
nsp14 <- df %>% dplyr::filter(gene == "nsp14")
nsp15 <- df %>% dplyr::filter(gene == "nsp15")
nsp16 <- df %>% dplyr::filter(gene == "nsp16")
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

zamowienie_corona_gene <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum((nsp1$n_enriched)[1]/dim(nsp1)[1]))/(805-266+1)), (sum((nsp2$n_enriched)[1]/dim(nsp2)[1]))/(2719-806+1), (sum((nsp3$n_enriched)[1]/dim(nsp3)[1]))/(8554-2720+1), (sum((nsp4$n_enriched)[1]/dim(nsp4)[1]))/(10054-8555+1), (sum((nsp5$n_enriched)[1]/dim(nsp5)[1]))/(10972-10055+1), (sum((nsp6$n_enriched)[1]/dim(nsp6)[1]))/(11842-10973+1), (sum((nsp7$n_enriched)[1]/dim(nsp7)[1]))/(12091-11843+1), (sum((nsp8$n_enriched)[1]/dim(nsp8)[1]))/(12685-12092+1), (sum((nsp9$n_enriched)[1]/dim(nsp9)[1]))/(13024-12686+1), (sum((nsp10$n_enriched)[1]/dim(nsp10)[1]))/(13441-13025+1), (sum((nsp11$n_enriched)[1]/dim(nsp11)[1]))/(13480-13442+1), (sum((nsp12$n_enriched)[1]/dim(nsp12)[1]))/(16236-13481+1), (sum((nsp13$n_enriched)[1]/dim(nsp13)[1]))/(18039-16237+1), (sum((nsp14$n_enriched)[1]/dim(nsp14)[1]))/(19620-18040+1), (sum((nsp15$n_enriched)[1]/dim(nsp15)[1]))/(20658-19621+1), (sum((nsp16$n_enriched)[1]/dim(nsp16)[1]))/(21552-20659+1), (sum((S$n_enriched)[1]/dim(S)[1]))/(25384-21563+1), (sum((ORF3a$n_enriched)[1]/dim(ORF3a)[1]))/(26220-25393+1), (sum((E$n_enriched)[1]/dim(E)[1]))/(26472-26245+1), (sum((M$n_enriched)[1]/dim(M)[1]))/(27191-26523+1), (sum((ORF6$n_enriched)[1]/dim(ORF6)[1]))/(27387-27202+1),
(sum((ORF7a$n_enriched)[1]/dim(ORF7a)[1]))/(27759-27394+1), (sum((ORF7b$n_enriched)[1]/dim(ORF7b)[1]))/(27887-27756+1), (sum((ORF8$n_enriched)[1]/dim(ORF8)[1]))/(28259-27894+1),
(sum((N$n_enriched)[1]/dim(N)[1]))/(29533-28274+1),
(sum((ORF10$n_enriched)[1]/dim(ORF10)[1]))/(29674-29558+1))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_animal_favourable_bd <- function(df) {

nsp1 <- df %>% dplyr::filter(gene == "nsp1")
nsp2 <- df %>% dplyr::filter(gene == "nsp2")
nsp3 <- df %>% dplyr::filter(gene == "nsp3")
nsp4 <- df %>% dplyr::filter(gene == "nsp4")
nsp5 <- df %>% dplyr::filter(gene == "nsp5")
nsp6 <- df %>% dplyr::filter(gene == "nsp6")
nsp7 <- df %>% dplyr::filter(gene == "nsp7")
nsp8 <- df %>% dplyr::filter(gene == "nsp8")
nsp9 <- df %>% dplyr::filter(gene == "nsp9")
nsp10 <- df %>% dplyr::filter(gene == "nsp10")
nsp11 <- df %>% dplyr::filter(gene == "nsp11")
nsp12 <- df %>% dplyr::filter(gene == "nsp12")
nsp13 <- df %>% dplyr::filter(gene == "nsp13")
nsp14 <- df %>% dplyr::filter(gene == "nsp14")
nsp15 <- df %>% dplyr::filter(gene == "nsp15")
nsp16 <- df %>% dplyr::filter(gene == "nsp16")
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

zamowienie_corona_gene <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c((((sum(nsp1$n_enriched)[1]/dim(nsp1)[1]))/(805-266+1))-(((sum(nsp1$n_depleted)[1]/dim(nsp1)[1]))/(805-266+1)),
(((sum(nsp2$n_enriched)[1]/dim(nsp2)[1]))/(2719-806+1))-(((sum(nsp2$n_depleted)[1]/dim(nsp2)[1]))/(2719-806+1)),
(((sum(nsp3$n_enriched)[1]/dim(nsp3)[1]))/(8554-2720+1))-(((sum(nsp3$n_depleted)[1]/dim(nsp3)[1]))/(8554-2720+1)),
(((sum(nsp4$n_enriched)[1]/dim(nsp4)[1]))/(10054-8555+1))-(((sum(nsp4$n_depleted)[1]/dim(nsp4)[1]))/(10054-8555+1)),
(((sum(nsp5$n_enriched)[1]/dim(nsp5)[1]))/(10972-10055+1))-(((sum(nsp5$n_depleted)[1]/dim(nsp5)[1]))/(10972-10055+1)),
(((sum(nsp6$n_enriched)[1]/dim(nsp6)[1]))/(11842-10973+1))-(((sum(nsp6$n_depleted)[1]/dim(nsp6)[1]))/(11842-10973+1)),
(((sum(nsp7$n_enriched)[1]/dim(nsp7)[1]))/(12091-11843+1))-(((sum(nsp7$n_depleted)[1]/dim(nsp7)[1]))/(12091-11843+1)),
(((sum(nsp8$n_enriched)[1]/dim(nsp8)[1]))/(12685-12092+1))-(((sum(nsp8$n_depleted)[1]/dim(nsp8)[1]))/(12685-12092+1)),
(((sum(nsp9$n_enriched)[1]/dim(nsp9)[1]))/(13024-12686+1))-(((sum(nsp9$n_depleted)[1]/dim(nsp9)[1]))/(13024-12686+1)),
(((sum(nsp10$n_enriched)[1]/dim(nsp10)[1]))/(13441-13025+1))-(((sum(nsp10$n_depleted)[1]/dim(nsp10)[1]))/(13441-13025+1)),
(((sum(nsp11$n_enriched)[1]/dim(nsp11)[1]))/(13480-13442+1))-(((sum(nsp11$n_depleted)[1]/dim(nsp11)[1]))/(13480-13442+1)),
(((sum(nsp12$n_enriched)[1]/dim(nsp12)[1]))/(16236-13481+1))-(((sum(nsp12$n_depleted)[1]/dim(nsp12)[1]))/(16236-13481+1)),
(((sum(nsp13$n_enriched)[1]/dim(nsp13)[1]))/(18039-16237+1))-(((sum(nsp13$n_depleted)[1]/dim(nsp13)[1]))/(18039-16237+1)),
(((sum(nsp14$n_enriched)[1]/dim(nsp14)[1]))/(19620-18040+1))-(((sum(nsp14$n_depleted)[1]/dim(nsp14)[1]))/(19620-18040+1)),
(((sum(nsp15$n_enriched)[1]/dim(nsp15)[1]))/(20658-19621+1))-(((sum(nsp15$n_depleted)[1]/dim(nsp15)[1]))/(20658-19621+1)),
(((sum(nsp16$n_enriched)[1]/dim(nsp16)[1]))/(21552-20659+1))-(((sum(nsp16$n_depleted)[1]/dim(nsp16)[1]))/(21552-20659+1)), (((sum(S$n_enriched)[1]/dim(S)[1]))/(25384-21563+1))-(((sum(S$n_depleted)[1]/dim(S)[1]))/(25384-21563+1)), (((sum(ORF3a$n_enriched)[1]/dim(ORF3a)[1]))/(26220-25393+1))-(((sum(ORF3a$n_depleted)[1]/dim(ORF3a)[1]))/(26220-25393+1)), (((sum(E$n_enriched)[1]/dim(E)[1]))/(26472-26245+1))-(((sum(E$n_depleted)[1]/dim(E)[1]))/(26472-26245+1)), (((sum(M$n_enriched)[1]/dim(M)[1]))/(27191-26523+1))-(((sum(M$n_depleted)[1]/dim(M)[1]))/(27191-26523+1)), (((sum(ORF6$n_enriched)[1]/dim(ORF6)[1]))/(27387-27202+1))-(((sum(ORF6$n_depleted)[1]/dim(ORF6)[1]))/(27387-27202+1)),
(((sum(ORF7a$n_enriched)[1]/dim(ORF7a)[1]))/(27759-27394+1))-(((sum(ORF7a$n_depleted)[1]/dim(ORF7a)[1]))/(27759-27394+1)), (((sum(ORF7b$n_enriched)[1]/dim(ORF7b)[1]))/(27887-27756+1))-(((sum(ORF7b$n_depleted)[1]/dim(ORF7b)[1]))/(27887-27756+1)), (((sum(ORF8$n_enriched)[1]/dim(ORF8)[1]))/(28259-27894+1))-(((sum(ORF8$n_depleted)[1]/dim(ORF8)[1]))/(28259-27894+1)),
(((sum(N$n_enriched)[1]/dim(N)[1]))/(29533-28274+1))-(((sum(N$n_depleted)[1]/dim(N)[1]))/(29533-28274+1)),
(((sum(ORF10$n_enriched)[1]/dim(ORF10)[1]))/(29674-29558+1))-(((sum(ORF10$n_depleted)[1]/dim(ORF10)[1]))/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_enriched_bd <- function(df) {

nsp1 <- df %>% dplyr::filter(gene == "nsp1")
nsp2 <- df %>% dplyr::filter(gene == "nsp2")
nsp3 <- df %>% dplyr::filter(gene == "nsp3")
nsp4 <- df %>% dplyr::filter(gene == "nsp4")
nsp5 <- df %>% dplyr::filter(gene == "nsp5")
nsp6 <- df %>% dplyr::filter(gene == "nsp6")
nsp7 <- df %>% dplyr::filter(gene == "nsp7")
nsp8 <- df %>% dplyr::filter(gene == "nsp8")
nsp9 <- df %>% dplyr::filter(gene == "nsp9")
nsp10 <- df %>% dplyr::filter(gene == "nsp10")
nsp11 <- df %>% dplyr::filter(gene == "nsp11")
nsp12 <- df %>% dplyr::filter(gene == "nsp12")
nsp13 <- df %>% dplyr::filter(gene == "nsp13")
nsp14 <- df %>% dplyr::filter(gene == "nsp14")
nsp15 <- df %>% dplyr::filter(gene == "nsp15")
nsp16 <- df %>% dplyr::filter(gene == "nsp16")
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

zamowienie_corona_gene <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c(((sum((nsp1$n_depleted)[1]/dim(nsp1)[1]))/(805-266+1)), (sum((nsp2$n_depleted)[1]/dim(nsp2)[1]))/(2719-806+1), (sum((nsp3$n_depleted)[1]/dim(nsp3)[1]))/(8554-2720+1), (sum((nsp4$n_depleted)[1]/dim(nsp4)[1]))/(10054-8555+1), (sum((nsp5$n_depleted)[1]/dim(nsp5)[1]))/(10972-10055+1), (sum((nsp6$n_depleted)[1]/dim(nsp6)[1]))/(11842-10973+1), (sum((nsp7$n_depleted)[1]/dim(nsp7)[1]))/(12091-11843+1), (sum((nsp8$n_depleted)[1]/dim(nsp8)[1]))/(12685-12092+1), (sum((nsp9$n_depleted)[1]/dim(nsp9)[1]))/(13024-12686+1), (sum((nsp10$n_depleted)[1]/dim(nsp10)[1]))/(13441-13025+1), (sum((nsp11$n_depleted)[1]/dim(nsp11)[1]))/(13480-13442+1), (sum((nsp12$n_depleted)[1]/dim(nsp12)[1]))/(16236-13481+1), (sum((nsp13$n_depleted)[1]/dim(nsp13)[1]))/(18039-16237+1), (sum((nsp14$n_depleted)[1]/dim(nsp14)[1]))/(19620-18040+1), (sum((nsp15$n_depleted)[1]/dim(nsp15)[1]))/(20658-19621+1), (sum((nsp16$n_depleted)[1]/dim(nsp16)[1]))/(21552-20659+1), (sum((S$n_depleted)[1]/dim(S)[1]))/(25384-21563+1), (sum((ORF3a$n_depleted)[1]/dim(ORF3a)[1]))/(26220-25393+1), (sum((E$n_depleted)[1]/dim(E)[1]))/(26472-26245+1), (sum((M$n_depleted)[1]/dim(M)[1]))/(27191-26523+1), (sum((ORF6$n_depleted)[1]/dim(ORF6)[1]))/(27387-27202+1),
(sum((ORF7a$n_depleted)[1]/dim(ORF7a)[1]))/(27759-27394+1), (sum((ORF7b$n_depleted)[1]/dim(ORF7b)[1]))/(27887-27756+1), (sum((ORF8$n_depleted)[1]/dim(ORF8)[1]))/(28259-27894+1),
(sum((N$n_depleted)[1]/dim(N)[1]))/(29533-28274+1),
(sum((ORF10$n_depleted)[1]/dim(ORF10)[1]))/(29674-29558+1))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}

calculate_proportion_kmers_over_coronavirus_genome_pc_wszystko_normalization_human_favourable_bd <- function(df) {

nsp1 <- df %>% dplyr::filter(gene == "nsp1")
nsp2 <- df %>% dplyr::filter(gene == "nsp2")
nsp3 <- df %>% dplyr::filter(gene == "nsp3")
nsp4 <- df %>% dplyr::filter(gene == "nsp4")
nsp5 <- df %>% dplyr::filter(gene == "nsp5")
nsp6 <- df %>% dplyr::filter(gene == "nsp6")
nsp7 <- df %>% dplyr::filter(gene == "nsp7")
nsp8 <- df %>% dplyr::filter(gene == "nsp8")
nsp9 <- df %>% dplyr::filter(gene == "nsp9")
nsp10 <- df %>% dplyr::filter(gene == "nsp10")
nsp11 <- df %>% dplyr::filter(gene == "nsp11")
nsp12 <- df %>% dplyr::filter(gene == "nsp12")
nsp13 <- df %>% dplyr::filter(gene == "nsp13")
nsp14 <- df %>% dplyr::filter(gene == "nsp14")
nsp15 <- df %>% dplyr::filter(gene == "nsp15")
nsp16 <- df %>% dplyr::filter(gene == "nsp16")
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

zamowienie_corona_gene <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
enrichment <- c((((sum(nsp1$n_depleted)[1]/dim(nsp1)[1]))/(805-266+1))-(((sum(nsp1$n_enriched)[1]/dim(nsp1)[1]))/(805-266+1)),
(((sum(nsp2$n_depleted)[1]/dim(nsp2)[1]))/(2719-806+1))-(((sum(nsp2$n_enriched)[1]/dim(nsp2)[1]))/(2719-806+1)),
(((sum(nsp3$n_depleted)[1]/dim(nsp3)[1]))/(8554-2720+1))-(((sum(nsp3$n_enriched)[1]/dim(nsp3)[1]))/(8554-2720+1)),
(((sum(nsp4$n_depleted)[1]/dim(nsp4)[1]))/(10054-8555+1))-(((sum(nsp4$n_enriched)[1]/dim(nsp4)[1]))/(10054-8555+1)),
(((sum(nsp5$n_depleted)[1]/dim(nsp5)[1]))/(10972-10055+1))-(((sum(nsp5$n_enriched)[1]/dim(nsp5)[1]))/(10972-10055+1)),
(((sum(nsp6$n_depleted)[1]/dim(nsp6)[1]))/(11842-10973+1))-(((sum(nsp6$n_enriched)[1]/dim(nsp6)[1]))/(11842-10973+1)),
(((sum(nsp7$n_depleted)[1]/dim(nsp7)[1]))/(12091-11843+1))-(((sum(nsp7$n_enriched)[1]/dim(nsp7)[1]))/(12091-11843+1)),
(((sum(nsp8$n_depleted)[1]/dim(nsp8)[1]))/(12685-12092+1))-(((sum(nsp8$n_enriched)[1]/dim(nsp8)[1]))/(12685-12092+1)),
(((sum(nsp9$n_depleted)[1]/dim(nsp9)[1]))/(13024-12686+1))-(((sum(nsp9$n_enriched)[1]/dim(nsp9)[1]))/(13024-12686+1)),
(((sum(nsp10$n_depleted)[1]/dim(nsp10)[1]))/(13441-13025+1))-(((sum(nsp10$n_enriched)[1]/dim(nsp10)[1]))/(13441-13025+1)),
(((sum(nsp11$n_depleted)[1]/dim(nsp11)[1]))/(13480-13442+1))-(((sum(nsp11$n_enriched)[1]/dim(nsp11)[1]))/(13480-13442+1)),
(((sum(nsp12$n_depleted)[1]/dim(nsp12)[1]))/(16236-13481+1))-(((sum(nsp12$n_enriched)[1]/dim(nsp12)[1]))/(16236-13481+1)),
(((sum(nsp13$n_depleted)[1]/dim(nsp13)[1]))/(18039-16237+1))-(((sum(nsp13$n_enriched)[1]/dim(nsp13)[1]))/(18039-16237+1)),
(((sum(nsp14$n_depleted)[1]/dim(nsp14)[1]))/(19620-18040+1))-(((sum(nsp14$n_enriched)[1]/dim(nsp14)[1]))/(19620-18040+1)),
(((sum(nsp15$n_depleted)[1]/dim(nsp15)[1]))/(20658-19621+1))-(((sum(nsp15$n_enriched)[1]/dim(nsp15)[1]))/(20658-19621+1)),
(((sum(nsp16$n_depleted)[1]/dim(nsp16)[1]))/(21552-20659+1))-(((sum(nsp16$n_enriched)[1]/dim(nsp16)[1]))/(21552-20659+1)), (((sum(S$n_depleted)[1]/dim(S)[1]))/(25384-21563+1))-(((sum(S$n_enriched)[1]/dim(S)[1]))/(25384-21563+1)), (((sum(ORF3a$n_depleted)[1]/dim(ORF3a)[1]))/(26220-25393+1))-(((sum(ORF3a$n_enriched)[1]/dim(ORF3a)[1]))/(26220-25393+1)), (((sum(E$n_depleted)[1]/dim(E)[1]))/(26472-26245+1))-(((sum(E$n_enriched)[1]/dim(E)[1]))/(26472-26245+1)), (((sum(M$n_depleted)[1]/dim(M)[1]))/(27191-26523+1))-(((sum(M$n_enriched)[1]/dim(M)[1]))/(27191-26523+1)), (((sum(ORF6$n_depleted)[1]/dim(ORF6)[1]))/(27387-27202+1))-(((sum(ORF6$n_enriched)[1]/dim(ORF6)[1]))/(27387-27202+1)),
(((sum(ORF7a$n_depleted)[1]/dim(ORF7a)[1]))/(27759-27394+1))-(((sum(ORF7a$n_enriched)[1]/dim(ORF7a)[1]))/(27759-27394+1)), (((sum(ORF7b$n_depleted)[1]/dim(ORF7b)[1]))/(27887-27756+1))-(((sum(ORF7b$n_enriched)[1]/dim(ORF7b)[1]))/(27887-27756+1)), (((sum(ORF8$n_depleted)[1]/dim(ORF8)[1]))/(28259-27894+1))-(((sum(ORF8$n_enriched)[1]/dim(ORF8)[1]))/(28259-27894+1)),
(((sum(N$n_depleted)[1]/dim(N)[1]))/(29533-28274+1))-(((sum(N$n_enriched)[1]/dim(N)[1]))/(29533-28274+1)),
(((sum(ORF10$n_depleted)[1]/dim(ORF10)[1]))/(29674-29558+1))-(((sum(ORF10$n_enriched)[1]/dim(ORF10)[1]))/(29674-29558+1)))

df_pool <- data.frame(zamowienie_corona_gene, enrichment)
df_pool$zamowienie_corona_gene <- factor(df_pool$zamowienie_corona_gene, levels = zamowienie_corona_gene)

return(df_pool)
}
