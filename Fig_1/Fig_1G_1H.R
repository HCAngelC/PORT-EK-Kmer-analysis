# Categorize enriched k-mers in deer in different intrinsic properties
specific_15mers_overlaid_coronavirus <- read.table("deer_specific_15mers_over_corona.txt", stringsAsFactors = F, header = F) %>% dplyr::select(V2, V3, V4, V5, V6, V7, V8, V10, V11, V12) %>% dplyr::rename(kmer_start = V2, kmer_end = V3, n_enriched = V4, n_depleted = V5, group_JW = V6, gene_JW = V7, series_id = V8, gene_start = V10, gene_end = V11, gene = V12)

#protein coding; all k-mers
sp_15_pc_enrich_wszystko <- specific_15mers_overlaid_coronavirus %>% dplyr::filter(n_enriched != 0 & n_depleted == 0) %>% dplyr::mutate(group = "deer-exclusive", kmer = "all", region = "pc")
sp_15_pc_deplete_wszystko <- specific_15mers_overlaid_coronavirus %>% dplyr::filter(n_enriched == 0 & n_depleted != 0) %>% dplyr::mutate(group = "human-exclusive", kmer = "all", region = "pc")

pc_pool_enrich_deplete <- dplyr::bind_rows(sp_15_pc_enrich_wszystko, sp_15_pc_deplete_wszystko)
unclarify_kmers_pc <- dplyr::anti_join(specific_15mers_overlaid_coronavirus, pc_pool_enrich_deplete, by = "series_id")

pc_15_fav_enrich <- unclarify_kmers_pc %>% dplyr::filter(n_enriched > n_depleted) %>% dplyr::mutate(group = "deer-favorable", kmer = "all", region = "pc")
pc_15_fav_deplete <- unclarify_kmers_pc %>% dplyr::filter(n_enriched < n_depleted) %>% dplyr::mutate(group = "human-favorable", kmer = "all", region = "pc")
pc_15_compariable <- unclarify_kmers_pc %>% dplyr::filter(n_enriched != 0 & n_depleted !=0) %>% dplyr::filter(n_enriched == n_depleted) %>% dplyr::mutate(group = "E", kmer = "all", region = "PC")

#Non-protein coding; all k-mers
sp_15_npc <- dplyr::anti_join(Specific_15mers_bed, specific_15mers_overlaid_coronavirus, by = "series_id") %>% dplyr::rename(kmer_start = pos, kmer_end = pos_1, group_JW = group) %>% dplyr::mutate(gene_start = 0, gene_end = 0) %>% dplyr::select(kmer_start, kmer_end, n_enriched, n_depleted, group_JW, series_id, gene_start, gene_end, gene)

sp_15_npc_enrich <- sp_15_npc %>% dplyr::filter(n_enriched != 0 & n_depleted == 0) %>% dplyr::mutate(group = "deer-exclusive", kmer = "all", region = "npc")
sp_15_npc_deplete <- sp_15_npc %>% dplyr::filter(n_enriched == 0 & n_depleted != 0) %>% dplyr::mutate(group = "human-exclusive", kmer = "all", region = "npc")

pool_enrich_deplete <- dplyr::bind_rows(sp_15_npc_enrich, sp_15_npc_deplete)
unclarify_kmers_npc <- dplyr::anti_join(sp_15_npc, pool_enrich_deplete, by = "series_id")

npc_15_fav_enrich <- unclarify_kmers_npc %>% dplyr::filter(n_enriched > n_depleted) %>% dplyr::mutate(group = "deer-favorable", kmer = "all", region = "npc")
npc_15_fav_deplete <- unclarify_kmers_npc %>% dplyr::filter(n_enriched < n_depleted) %>% dplyr::mutate(group = "human-favorable", kmer = "all", region = "npc")
npc_15_compariable <- unclarify_kmers_npc %>% dplyr::filter(n_enriched != 0 & n_depleted !=0)  %>% dplyr::filter(n_enriched == n_depleted) %>% dplyr::mutate(group = "dh-comparable", kmer = "all", region = "npc")

df_matrix_wszystko <- dplyr::bind_rows(sp_15_pc_deplete_wszystko, pc_15_fav_deplete, pc_15_compariable, sp_15_npc_enrich, sp_15_npc_deplete, npc_15_fav_enrich, npc_15_fav_deplete, npc_15_compariable)
df_matrix_pc_wszystko <- dplyr::bind_rows(sp_15_pc_deplete_wszystko, pc_15_fav_deplete, pc_15_compariable)
df_matrix_npc_wszystko <- dplyr::bind_rows(sp_15_npc_enrich, sp_15_npc_deplete, npc_15_fav_enrich, npc_15_fav_deplete, npc_15_compariable)

ratio <- c(dim(sp_15_pc_deplete_wszystko)[1]/dim(df_matrix_pc_wszystko)[1], dim(pc_15_fav_deplete)[1]/dim(df_matrix_pc_wszystko)[1], dim(pc_15_compariable)[1]/dim(df_matrix_pc_wszystko)[1], dim(sp_15_npc_enrich)[1]/dim(df_matrix_npc_wszystko)[1], dim(sp_15_npc_deplete)[1]/dim(df_matrix_npc_wszystko)[1], dim(npc_15_fav_enrich)[1]/dim(df_matrix_npc_wszystko)[1], dim(npc_15_fav_deplete)[1]/dim(df_matrix_npc_wszystko)[1], dim(npc_15_compariable)[1]/dim(df_matrix_npc_wszystko)[1])

group <- c("human-exclusive", "human-favorable", "dh-comparable", "deer-exclusive", "human-exclusive", "deer-favorable", "human-favorable", "dh-comparable")
region <- c(rep("pc", 3), rep("npc", 5))

df_stacked_bar <- data.frame(ratio, group, region)
Kierunek_region <- c("pc", "npc")
df_stacked_bar$region <- factor(df_stacked_bar$region, levels = Kierunek_region)

# Figure 1G; deer
> df_deer_stacked_bar
       ratio group region
1 0.88042948  "human-exclusive"  pc
2 0.01952172   "human-favorable"  pc
3 0.10004880   "dh-comparable"  pc
4 0.19337017   "deer-exclusive"  npc
5 0.61325967   "human-exclusive"  npc
6 0.03867403  "deer-favorable"  npc
7 0.07734807  "human-favorable"  npc
8 0.07734807  "dh-comparable"  npc

ggplot(df_deer_stacked_bar, aes(x = region, y = ratio, fill = group))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+xlab("Viral genomic regions")+ylab("Propertion")+ylim(0, 1.05)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

----
# Categorize enriched k-mers in bats in different intrinsic properties
specific_15mers_bat_overlaid_coronavirus <- read.table("bat_specific_15mers_over_corona.txt", stringsAsFactors = F, header = F) %>% dplyr::select(V2, V3, V4, V5, V6, V7, V8, V10, V11, V12) %>% dplyr::rename(kmer_start = V2, kmer_end = V3, n_enriched = V4, n_depleted = V5, group_JW = V6, gene_JW = V7, series_id = V8, gene_start = V10, gene_end = V11, gene = V12)

#protein coding; all k-mers
sp_15_pc_enrich_bat_wszystko <- specific_15mers_bat_overlaid_coronavirus %>% dplyr::filter(n_enriched != 0 & n_depleted == 0) %>% dplyr::mutate(group = "bat-exclusive", kmer = "all", region = "pc")
sp_15_pc_deplete_bat_wszystko <- specific_15mers_bat_overlaid_coronavirus %>% dplyr::filter(n_enriched == 0 & n_depleted != 0) %>% dplyr::mutate(group = "human-exclusive", kmer = "all", region = "pc")

pc_pool_bat_enrich_deplete <- dplyr::bind_rows(sp_15_pc_enrich_bat_wszystko, sp_15_pc_deplete_bat_wszystko)
unclarify_kmers_pc_bat <- dplyr::anti_join(specific_15mers_bat_overlaid_coronavirus, pc_pool_bat_enrich_deplete, by = "series_id")

pc_15_fav_enrich_bat <- unclarify_kmers_pc_bat %>% dplyr::filter(n_enriched > n_depleted) %>% dplyr::mutate(group = "bat-favorable", kmer = "all", region = "pc")
pc_15_fav_deplete_bat <- unclarify_kmers_pc_bat %>% dplyr::filter(n_enriched < n_depleted) %>% dplyr::mutate(group = "human-favorable", kmer = "all", region = "pc")
pc_15_compariable_bat <- unclarify_kmers_pc_bat %>% dplyr::filter(n_enriched != 0 & n_depleted !=0) %>% dplyr::filter(n_enriched == n_depleted) %>% dplyr::mutate(group = "E", kmer = "all", region = "pc")

#Non-protein coding; all k-mers
sp_15_npc_bat <- dplyr::anti_join(Specific_15mers_bat_bed, specific_15mers_bat_overlaid_coronavirus, by = "series_id") %>% dplyr::rename(kmer_start = pos, kmer_end = pos_1, group_JW = group) %>% dplyr::mutate(gene_start = 0, gene_end = 0) %>% dplyr::select(kmer_start, kmer_end, n_enriched, n_depleted, group_JW, series_id, gene_start, gene_end, gene)

sp_15_npc_enrich_bat <- sp_15_npc_bat %>% dplyr::filter(n_enriched != 0 & n_depleted == 0) %>% dplyr::mutate(group = "bat-exclusive", kmer = "all", region = "npc")
sp_15_npc_deplete_bat <- sp_15_npc_bat %>% dplyr::filter(n_enriched == 0 & n_depleted != 0) %>% dplyr::mutate(group = "human-exclusive", kmer = "all", region = "npc")

pool_enrich_deplete_bat <- dplyr::bind_rows(sp_15_npc_enrich_bat, sp_15_npc_deplete_bat)
unclarify_kmers_npc_bat <- dplyr::anti_join(sp_15_npc_bat, pool_enrich_deplete_bat, by = "series_id")

npc_15_fav_enrich_bat <- unclarify_kmers_npc_bat %>% dplyr::filter(n_enriched > n_depleted) %>% dplyr::mutate(group = "bat-favorable", kmer = "all", region = "npc")
npc_15_fav_deplete_bat <- unclarify_kmers_npc_bat %>% dplyr::filter(n_enriched < n_depleted) %>% dplyr::mutate(group = "human-favorable", kmer = "all", region = "npc")
npc_15_compariable_bat <- unclarify_kmers_npc_bat %>% dplyr::filter(n_enriched != 0 & n_depleted !=0)  %>% dplyr::filter(n_enriched == n_depleted) %>% dplyr::mutate(group = "dh-comparable", kmer = "all", region = "npc")

df_matrix_wszystko_bat <- dplyr::bind_rows(sp_15_pc_enrich_bat_wszystko, sp_15_pc_deplete_bat_wszystko, pc_15_fav_enrich_bat, pc_15_fav_deplete_bat, pc_15_compariable_bat, sp_15_npc_deplete_bat, npc_15_fav_deplete_bat, npc_15_compariable_bat)
df_matrix_pc_bat_wszystko <- dplyr::bind_rows(sp_15_pc_enrich_bat_wszystko, sp_15_pc_deplete_bat_wszystko, pc_15_fav_enrich_bat, pc_15_fav_deplete_bat, pc_15_compariable_bat)
df_matrix_npc_bat_wszystko <- dplyr::bind_rows(sp_15_npc_deplete_bat, npc_15_fav_deplete_bat, npc_15_compariable_bat)

ratio_bat <- c(dim(sp_15_pc_enrich_bat_wszystko)[1]/dim(df_matrix_pc_bat_wszystko)[1], dim(sp_15_pc_deplete_bat_wszystko)[1]/dim(df_matrix_pc_bat_wszystko)[1], dim(pc_15_fav_enrich_bat)[1]/dim(df_matrix_pc_bat_wszystko)[1], dim(pc_15_fav_deplete_bat)[1]/dim(df_matrix_pc_bat_wszystko)[1], dim(pc_15_compariable_bat)[1]/dim(df_matrix_pc_bat_wszystko)[1],  dim(sp_15_npc_deplete_bat)[1]/dim(df_matrix_npc_bat_wszystko)[1], dim(npc_15_fav_deplete_bat)[1]/dim(df_matrix_npc_bat_wszystko)[1], dim(npc_15_compariable_bat)[1]/dim(df_matrix_npc_bat_wszystko)[1])

group_bat <-c("bat-exclusive", "human-exclusive", "bat-favorable", "human-favorable", "dh-comparable", "human-exclusive", "human-favorable", "dh-comparable")
region_bat <- c(rep("pc", 5), rep("npc", 3))

df_bat_stacked_bar <- data.frame(ratio_bat, group_bat, region_bat)
Kierunek_region <- c("pc", "npc")
df_bat_stacked_bar$region_bat <- factor(df_bat_stacked_bar$region, levels = Kierunek_region)

# Figure 1H; bats
> df_bat_stacked_bar
  ratio_bat group_bat region_bat
1 0.0019775649   "bat-exclusive"  pc
2 0.7957993795   "human-exclusive"  pc
3 0.0006819189   "bat-favorable" pc
4 0.1969040881   "human-favorable"  pc
5 0.0046370487   "dh-comparable"  pc
6 0.6263498920   "human-exclusive"  npc
7 0.3714902808   "human-favorable"  npc
8 0.0021598272   "dh-comparable" npc

ggplot(df_bat_stacked_bar, aes(x = region_bat, y = ratio_bat, fill = group_bat))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+xlab("Viral genomic regions")+ylab("Propertion")+ylim(0, 1.05)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
