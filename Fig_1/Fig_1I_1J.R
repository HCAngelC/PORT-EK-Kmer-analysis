# Figure 1I; deer

Fig_1I <- read.table("Fig_1I.txt")

ggplot(Fig_1I, aes(x = kmer_start, xend = kmer_end, y = series_id_3, yend = series_id_3, color = group))+geom_segment(linewidth = 3)+theme_bw()+xlab("Wuhan-Hu-1 NC_045512.2 genome")+ylab("K-mers")+xlim(0, 30000)+scale_color_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 0, vjust = 0), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# Figure 1J; bats

Fig_1J <- read.table("Fig_1J.txt")

ggplot(Fig_1J, aes(x = kmer_start, xend = kmer_end, y = series_id_3, yend = series_id_3, color = group))+geom_segment(linewidth = 3)+theme_bw()+xlab("Wuhan-Hu-1 NC_045512.2 genome")+ylab("K-mers")+xlim(0, 30000)+scale_color_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 0, vjust = 0), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
