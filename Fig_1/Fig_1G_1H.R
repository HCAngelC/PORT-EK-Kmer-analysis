# Figure 1G; deer
> df_deer_stacked_bar
       ratio group region
1 0.88042948  human-exclusive  pc
2 0.01952172  human-favorable  pc
3 0.10004880  dh-comparable  pc
4 0.19337017  deer-exclusive  npc
5 0.61325967  human-exclusive  npc
6 0.03867403 deer-favorable  npc
7 0.07734807 human-favorable  npc
8 0.07734807 dh-comparable  npc

ggplot(df_deer_stacked_bar, aes(x = region, y = ratio, fill = group))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+xlab("Viral genomic regions")+ylab("Propertion")+ylim(0, 1.05)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# Figure 1H; bats
> df_bat_stacked_bar
  ratio_bat group_bat region_bat
1 0.0019775649  bat-exclusive  pc
2 0.7957993795  human-exclusive  pc
3 0.0006819189  bat-favorable pc
4 0.1969040881  human-favorable  pc
5 0.0046370487  dh-comparable  pc
6 0.6263498920  human-exclusive  npc
7 0.3714902808  human-favorable  npc
8 0.0021598272  dh-comparable npc

ggplot(df_bat_stacked_bar, aes(x = region_bat, y = ratio_bat, fill = group_bat))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = c("#ff33cc", "#ff99cc","#3366ff", "#99ccff", "#999999"))+xlab("Viral genomic regions")+ylab("Propertion")+ylim(0, 1.05)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
