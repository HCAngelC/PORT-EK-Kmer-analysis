# deers
subsets_d <- c(rep("common", 4), rep("enriched", 2), rep("unique", 2))
count_d <- c(7154, 8993, 10379, 15420, 146, 1324, 548, 458)
cat_d <- c("deer over-represented", "humans over-represented", "single time period", "no significance", "deer enriched", "humans enriched", "deer unique", "human unique")

df_count_order_d <- data.frame(subsets_d, count_d, cat_d)

#ggplot(df_count_order_d, aes(area = "count_d", fill = "cat_d", label = "subsets_d", subgroup = "cat_d"))+geom_treemap(stat = "identity")

treemap(dtf = df_count_order_d,
        index = c("cat_d", "subsets_d"),
        vSize = "count_d",
        vColor = "cat_d")

# bats

subsets_b <- c(rep("common", 2), rep("enriched", 2), rep("unique", 2))
count_b <- c(2744, 29822, 2984, 29933, 29, 2676)
cat_b <- rep("bat common", "human over-represented: common", "bat enriched", "human over-represented: enriched", "bat unique", "human unique") 

df_count_order_b <- data.frame(subsets_b, count_b, cat_b)

treemap(dtf = df_count_order_b,
        index = c("cat_b", "subsets_b"),
        vSize = "count_b",
        vColor = "cat_b")
