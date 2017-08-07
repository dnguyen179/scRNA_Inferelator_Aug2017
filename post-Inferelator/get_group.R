## get list of TFs according to predictor groups 

load("~/working_path/params_and_input.RData")


group_1 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.1
group_2 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.2
group_3 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.3
group_4 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.4
group_5 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.5
group_6 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.6
group_7 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.7
# group_8 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.8
# group_9 <- IN$grouped.pred$`frac_tp_100_perm_1--frac_fp_0_perm_1`$pred.groups$pred.group.9

group_1 <- t(matrix(c(group_1)))
group_2 <- t(matrix(c(group_2)))
group_3 <- t(matrix(c(group_3)))
group_4 <- t(matrix(c(group_4)))
group_5 <- t(matrix(c(group_5)))
group_6 <- t(matrix(c(group_6)))
group_7 <- t(matrix(c(group_7)))
# group_8 <- t(matrix(c(group_8)))
# group_9 <- t(matrix(c(group_9)))


write.table(group_1, file = "pred_group_1.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_2, file = "pred_group_2.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_3, file = "pred_group_3.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_4, file = "pred_group_4.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_5, file = "pred_group_5.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_6, file = "pred_group_6.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(group_7, file = "pred_group_7.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(group_8, file = "pred_group_8.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(group_9, file = "pred_group_9.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
