#load dataset on data_df
data_df = read.table("maize_bivariate_pvals.txt", header= TRUE, sep="")

#filter dataset using cutoff points
#pvalsB and pvalsM are column names (for pvalues) of the data set for this test(to be changed to column names of the loaded data set)
p_vals = data_df  %>% filter(pvalsB >=c1, pvalsM>=c2)

# convert selected pvalues to z values
z_val = as.data.frame(qnorm(as.matrix(p_vals), lower.tail = TRUE))
colnames(z_val) = c("zvals1", "zvals2")

# convert lambda(truncation points) to z values
z_val_extremums = as.data.frame(qnorm(as.matrix(cbind(c(c1,1),c(c2,1))), lower.tail = TRUE))

#selecting the lower bounds (minimum z values) of the bivariate z values
min_z1 <- z_val_extremums[1,1]
min_z2 <- z_val_extremums[1,2]