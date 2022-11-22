
Nout = 121

wld_1 = read.table("./OUTPUT/Check_seed/wld_sd_1111.txt", sep="\t")[,1:Nout]
mnt_1 = read.table("./OUTPUT/Check_seed/mnt_sd_1111.txt", sep="\t")[,1:Nout]
wld_2 = read.table("./OUTPUT/Check_seed/wld_sd_1212.txt", sep="\t")[,1:Nout]
mnt_2 = read.table("./OUTPUT/Check_seed/mnt_sd_1212.txt", sep="\t")[,1:Nout]


cn_1 = wld_1+mnt_1
ml_1 = mnt_1 / cn_1
cn_2 = wld_2+mnt_2
ml_2 = mnt_2 / cn_2

dim(cn_1)
dim(cn_2)

sum(abs(cn_1-cn_2))
sum(ml_1-ml_2)
range(ml_1)
range(ml_2)

