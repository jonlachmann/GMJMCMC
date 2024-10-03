#Read Sanger data
#Download first from https://ftp.sanger.ac.uk/pub/genevar/ the files
 #https://ftp.sanger.ac.uk/pub/genevar/CEU_children_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/CHB_unrelated_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/JPT_unrelated_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/YRI_parents_norm_march2007.zip
#and unzip these

#Specify right path
path = "/mn/sarpanitu/ansatte-u2/geirs/prj/FBMS/data/"
x1 = t(read.table(paste(path,"CEU_parents_norm_march2007.txt",sep=""),header=T))
x2 = t(read.table(paste(path,"CHB_unrelated_norm_march2007.txt",sep=""),header=T))
x3 = t(read.table(paste(path,"JPT_unrelated_norm_march2007.txt",sep=""),header=T))
x4 = t(read.table(paste(path,"YRI_parents_norm_march2007.txt",sep=""),header=T))

rnames = c(rownames(x1)[-1],rownames(x2)[-1],rownames(x3)[-1],rownames(x4)[-1])
nam = x1[1,]
x1 = apply(x1[-1,],2,as.numeric)
x2 = apply(x2[-1,],2,as.numeric)
x3 = apply(x3[-1,],2,as.numeric)
x4 = apply(x4[-1,],2,as.numeric)
x = rbind(x1,x2,x3,x4)
colnames(x) = nam
rownames(x) = rnames

#Choose response variable
SangerData = x
usethis::use_data(SangerData,overwrite=TRUE)

#Reduced dataset 
df <- SangerData[,c(24266,1:24265,24267:ncol(SangerData))]
#Rename columns
colnames(df) = c("y",paste0("x",1:47292))

# Candidates  based on marginal p values
p.vec = unlist(mclapply(2:47293, function(x)cor.test(df[,1],df[,x])$p.value))
ids = sort(order(p.vec)[1:100])

#Reduce data
df = df[,c(1,1+ids)]
SangerData2 = df
usethis::use_data(SangerData2)
