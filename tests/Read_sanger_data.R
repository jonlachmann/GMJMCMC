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
x4 = t(read.table(paste("YRI_parents_norm_march2007.txt",sep=""),header=T))

x = data.frame(y=rep(NA,nrow(df)),x=matrix(NA,nrow=nrow(df),ncol=ncol(df)-1))
x = matrix(0,nrow=nrow(df),ncol=ncol(df))
rownames(x) = rownames(df)
ind = pmatch(rownames(x1)[-1],rownames(df))
n = 24266
#n = ncol(df)
x[ind,-1] = as.numeric(as.matrix(x1[-1,-n]))
x[ind,1] = as.numeric(x1[-1,n])
ind = pmatch(rownames(x2)[-1],rownames(df))
x[ind,-1] = as.numeric(x2[-1,-n])
x[ind,1] = as.numeric(x2[-1,n])
ind = pmatch(rownames(x3)[-1],rownames(df))
x[ind,-1] = as.numeric(x3[-1,-n])
x[ind,1] = as.numeric(x3[-1,n])
ind = pmatch(rownames(x4)[-1],rownames(df))
x[ind,-1] = as.numeric(x4[-1,-n])
x[ind,1] = as.numeric(x4[-1,n])

show(range(df-x))
