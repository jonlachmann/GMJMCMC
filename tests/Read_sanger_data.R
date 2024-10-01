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

nam = x1[1,]
x1 = apply(x1[-1,],2,as.numeric)
x2 = apply(x2[-1,],2,as.numeric)
x3 = apply(x3[-1,],2,as.numeric)
x4 = apply(x4[-1,],2,as.numeric)
x = rbind(x1,x2,x3,x4)
colnames(x) = nam

#Choose response variable
n = 24266
x = cbind(x[,n],x[,-n])

SangerData = x
