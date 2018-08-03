####
# This script is specific for cacao and for the use with Matina reference. 
# chromosome names are specified as "scaffold_#" where # corresponds to the given number
####

args<-commandArgs(TRUE) 
nchrom=args[[1]]
mydirectory=eval(args[[2]])

CHRLEN = c(38988864, 42436413, 34397752, 33492547, 40442896,
27290986, 24379470, 21543242, 42035188, 25448839);
data=c()
#reduced_len=c()
#reduced_len2=c()


for (i in 1:nchrom) {
	nam <- paste("chr", i,sep = "")
	assign(nam,read.table(paste(mydirectory,"/scaffold_",i,".tab",sep=""),header=TRUE))
	pos_chr <- rep(NA,dim(eval(parse(text=paste("chr",i,sep=""))))[1] + 2)
	pos_chr[1] <- 1
	pos_chr[length(pos_chr)] <- CHRLEN[i]
	for(j in 1:dim(eval(parse(text=paste("chr",i,sep=""))))[1]){pos_chr[j + 1] <- eval(parse(text=paste("chr",i,sep="")))[j,2] + (eval(parse(text=paste("chr",i,sep="")))[j,3] - eval(parse(text=paste("chr",i,sep="")))[j,2])/2}
	mid_pos_chr <- pos_chr[3:length(pos_chr) - 1]
	frag_len_chr <- rep(NA,length(mid_pos_chr))
	for(j in 1:length(mid_pos_chr)){frag_len_chr[j] <- pos_chr[j + 1] - pos_chr[j]}
	frag_len_chr[which(frag_len_chr == max(frag_len_chr))] <- 0
	chr_len <- cbind(mid_pos_chr,frag_len_chr)
	temp1 <- chr_len[which(chr_len[,2] >= 200 & chr_len[,2] <= 700),]
    temp2 <- chr_len[which(chr_len[,2] >= 200 & chr_len[,2] <= 1000),]
	lab_chr <- rep(i,dim(temp1)[1])
	lab_chr_2 <- rep(i,dim(temp2)[1])
	nam1 <- paste("reduced_chr",i,"_len",sep="")
	assign(nam1,cbind(lab_chr,temp1))
	nam2 <- paste("reduced_chr",i,"_len2",sep="")
	assign(nam2,cbind(lab_chr_2,temp2))
}
data <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10)

###
#To plot:
####


combined_200_700 <- rbind(reduced_chr1_len,reduced_chr2_len,reduced_chr3_len,reduced_chr4_len,reduced_chr5_len,reduced_chr6_len,reduced_chr7_len,reduced_chr8_len,reduced_chr9_len,reduced_chr10_len)
combined_200_1000 <- rbind(reduced_chr1_len2,reduced_chr2_len2,reduced_chr3_len2,reduced_chr4_len2,reduced_chr5_len2,reduced_chr6_len2,reduced_chr7_len2,reduced_chr8_len2,reduced_chr9_len2,reduced_chr10_len2)
Total_frags200_700 <- dim(combined_200_700)[1]
cat("total number of fragments in 200-700 bp range is",Total_frags200_700,"\n")
Total_frags200_1000 <- dim(combined_200_1000)[1]
cat("total number of fragments in 200-1000 bp range is",Total_frags200_1000,"\n")
Length_Total_frags200_700 <- sum(combined_200_700[,3])
cat("total length of captured genome with 200-700 bp fragments is",Length_Total_frags200_700/1000," KB\n")
Length_Total_frags200_1000 <- sum(combined_200_1000[,3])
cat("total length of captured genome with 200-1000 bp fragments is",Length_Total_frags200_1000/1000," KB\n")

pdf("distributionFragments_afterRestrict.pdf", width=10,height=6);
par(mfrow=c(1,2))
hist(combined_200_700[,3], main=paste("restriction analysis BsaXI (",Length_Total_frags200_700," bp)",sep=""),xlab="fragment sizes",ylab="number of fragments")
hist(combined_200_1000[,3], main=paste("restriction analysis BsaXI (",Length_Total_frags200_1000," bp)",sep=""),xlab="fragment sizes",ylab="number of fragments")
#dev.off()
graphics.off()

###
# to plot positions
###

jpeg("enzymes_cuts_chr_location.jpeg",height=6,width=8,pointsize=15,units="in",res=300);
plot(1,1,type='n',xlim=c(0,max(CHRLEN)/1e6),ylim=c(1,10),xlab="genome position (Mb)", ylab="chromosome", main="Location enzymes digestion", yaxt="n")
at <- seq(from = 1, to = 10, by = 1)
axis(side = 2, at = at, las=2)
for(c in 1:10){
       segments(x0=1,x1=CHRLEN[c]/1e6,y0=c,y1=c) #plot chromosome
 		}
for(i in 1:dim(reduced_chr1_len)[1]){
	points(reduced_chr1_len[i,2]/1e6,1, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr2_len)[1]){
	points(reduced_chr2_len[i,2]/1e6,2, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr3_len)[1]){
	points(reduced_chr3_len[i,2]/1e6,3, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr4_len)[1]){
	points(reduced_chr4_len[i,2]/1e6,4, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr5_len)[1]){
	points(reduced_chr5_len[i,2]/1e6,5, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr6_len)[1]){
	points(reduced_chr6_len[i,2]/1e6,6, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr7_len)[1]){
	points(reduced_chr7_len[i,2]/1e6,7, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr8_len)[1]){
	points(reduced_chr8_len[i,2]/1e6,8, col="darkred",pch=18,cex=0.5)
	}
for(i in 1:dim(reduced_chr9_len)[1]){
	points(reduced_chr9_len[i,2]/1e6,9, col="darkred",pch=18,cex=0.5)
	}
	for(i in 1:dim(reduced_chr10_len)[1]){
	points(reduced_chr10_len[i,2]/1e6,10, col="darkred",pch=18,cex=0.5)
	}
#dev.off()
graphics.off()