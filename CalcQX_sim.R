#!/usr/bin/env Rscript
library("optparse")

option_list = list(
              make_option(c("-w", "--gwasfile"), type="character", default=NULL, help="GWAS input file name"),
              make_option(c("-e", "--neutfile"), type="character", default=NULL, help="Neutral input file name"),
              make_option(c("-o", "--outfile"), type="character", default="output.txt", help="Output file name"),
              make_option(c("-s", "--scorefile"), type="character", default=NULL, help="Score file name")
    );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gwasfile)){
    print <- help(opt_parser)
    stop("GWAS file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$neutfile)){
    print <- help(opt_parser)
    stop("Neutral file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$outfile)){
    print <- help(opt_parser)
    stop("Output file name must be supplied.n", call.=FALSE)
}

if (is.null(opt$scorefile)){
    print <- help(opt_parser)
    stop("Score file name must be supplied.n", call.=FALSE)
}


gwasfile <- opt$gwasfile
neutfile <- opt$neutfile
outfile <- opt$outfile
scorefile <- opt$scorefile
pseudorep <- opt$numrep
popmatch <- opt$popmatch
emprep <- opt$emprep

ComputeFmat <- function(neut_leaves_freqs){
    leaves <- colnames(neut_leaves_freqs)

    checksegneut <- which( apply(neut_leaves_freqs,1,sum)/length(leaves) < 0.95  & 
                             apply(neut_leaves_freqs,1,sum)/length(leaves) > 0.05 )
    neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]
    print('Neutral SNPs')
    print(nrow(neut_leaves_freqs))
    
    neut_leaves_freqs_means <- apply(neut_leaves_freqs, 1, mean)
    mean_hetero <- neut_leaves_freqs_means*(1-neut_leaves_freqs_means)
    numSNPs <- length(neut_leaves_freqs_means)
    
    Fmat <- sapply(seq(1,dim(neut_leaves_freqs)[2]),function(x){
        sapply(seq(1,dim(neut_leaves_freqs)[2]),function(y){
            cov(neut_leaves_freqs[,x]/sqrt(mean_hetero), neut_leaves_freqs[,y]/sqrt(mean_hetero))
        })
    })

    colnames(Fmat) <- colnames(neut_leaves_freqs)
    rownames(Fmat) <- colnames(neut_leaves_freqs)

    #print(Fmat)
    
    return(Fmat)
}


ChiSquared <- function(leaves_freqs,effects,F_mat,randomize=FALSE){

    leaves <- colnames(leaves_freqs)

    # Average between 0.05 and 0.95
    checkseg <- which( apply(leaves_freqs,1,sum)/length(leaves) < 0.95  & apply(leaves_freqs,1,sum)/length(leaves) > 0.05 )

    leaves_freqs <- leaves_freqs[checkseg,]
    
    effects <- effects[checkseg]
    print('Candidates:')
    print(nrow(leaves_freqs))
    
    # Randomize effects if necessary
    if(randomize == TRUE){effects <- effects * sample(c(-1,1),length(effects),replace=TRUE)}
    
    # Compute mean genetic values
    meangen <- apply(leaves_freqs * effects, 2, function(x){sum(x)})
    
    # Scale by average genetic value
    meangen <- (meangen - mean(meangen))
    
    # Compute the estimated ancestral genetic variance over all populations
    meanfreqs <- apply(leaves_freqs,1,mean)
    
    varmean <- sum(sapply(seq(1,length(meanfreqs)),function(i){
        score = meanfreqs[i]*(1-meanfreqs[i])*effects[i]^2
        return(score)
    }))

    meangenvec <- 2*meangen / sqrt( 4*varmean)

    # Compute Q_X statistic
    numerator <- t(meangen) %*% solve(Fmat) %*% meangen
    denominator <- varmean
    Qteststat <- numerator / denominator
    
    Pval <- 1 - pchisq(Qteststat,qr(Fmat)$rank)
    allstats <- c(Qteststat,Pval)

    return(list(allstats,meangenvec, nrow(leaves_freqs)))
}


print("Loading data...")
data <- read.table(gwasfile, h=T)
#data <- read.table('tester.txt', h=T)

leaves_freqs <- as.matrix(data[,c(4,5)])
effects <- as.numeric(data[,3])

# Load Neutral data
print('Loading neutral...')
neutdata <- read.table(neutfile, h=T)
#neutdata <- read.table('NEUT_tester.txt', h=T)
neut_leaves_freqs <- as.matrix(neutdata[,c(4, 5)])

# Calculate covariance matrix
print("Computing covariance matrix...")
Fmat <- ComputeFmat(neut_leaves_freqs)
print(Fmat)

# Calculate chi-squared statistics
print("Computing Q_X statistic...")
totaltest <- ChiSquared(leaves_freqs,effects,Fmat,randomize=FALSE)
totalstat <- totaltest[[1]]
meangenvec <- totaltest[[2]]
qtab <- cbind(round(totalstat[1],3),totalstat[2])
colnames(qtab) <- c("Q_X","Pval")

write("Genetic scores",file=outfile,append=FALSE)
write(paste(names(meangenvec),collapse="\t"),file=outfile,append=TRUE)
write(paste(meangenvec,collapse="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

printgenvec <- cbind(names(meangenvec),meangenvec)
rownames(printgenvec) <- c()
colnames(printgenvec) <- c("#POP","SCORE")
write.table(printgenvec,file=scorefile,sep="\t",row.names=FALSE,col.names=TRUE, quote=FALSE,append=FALSE)

write(paste("Q_X:",qtab[1],sep="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

write(paste("Chi-squared distribution P-value:",qtab[2],sep="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

write(paste("Nr of SNPs used:",totaltest[[3]],sep="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

