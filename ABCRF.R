setwd("C:/Users/momiglia/Dropbox/LocalDocuments/Finland work/Bird_paper/FSC/")
# load required libraries

library(abc)    # for approximate Bayesian computation functions
require(hexbin) # for plots of PCA
require(grid)   # for plots of PCA
require(abcrf)           # abc random forest: model choice
require(quantregForest)  # abc random forest: parameter estimation
library(plyr)
library(faraway)
library(scales)


num_of_threads=3 # This cause i have a 4 core , need one for other stuff



# Load data and prepare reference tables and target summary stats ---------------------------------------------------------------

# Load target summary statistics, and their names (names needed to split reference tabkes later)
MSAT_targetStats <-read.table ("./MSAT/OBS_DAT/OBSoutSumStats.txt", header = T)
mtDNA_targetStats <-read.table ("./mtDNA/OBS_STATS/outSumStats.txt", header = T)
MSAT_SS_names <- scan(file="./MSAT/OBS_DAT/OBSoutSumStats.txt",what=character(),nlines=1)
mtDNA_SS_names <- scan(file="./mtDNA/OBS_STATS/outSumStats.txt", what=character(),nlines=1)


# Load reference table for every simulation model (IM, SI, RI) and both mtDNA and MSAT data
SI_REF_MSAT<-read.table("./MSAT/SIM_STATS/SI_reftable.txt", header = T)
SC_REF_MSAT<-read.table("./MSAT/SIM_STATS/SC_reftable.txt", header = T)
RI_REF_MSAT<-read.table("./MSAT/SIM_STATS/RI_reftable.txt", header = T)
SI_REF_mtDNA<-read.table("./mtDNA/SIM_STATS/SI_reftable.txt", header = T)
SC_REF_mtDNA<-read.table("./mtDNA/SIM_STATS/SC_reftable.txt", header = T)
RI_REF_mtDNA<-read.table("./mtDNA/SIM_STATS/RI_reftable.txt", header = T)

# Add model codes
RI_REF_mtDNA$model<-"RI"
SI_REF_mtDNA$model<-"SI"
SC_REF_mtDNA$model<-"SC"
RI_REF_MSAT$model<-"RI"
SI_REF_MSAT$model<-"SI"
SC_REF_MSAT$model<-"SC"

# Provide color codes to plot 
RI_REF_mtDNA$col<-"red"
SI_REF_mtDNA$col<-"blue"
SC_REF_mtDNA$col<-"green"
RI_REF_MSAT$col<-"red"
SI_REF_MSAT$col<-"blue"
SC_REF_MSAT$col<-"green"


# Merge reference tables from single models into one ref table for MSAT and one for mtDNA 
MSAT_ref_table<-rbind.fill(SI_REF_MSAT, RI_REF_MSAT, SC_REF_MSAT)
mtDNA_ref_table<-rbind.fill(SI_REF_mtDNA, RI_REF_mtDNA, SC_REF_mtDNA)

modindex <- as.factor(mtDNA_ref_table$model)
mtDNA_sumstats   <- mtDNA_ref_table[,mtDNA_SS_names]
MSAT_sumstats   <- MSAT_ref_table[,MSAT_SS_names]



# ABC-RF Model Choice and Checking ----------------------------------------

# Constructs a random forest from a reference table towards performing an ABC model choice
# 2000 trees should be more than enough

# MSAT datatset
# Nref vector = specify the number of simulations in your reftable.txt and paste as many times as you want iterations. Here 10 iterations.
nscenarios = 3
output_MSAT <- data.frame(matrix(ncol=nscenarios+4))
output_file_MSAT <- "RF_MSAT.txt"
for (i in 1:10) {
        
MSAT_mc.rf <- abcrf(modindex~.,
                  data=data.frame(modindex, MSAT_sumstats),
                  ntree=1000,
                  paral=T,
                  ncores=num_of_threads)
MSAT_pred.rf        <- predict(object         = MSAT_mc.rf,
                             obs            = MSAT_targetStats,
                             training       = data.frame(modindex, MSAT_targetStats),
                             ntree          = 1000,
                             paral          = T,
                             ncores         = num_of_threads,
                             paral.predict  = T,
                             ncores.predict = num_of_threads)

output_MSAT <- rbind(output_MSAT, c(30000,
                          round(as.numeric(MSAT_pred.rf[length(MSAT_pred.rf)]),4),
                          unlist(MSAT_pred.rf[2:(length(MSAT_pred.rf)-1)]),
                          round(as.numeric(MSAT_mc.rf$prior.err),4)))

}
output_MSAT <- output_MSAT[-1,]
#on met les bons noms de colonnes à output:
colnames(output_MSAT) <- c("Nref",  "post_prob", rownames(as.data.frame(unlist(MSAT_pred.rf[2:(length(MSAT_pred.rf)-1)]))), "prior_err")

write.table(output_MSAT, output_file_MSAT, quote=FALSE, row.names=FALSE)


# Now the same, but using the mtDNA
nscenarios = 3
output_mtDNA <- data.frame(matrix(ncol=nscenarios+4))
output_file_mtDNA <- "RF_mtDNA.txt"
for (i in 1:10) {
        
        mtDNA_mc.rf <- abcrf(modindex~.,
                            data=data.frame(modindex, mtDNA_sumstats),
                            ntree=1000,
                            paral=T,
                            ncores=num_of_threads)
        mtDNA_pred.rf        <- predict(object         = mtDNA_mc.rf,
                                       obs            = mtDNA_targetStats,
                                       training       = data.frame(modindex, mtDNA_targetStats),
                                       ntree          = 1000,
                                       paral          = T,
                                       ncores         = num_of_threads,
                                       paral.predict  = T,
                                       ncores.predict = num_of_threads)
        
        output_mtDNA <- rbind(output_mtDNA, c(30000,
                                            round(as.numeric(mtDNA_pred.rf[length(mtDNA_pred.rf)]),4),
                                            unlist(mtDNA_pred.rf[2:(length(mtDNA_pred.rf)-1)]),
                                            round(as.numeric(mtDNA_mc.rf$prior.err),4)))
        
}
output_mtDNA <- output_mtDNA[-1,]
#on met les bons noms de colonnes à output:
colnames(output_mtDNA) <- c("Nref",  "post_prob", rownames(as.data.frame(unlist(mtDNA_pred.rf[2:(length(mtDNA_pred.rf)-1)]))), "prior_err")

write.table(output_mtDNA, output_file_mtDNA, quote=FALSE, row.names=FALSE)


# Now put mean results and SD in a single table



output_all <- data.frame(matrix(ncol=nscenarios+5, nrow = 2))

colnames(output_all)<-c("Data","Nref",  "post_prob", rownames(as.data.frame(unlist(mtDNA_pred.rf[2:(length(mtDNA_pred.rf)-1)]))), "prior_err")
output_all[1,1]<-"MSAT"
output_all[2,1]<-"mtDNA"
output_all[,2]<-30000
output_all[1,3]<-mean(output_MSAT$post_prob)
output_all[2,3]<-mean(output_mtDNA$post_prob)
output_all[1,4]<-mean(output_MSAT$allocation)
output_all[2,4]<-mean(output_mtDNA$allocation)
output_all[1,5]<-mean(output_MSAT$vote1)
output_all[2,5]<-mean(output_mtDNA$vote1)
output_all[1,6]<-mean(output_MSAT$vote2)
output_all[2,6]<-mean(output_mtDNA$vote2)
output_all[1,7]<-mean(output_MSAT$vote3)
output_all[2,7]<-mean(output_mtDNA$vote3)
output_all[1,8]<-mean(output_MSAT$prior_err)
output_all[2,8]<-mean(output_mtDNA$prior_err)
output_all$postSD[1]<-sd(output_MSAT$post_prob)
output_all$postSD[2]<-sd(output_mtDNA$post_prob)
output_all$priorSD[1]<-sd(output_MSAT$prior_err)
output_all$priorSD[2]<-sd(output_mtDNA$prior_err)
output_all$allocation[which(output_all$allocation==1)]<-"RI"
write.table(output_all, "output_all.txt", quote=FALSE, row.names=FALSE)






# Bla ---------------------------------------------------------------------

output_MSAT<- read.table("./RF_MSAT.txt", header = T)

sd(output_MSAT$)
# Print the prior error rate
MSAT_model_RF$prior.err

# Make a figure showing a) variables' relative 
# importance as well as b) LDA showing simulations and observed data
pdf("LDA_MSAT.pdf", width= 6, height = 5)
plot(MSAT_mc.rf,
     training=data.frame(modindex, MSAT_sumstats),
     obs=MSAT_targetStats)
dev.off()

pdf("Prior_ErrRate_MSAT.pdf", width= 6, height = 4)
err.abcrf(MSAT_mc.rf,
          training=data.frame(modindex, MSAT_sumstats),
          paral=T,
          ncores=num_of_threads)
dev.off()


MSAT_ModSelRes_RF

# LDA and error plto for mtDNA

mtDNA_model_RF <- abcrf(modindex~.,
                       data=data.frame(modindex, mtDNA_sumstats),
                       ntree=2000,
                       paral=T,
                       ncores=num_of_threads)

# Print the prior error rate
mtDNA_model_RF$prior.err

# Make a figure showing a) variables' relative 
# importance as well as b) LDA showing simulations and observed data
pdf("LDA_mtDNA.pdf", width= 6, height = 5)
plot(mtDNA_model_RF,
     training=data.frame(modindex, mtDNA_sumstats),
     obs=mtDNA_targetStats)
dev.off()

pdf("Prior_ErrRate_mtDNA.pdf", width= 6, height = 4)
err.abcrf(mtDNA_model_RF,
          training=data.frame(modindex, mtDNA_sumstats),
          paral=T,
          ncores=num_of_threads)
dev.off()


# LDA and error plto for MSAT


MSAT_model_RF <- abcrf(modindex~.,
                        data=data.frame(modindex, MSAT_sumstats),
                        ntree=2000,
                        paral=T,
                        ncores=num_of_threads)

# Print the prior error rate
MSAT_model_RF$prior.err

# Make a figure showing a) variables' relative 
# importance as well as b) LDA showing simulations and observed data
pdf("LDA_MSAT.pdf", width= 6, height = 5)
plot(MSAT_model_RF,
     training=data.frame(modindex, MSAT_sumstats),
     obs=MSAT_targetStats)
dev.off()

pdf("Prior_ErrRate_MSAT.pdf", width= 6, height = 4)
err.abcrf(MSAT_model_RF,
          training=data.frame(modindex, MSAT_sumstats),
          paral=T,
          ncores=num_of_threads)
dev.off()












