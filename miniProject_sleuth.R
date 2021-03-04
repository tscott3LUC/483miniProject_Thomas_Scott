#load in the sleuth package
library(sleuth)

#read in the table describing the samples
stab<-read.table("miniProject_table.txt",header=TRUE,stringsAsFactors=FALSE)

#initialize sleuth object
so<-sleuth_prep(stab)

#fit a model comparing the two conditions
so<-sleuth_fit(so,~condition,'full')

#fit the reduced model to compare the likelihood ratio test
so<-sleuth_fit(so,~1,'reduced')

#perform the likelihood ratio test for differential expression
so<-sleuth_lrt(so,'reduced','full')

#load in dplyr for filtering data.frames
library(dplyr)
rm 
#extract results from the sleuth object
sleuth_table<-sleuth_results(so,'reduced:full','lrt',show_all=FALSE)

#filter significant (FDR<0.05) results, sorting by pval
sleuth_significant<-dplyr::filter(sleuth_table,qval<=0.05)%>%dplyr::arrange(pval)

#select by column header names to get the requested details:
#target_id, test_stat, pval, and qval
sleuth_significant_selected<-dplyr::select(sleuth_significant, target_id, test_stat, pval, qval)

#Write the target_id, test_stat, pval, and qval for each significant transcript to the log file (tab delimited)
write.table(sleuth_significant_selected,file='sleuthResults.txt',quote=FALSE,row.names=FALSE,sep='\t')