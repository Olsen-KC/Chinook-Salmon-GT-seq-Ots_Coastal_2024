library(tidyverse)
library(magrittr)
library(adegenet)

##read in original genotypes for 2023 samples
genos_0.1<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/Ots_coastal_chinook_2023_GTs_0.1.csv")

##Problem with how the command above reads the file in. Adds a comma to the sex genotype


unique(genos_0.1$`Ots_SEXY3-1`)

##fix 

genos_0.1 %<>%
  mutate(`Ots_SEXY3-1` = case_when(`Ots_SEXY3-1` == "XX," ~ 'XX',
                                   `Ots_SEXY3-1` == "XY," ~ 'XY', 
                                   `Ots_SEXY3-1` == "00," ~ '00'
  ))

##read in marker info for 2023 samples

marker_info<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/marker_info.csv")

marker_info %<>%
  mutate(a1_count =  as.numeric(substr(a1_count, 3, nchar(a1_count)))) %>%
  mutate(a2_count =  as.numeric(substr(a2_count, 3, nchar(a2_count)))) %>%
  mutate(ind = str_remove(ind, "^\\./")) %>%
  mutate(ind = str_remove(ind, "\\.genos"))



##create new column identifying positive and negative controls
##clean up sample name 
##identify duplicate positive controls with adapter in name
genos_0.1 %<>%
  mutate(control = case_when(str_detect(Sample, "negative") ~ 'negative',
                             str_detect(Sample, "positive") ~ 'positive', 
                                 TRUE ~ 'sample')) %>%
  mutate(adapter = str_extract(Sample, "[ATCG]{6}-[ATCG]{6}")) %>%
  mutate(sample_simple = str_extract(Sample, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")) %>%
  relocate(Sample, sample_simple, adapter, control)

##plot on target reads for negative, positive controls and samples##
ggplot()+geom_histogram(data = genos_0.1, aes(x = `On-Target Reads`, fill= control)) + theme_classic()+scale_fill_viridis_d()

##pot % of reads on target
ggplot()+geom_histogram(data = genos_0.1, aes(x = `%On-Target`, fill = control)) + theme_classic()+scale_fill_viridis_d()


##remove negative and positive controls and control column

genos_0.11<-genos_0.1 %>%
  filter(control == "sample") %>%
  select(-control)

##grab duplicated samples 
genos_0.11 %<>% 
  group_by(sample_simple) %>%
  slice_max(order_by = `On-Target Reads`, n = 2)

dups <- genos_0.11[genos_0.11$sample_simple %in% genos_0.11$sample_simple[duplicated(genos_0.11$sample_simple)],]
dups <- dups[order(dups$sample_simple),]


##calculate concordance between duplicate pairs
dups_genos <- dups[,9:ncol(dups)] 
rep_info <- matrix(ncol=ncol(dups_genos), nrow=nrow(dups_genos)/2)
colnames(rep_info) <- colnames(dups_genos)
for (j in 1:(nrow(dups_genos)/2)) {
  for (i in 1:ncol(dups_genos)) {
    rep_info[j,i] <- sum(dups_genos[(j*2)-1,i]==dups_genos[(j*2),i])
  }
}

geno_concordance <- as.data.frame(as.matrix(rep_info)) %>%
  rowMeans()

rep_data <- as.data.frame(cbind(dups[c(1:length(geno_concordance))*2,1], geno_concordance))

#plot concordance between duplicate pairs
ggplot(data=rep_data)+geom_histogram(aes(x=geno_concordance))+theme_classic()

##most have high concordance but some pairs do not 
##look at pairs with less than 80% concordance 

low_concordance<-rep_data[which(rep_data$geno_concordance < 0.8), ]

##27 pairs with concordance below 80% 
##of those 27, 25 are carcass samples so low concordance is likely due to degradation

#this writes a new dataset (0.2) by choosing the samples within duplicates and keeping the one with the highest genotyping success
genos_0.2 <- genos_0.11 %>%
  group_by(sample_simple) %>%
  filter(`On-Target Reads` == max(`On-Target Reads`))


##look at sex marker

sex_marker_info <- marker_info %>%
  filter(str_detect(marker, "SEXY")) %>%
  mutate(called_geno = replace_na(called_geno, "00")) %>%
  mutate(called_geno = case_when(called_geno == "A1HOM" ~ "XX",
                                 called_geno == "HET" ~ "XY",
                                 called_geno == "00" ~ "00"))

##plot XY depth and genotypic sex call
ggplot(data = sex_marker_info)+geom_point(aes(a2_count, a1_count, color = called_geno))+scale_colour_viridis_d(name = "Called Sex Genotype")+theme_classic()+
  xlab("Y-specific probe count")+ylab("Theoretical X chromosome count")+geom_abline(aes(intercept = 0, slope = 0.1))+geom_abline(aes(intercept = 0, slope = 5))+
  geom_abline(aes(intercept = 0, slope = 10))+geom_abline(aes(intercept = 0, slope = 0.2))+xlim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+
  ylim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+geom_abline(aes(intercept = 0, slope = 1), color = "darkred")

##single sample has very high allele 1 and allele 2 counts
##adjust x and y limits 
ggplot(data = sex_marker_info)+geom_point(aes(a2_count, a1_count, color = called_geno))+scale_colour_viridis_d(name = "Called Sex Genotype")+theme_classic()+
  xlab("Y-specific probe count")+ylab("Theoretical X chromosome count")+geom_abline(aes(intercept = 0, slope = 0.1))+geom_abline(aes(intercept = 0, slope = 5))+
  geom_abline(aes(intercept = 0, slope = 10))+geom_abline(aes(intercept = 0, slope = 0.2))+xlim(0,max(5000))+
  ylim(0,max(5000))+geom_abline(aes(intercept = 0, slope = 1), color = "darkred")

##Ots sex marker looks good
##does not need correcting


##visualize IFI scores, individual call rates, and locus call rates##
ggplot(genos_0.2)+geom_histogram(aes(x=IFI))+geom_vline(aes(xintercept= 2.5), color="red")+theme_classic()

ggplot(genos_0.2)+geom_histogram(aes(x=`%GT`))+geom_vline(aes(xintercept= 90), color="red")+theme_classic()

##quite a few individuals with high IFI and low genotyping call rates##

##calculate marker call rates##
##grab genotypes

genotypes<-as.matrix(genos_0.2[, c(9:ncol(genos_0.2))])
rownames(genotypes)<-genos_0.2$sample_simple

##sum proportion of missing genotypes for each marker
marker_miss<-as.data.frame(apply(genotypes, MARGIN = 2, function(x){ sum(x == "00" | x == "0")}) / nrow(genotypes))
colnames(marker_miss)<-"proportion_missing"

ggplot(marker_miss) + geom_histogram(aes(x=proportion_missing))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")

##several markers with low call rates


##remove samples with call rates below 70%
##98 samples have call rates below 70%
genos_0.3 <- genos_0.2 %>%
  filter(`%GT` > 70)

##recalculate marker call rates##
##grab genotypes

genotypes<-as.matrix(genos_0.3[, c(9:ncol(genos_0.3))])
rownames(genotypes)<-genos_0.3$sample_simple

##sum proportion of missing genotypes for each marker
marker_miss<-as.data.frame(apply(genotypes, MARGIN = 2, function(x){ sum(x == "00" | x == "0")}) / nrow(genotypes))
colnames(marker_miss)<-"proportion_missing"

ggplot(marker_miss) + geom_histogram(aes(x=proportion_missing))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")

##drop two markers missing in more than 50% of samples (Ots17_1486479_C6, Ots_wenYhap_33126)

genos_0.3 %<>%
  select(-one_of(c("Ots17_1486479_C6", "Ots_wenYhap_33126")))


##recalculate IFI scores##

IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.3)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)



IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.3 <- genos_0.3 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))

range(genos_0.3$IFI)

##0 samples have IFI scores greater than 10##

##recalculate individual call rates##

genos_0.3 %<>% #dont forget to change to the correct dataset
  rowwise %>%
  mutate(`%GT` = 100*(1-sum((cur_data() == "00"), na.rm = TRUE)/(ncol(.)-8))) %>% #change number of columns if needed (8 meta columns in this file)
  ungroup() #always ungroup after rowwise

missing_ind<- which(genos_0.3$`%GT` < 90)
length(missing_ind)


##49 individuals have call rates below 90%##

#remove these 49 individuals
genos_0.4 <- genos_0.3 %>%
  filter(`%GT` > 90)

##recalculate marker call rates##
##grab genotypes

genotypes<-as.matrix(genos_0.4[, c(9:ncol(genos_0.4))])
rownames(genotypes)<-genos_0.4$sample_simple

##sum proportion of missing genotypes for each marker and divide by number of samples (rows)
marker_miss<-as.data.frame(apply(genotypes, MARGIN = 2, function(x){ sum(x == "00" | x == "0")}) / nrow(genotypes))
colnames(marker_miss)<-"proportion_missing"

ggplot(marker_miss) + geom_histogram(aes(x=proportion_missing))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")

##all markers called in greater than 80% of samples

##recalculate IFI scores##

IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)


IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.4 <- genos_0.4 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))

range(genos_0.4$IFI)
hist(genos_0.4$IFI)

##12 samples have IFI greater than 2.5
##Of these 12, 9 are carcass samples again pointing to degredation

##remove these 12 samples 

genos_0.4 %<>%
  filter(IFI < 2.5)

##recalculate marker call rates##
##grab genotypes

genotypes<-as.matrix(genos_0.4[, c(9:ncol(genos_0.4))])
rownames(genotypes)<-genos_0.4$sample_simple

##sum proportion of missing genotypes for each marker and divide by number of samples (rows)
marker_miss<-as.data.frame(apply(genotypes, MARGIN = 2, function(x){ sum(x == "00" | x == "0")}) / nrow(genotypes))
colnames(marker_miss)<-"proportion_missing"

marker_miss %<>%
  rownames_to_column("sample")

ggplot(marker_miss) + geom_histogram(aes(x=proportion_missing))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")

##At this point the dataset has 387 Individuals and 350 loci and a sex marker##
##Individual call rates range 91 - 100%; Marker call rates range 81 - 100%
##IFI ranges 0.10 - 2.44##


##Evaluate possible paralogs##

#get marker names of markers with 0.1 > missingness < 0.2
miss_mod <- marker_miss[which(marker_miss$proportion_missing > 0.1 & marker_miss$proportion_missing < 0.2), 1]


hets <- filter(marker_info, called_geno == "HET" | is.na(called_geno))

##because a lot of samples in the original marker_info file ended up being removed during quality filtering
##I added an additional filter here that also removed these individuals when assessing potential paralogs 
##otherwise the following script identifies an inordinate number of makrers (195) as potential paralogs when in reality it's just the poor samples making these markers look problematic


models <- hets %>%
  mutate(sample_simple = str_extract(ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")) %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(sample_simple %in% genos_0.4$sample_simple) %>%    ##remove samples that were removed during quality filtering for this calculation
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  group_map(~ lm(a1_count ~ a2_count, data= .))


lms <- lapply(models, coef)
ggplot()+geom_histogram(aes(x = sapply(lms,`[`,2)))+theme_classic()+ggtitle("allele ratios for all NA and HET calls")+geom_vline(aes(xintercept = 1.5), color = "red", linetype = 2)+geom_vline(aes(xintercept = (2/3)), color = "red", linetype = 2)+xlab("allele ratio (a1/a2)")+geom_vline(aes(xintercept = 1), color = "black")

lms_anova <- lapply(models, summary)

# collect info about each bad model
paralog_possible <- which(abs(sapply(lms,`[`,2)) > 1.5) #bad because a positively skewed allele ratio
paralog_possible2 <- which(abs(sapply(lms,`[`,2)) < (2/3)) # bad because a negative skewed allele ratio

paralog_possible3 <- which(sapply(lms_anova, function(x) x$coefficients[,4][2])> 0.01) # bad because too much variance in allele ratio, even if mean ratio is 1

paralog_possible <- c(paralog_possible, paralog_possible2, paralog_possible3)

plots <- marker_info %>%
  mutate(sample_simple = str_extract(ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")) %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(sample_simple %in% genos_0.4$sample_simple) %>%    ##remove samples that were removed during quality filtering for these plots
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  do(plots=ggplot(data=.)+geom_point(aes(a1_count, a2_count, color = called_geno))+theme_classic()+geom_abline(aes(slope=1, intercept=0))+geom_abline(aes(slope = 10, intercept=0), color = "green")+geom_abline(aes(slope = 0.1, intercept=0), color = "red")+geom_abline(aes(slope = 0.2, intercept=0), color = "blue")+geom_abline(aes(slope = 5, intercept=0), color = "blue")+coord_equal(ratio=1)+geom_abline(slope = -1, intercept = 10)+ggtitle(unique(.$marker)))

#plot all "bad markers"

#first add the missningness markers to the list to examine
mod_bad_plot_index <- which(plots$marker %in% miss_mod)
paralog_possible <- c(mod_bad_plot_index, paralog_possible)

# then loop through the plots by changing the index until you have looked at all your questionable markers
plots$plots[[paralog_possible[1]]]

to_filt <- c("Ots17_1066109_C6", "Ots_CHI06105101_16717", "Ots_110495-380", "Ots_GPDH-338") 

##remove 4 loci that are possible paralogs##
genos_0.5 <- genos_0.4 %>%
  dplyr::select(-one_of(to_filt))

##identify and drop monomorphic markers
##remove "dots" from column names##

geno_genind <- genos_0.5
colnames(geno_genind) <- gsub("\\.", "_", colnames(geno_genind))


geno_genind <- as.matrix(geno_genind[, c(9:(ncol(geno_genind)))])

row.names(geno_genind) <- genos_0.5$sample_simple


geno_genind<-df2genind(geno_genind, sep = "", ploidy = 2, NA.char = "00")


mono<-as.data.frame(which(geno_genind$loc.n.all < 2))
colnames(mono)<-"marker#"

mono %<>%
  rownames_to_column("marker") 

to_drop<-c(mono$marker[c(1, 3:7, 9:12)], "Ots_IGF-I.1-76", "Ots_u07-64.221")

##Normally would remove these but they may be polymorphic when integrating samples from other years so will leave for now and 
##reassess monomorphic markers after combining data across years. 



##identify duplicated sampling events with relatedness estimates in coancestry
##format for coancestry##

##add a delimiter to genotypes so you can split into 2 columns
##grab genotypes
replacedGT<-as.data.frame(genos_0.5[, c(9:(ncol(genos_0.5)-1))])

##add delimiter
delim_gtypes <- replacedGT[, ] %>%
  mutate_all(funs(case_when(. == "AA" ~ "A/A",
                            . == "AT" ~ "A/T", 
                            . == "AC" ~ "A/C", 
                            . == "AG" ~ "A/G", 
                            . == "TT" ~ "T/T", 
                            . == "TA" ~ "T/A", 
                            . == "TC" ~ "T/C", 
                            . == "TG" ~ "T/G",
                            . == "CC" ~ "C/C", 
                            . == "CA" ~ "C/A", 
                            . == "CT" ~ "C/T", 
                            . == "CG" ~ "C/G", 
                            . == "GG" ~ "G/G", 
                            . == "GA" ~ "G/A", 
                            . == "GT" ~ "G/T", 
                            . == "GC" ~ "G/C", 
                            . == "--" ~ "-/-", 
                            . == "A-" ~ "A/-", 
                            . == "T-" ~ "T/-", 
                            . == "C-" ~ "C/-", 
                            . == "G-" ~ "G/-", 
                            . == "-A" ~ "-/A",
                            . == "-T" ~ "-/T", 
                            . == "-C" ~ "-/C", 
                            . == "-G" ~ "-/G", 
                            . == "00" ~ "0/0", 
                            . == "0" ~ "0/0"
  )))

library(splitstackshape)


vrtnames<-c(colnames(delim_gtypes))

##split genotypes into 2 alleles

dfspltgeno<-as.data.frame(cSplit(delim_gtypes, splitCols = vrtnames, sep = "/", type.convert = FALSE))


##coancestry has character limits for sample names 
##short sample names
coancestrylbls<-data.frame(c(genos_0.5$sample_simple))

colnames(coancestrylbls)<-c("Ind")

coancestrylbls %<>%
  mutate(Ind = substr(Ind, 4, 16))

##add shortened sample names to rownames
rownames(dfspltgeno)<-as.character(coancestrylbls$Ind)

##replace nucleotides with numeric values for coancestry## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###

dfspltgeno %<>%
  mutate_all(funs(case_when(. == "A" ~ 1,
                            . == "T" ~ 2,
                            . == "C" ~ 3,
                            . == "G" ~ 4, 
                            . == "-" ~ 5,
                            . == "0" ~ 0, 
  )))


##write coancestry input file
write.table(dfspltgeno, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/Coastal chinook 2023 coancestry input.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)

##read in relatedness estimates after running coancestry##


rel<-read.table("C:/ZSL/Coancestry/Coastal chinook 2023/RelatednessEstimates.txt", sep = ",")

##grab lynchR estimator##

rel<-rel[, c(1:3, 8)]

colnames(rel)<-c("dyad", "Ind1", "Ind2", "lynchr")


##grab sample pairs with high relatedness and duplicate genotypes

high_rel<-rel[which(rel$lynchr > 0.95), ]

##add "Ots" prefix back to sample names
high_rel %<>%
  mutate(Ind1 = paste("Ots", Ind1, sep = "")) %>%
  mutate(Ind2 = paste("Ots", Ind2, sep = "")) 

##grab call rates for these samples
Ind1_GT<-genos_0.5[genos_0.5$sample_simple %in% high_rel$Ind1, c(2, 7)]
Ind2_GT<-genos_0.5[genos_0.5$sample_simple %in% high_rel$Ind2, c(2, 7)]

##put call rates in same order as high_rel dataframe
Ind1_GT<-Ind1_GT[match(high_rel$Ind1, Ind1_GT$sample_simple), ]
Ind2_GT<-Ind2_GT[match(high_rel$Ind2, Ind2_GT$sample_simple), ]

##check they are in same order
identical(high_rel$Ind1, Ind1_GT$sample_simple)
identical(high_rel$Ind2, Ind2_GT$sample_simple)

colnames(Ind1_GT)<-c("sample_simple_1", "%GT_Ind1")
colnames(Ind2_GT)<-c("sample_simple_2", "%GT_Ind2")


##add call rates to high_rel df
high_rel<-cbind(high_rel, Ind1_GT[, 2], Ind2_GT[, 2])

high_rel %<>%
  relocate(dyad, Ind1, `%GT_Ind1`, Ind2, `%GT_Ind2`, lynchr)

##check which Ind1 are also in Ind2 

intersect(high_rel$Ind1, high_rel$Ind2)

##1 sample is in both OtsAC23CHER_0045##

write.csv(high_rel, quote = FALSE, row.names = FALSE, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/coastal chinook 2023 high relatedness duplicates.csv")

##36 instances of high pairwise relatedness indicative of duplicate genotypes
##remove samples that make up more than one duplicate pair
##For each pair remove sample with lower genotype call rate

to_remove<-c("OtsAC23CHER_1035", "OtsAC23CHER_0004", "OtsAC23CHER_0007", "OtsAC23CHER_1044", "OtsAC23CHER_0013", "OtsAC23CHER_1001", "OtsAC23CHER_1049", "OtsAC23CHER_1038", "OtsAC23CHER_1009", 
             "OtsAC23CHER_1045", "OtsAC23CHER_1004", "OtsAC23CHER_1032", "OtsAC23CHER_1050", "OtsAC23CHER_0027", "OtsAC23CHER_0028", "OtsAC23CHER_1018", "OtsAC23CHER_0034", "OtsAC23CHER_1017", 
             "OtsAC23CHER_0036", "OtsAC23CHER_1023", "OtsAC23CHER_0039", "OtsAC23CHER_0041", "OtsAC23CHER_0044", "OtsAC23CHER_0045", "OtsAC23CHER_1007", "OtsAC23CHER_1015", "OtsAC23CHER_1011", 
             "OtsAC23CHER_0053", "OtsAC23CHER_1027", "OtsAC23CHER_1039", "OtsAC23CHER_0063", "OtsAC23CHER_1031", "OtsAC23CHER_1041", "OtsAC23CHER_1034", "OtsAC24ELKR_0016", "OtsCC23SUMP_0012")

##remove these 36 samples##
genos_2.0<-genos_0.5 %>%
  filter(! sample_simple %in% to_remove)


##check that all duplicates have been removed 


replacedGT<-as.data.frame(genos_2.0[, c(9:(ncol(genos_2.0)-1))])

##add delimiter
delim_gtypes <- replacedGT[, ] %>%
  mutate_all(funs(case_when(. == "AA" ~ "A/A",
                            . == "AT" ~ "A/T", 
                            . == "AC" ~ "A/C", 
                            . == "AG" ~ "A/G", 
                            . == "TT" ~ "T/T", 
                            . == "TA" ~ "T/A", 
                            . == "TC" ~ "T/C", 
                            . == "TG" ~ "T/G",
                            . == "CC" ~ "C/C", 
                            . == "CA" ~ "C/A", 
                            . == "CT" ~ "C/T", 
                            . == "CG" ~ "C/G", 
                            . == "GG" ~ "G/G", 
                            . == "GA" ~ "G/A", 
                            . == "GT" ~ "G/T", 
                            . == "GC" ~ "G/C", 
                            . == "--" ~ "-/-", 
                            . == "A-" ~ "A/-", 
                            . == "T-" ~ "T/-", 
                            . == "C-" ~ "C/-", 
                            . == "G-" ~ "G/-", 
                            . == "-A" ~ "-/A",
                            . == "-T" ~ "-/T", 
                            . == "-C" ~ "-/C", 
                            . == "-G" ~ "-/G", 
                            . == "00" ~ "0/0", 
                            . == "0" ~ "0/0"
  )))

library(splitstackshape)


vrtnames<-c(colnames(delim_gtypes))

##split genotypes into 2 alleles

dfspltgeno<-as.data.frame(cSplit(delim_gtypes, splitCols = vrtnames, sep = "/", type.convert = FALSE))


##coancestry has character limits for sample names 
##short sample names
coancestrylbls<-data.frame(c(genos_2.0$sample_simple))

colnames(coancestrylbls)<-c("Ind")

coancestrylbls %<>%
  mutate(Ind = substr(Ind, 4, 16))

##add shortened sample names to rownames
rownames(dfspltgeno)<-as.character(coancestrylbls$Ind)

##replace nucleotides with numeric values for coancestry## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###

dfspltgeno %<>%
  mutate_all(funs(case_when(. == "A" ~ 1,
                            . == "T" ~ 2,
                            . == "C" ~ 3,
                            . == "G" ~ 4, 
                            . == "-" ~ 5,
                            . == "0" ~ 0, 
  )))


##write coancestry input file
write.table(dfspltgeno, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/Coastal chinook 2023 coancestry input duplicates removed.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)

##read in relatedness estimates after running coancestry##


rel<-read.table("C:/ZSL/Coancestry/Coastal chinook 2023 duplicates removed/RelatednessEstimates.txt", sep = ",")

##grab lynchR estimator##

rel<-rel[, c(1:3, 8)]

colnames(rel)<-c("dyad", "Ind1", "Ind2", "lynchr")


##grab sample pairs with high relatedness and duplicate genotypes

high_rel<-rel[which(rel$lynchr > 0.95), ]

##no pairs with relatedness estimates greater than 0.95

##filtered dataset has 351 individuals genotyped at 346 markers and a sex marker

ggplot(genos_2.0)+geom_density(aes(x=`On-Target Reads`))+geom_vline(aes(xintercept=median(`On-Target Reads`)), color = "red") +theme_classic()

ggplot(genos_2.0)+geom_density(aes(x=`%On-Target`))+geom_vline(aes(xintercept=median(`%On-Target`)), color = "red") +theme_classic()

##Look at sequencing depth
marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  summarise(mean=mean(sumdepth, na.rm = TRUE), median=median(sumdepth, na.rm = TRUE), sd=sd(sumdepth, na.rm = TRUE))


save(genos_2.0, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/coastal_chinook_2023_QC_filtered_genotypes_2.0.R")

##now that the 2023 genotypes have been filtered 
##look at 2020 and 2021 genotyping 

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2020_genotyping_dayan/Ots_coastal_chinook_2020_2021_genos_2.2.R")

Ots_20_21<-genos_2.2
rm(genos_2.2)

dim(Ots_20_21)
colnames(Ots_20_21)

##977 individuals genotyped at 324 markers


load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/2023_genotyping/coastal_chinook_2023_QC_filtered_genotypes_2.0.R")

Ots_23<-genos_2.0

rm(genos_2.0)

##look at which columns/markers are in both datasets

intersect(colnames(Ots_23), colnames(Ots_20_21))

##321 markers in both datasets

##add Population column to 2023 dataset 

Ots_23 %<>%
  mutate(Population = substr(sample_simple, 8, 11)) %>%
  mutate(year = substr(sample_simple, 6, 7)) %>%
  relocate(Sample, sample_simple, adapter, `Raw Reads`, `On-Target Reads`, `%On-Target`, `%GT`,  IFI, Population)


##read in 2023 genotyping metadata 

coastal_2023_meta<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/Coastal Chinook 2023 metadata.csv")


##check if any of the samples in the genotype dataset are missing from the metadata##

which(! Ots_23$sample_simple %in% coastal_2023_meta$sample)

##remove samples from metadata that are not in filtered genotype dataset and drop the Population column which is already in the genotype dataset. 

coastal_2023_meta %<>%
  filter(sample %in% Ots_23$sample_simple) %>%
  select(! Population)


##put rows of metadata into same order as in the genotype dataframe

coastal_2023_meta<-coastal_2023_meta[order(match(coastal_2023_meta$sample, Ots_23$sample_simple)), ]

##check they are in the same order

identical(Ots_23$sample_simple, coastal_2023_meta$sample)

##now that the samples match, add the metadata to the genotypes 

Ots_23<-cbind(Ots_23, coastal_2023_meta)

Ots_23 %<>%
  select(! c(sample, adapter)) %>%
  relocate(sample_simple, Population, date, year, stream, location, run, sex, origin)

##combine 2020/2021 & 2024 genotyping efforts 

Ots_23 %<>%
  mutate(sample = sample_simple) %>%
  select(! c(Sample, sample_simple, year)) %>%
  relocate(sample, Population, date, stream, location, run, sex, meps, fl, origin)

Ots_21_23<-rbind(Ots_20_21[, 1:10], Ots_23[, 1:10])

colnames(Ots_21_23)<-c("sample", "Population", "sample_date", "stream", "location", "phenotypic_run", "phenotypic_sex", "meps_length", "fork_length", "origin")

Ots_21_23 %<>%
  mutate(sample_type = substr(sample, 4, 5))

write.csv(Ots_21_23, quote = FALSE, row.names = FALSE, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/Coastal chinook/Coastal Chinook quality filtered sample information DD and KCO genotyping efforts.csv")
