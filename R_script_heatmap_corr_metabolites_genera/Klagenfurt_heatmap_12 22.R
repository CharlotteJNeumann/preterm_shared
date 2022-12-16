#install.packages("ggpubr")
library("ggpubr")

#install.packages("tidyverse")
library("tidyverse")

# Set the working directory containing the files to import the data from and save plots into.
setwd("C:/Users/charl/Nextcloud/Charlotte/Preterm/Paper II/Submissions/Nature Communications/Major_revisions_02/04_heatmap")

## Metabolomic data.
# Import the data
metbabolomic_data <- read_csv("preterm_metabolom_tp7.csv")

# Replace "-" with "_" in the sample IDs.
metbabolomic_data$Sample <- gsub("-", "_", metbabolomic_data$Sample)
# Renname the Sample as Sample_ID
names(metbabolomic_data)[names(metbabolomic_data) == "Sample"] <- "Sample_ID"

# Replace "3,4-" with "" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "3,4-", replacement = "")
# Replace "5-" with "" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "5-", replacement = "")
# Replace "2'FL" with "Three" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "2'FL", replacement = "Two_FL")
# Replace "3FL" with "Three" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "3FL", replacement = "Three_FL")
# Replace "3SL" with "Three_SL" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "3SL", replacement = "Three_SL")
# Replace "6SL" with "Six_SL" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "6SL", replacement = "Six_SL")
# Replace "-" with "_" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = "-", replacement = "_")
# Replace "," with "_" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = ",", replacement = "_")
# Replace " " with "_" in the column headings.
colnames(metbabolomic_data) <- gsub(x = colnames(metbabolomic_data), pattern = " ", replacement = "_")




# Select the metbabolomic data to analyse.
metbabolomic_data <- metbabolomic_data %>%
                                   select(Sample_ID,
                                     Capric_acid:Methyl
                                   )

## Genus microbiota data.
# Import the genus data
genus_data <- read_csv("rel_ab_bac_amplicon_tp7_clr.csv")

# Make the "Features" column into rownames.
genus_data <- column_to_rownames(genus_data, var = "Features")
# Transpose the data.
genus_data_t <- t(genus_data)
# Convert the count data into percentages.
genus_data_t_prop <- genus_data_t/rowSums(genus_data_t)*100
# Make the new rownames into a column named "Sample".
genus_data_t_prop <- rownames_to_column(as.data.frame(genus_data_t_prop), "Sample")

# Import the data
genus_metadata <- read_csv("preterm_metadata_4correlation.csv")
# Rename column headings to match the other file.
names(genus_metadata)[names(genus_metadata) == "sample ID"] <- "Sample"
names(genus_metadata)[names(genus_metadata) == "label"] <- "Sample_ID"
names(genus_metadata)[names(genus_metadata) == "sample no"] <- "Sample_no"
# Replace "-" with "_" in the sample ID.
genus_metadata$Sample_ID <- gsub("-", "_", genus_metadata$Sample_ID)

# Select meta data and filter by hospital.
genus_metadata <- genus_metadata %>%
  select(Sample,
         Sample_ID,
         Sample_no,
         Group) %>%
  filter(Group == "Klagenfurt")


# Merge the genus data and metadata using the "Sample" column.
genus_data_t_merged <- merge(genus_metadata, genus_data_t_prop, by = c("Sample"))

# Select the data to only include the genus data and the Sample ID column.
genus_data_t_merged <- genus_data_t_merged %>%
  select(Sample_ID,
         Bifidobacterium, Lactobacillus, Staphylococcus)

# Merge the metabolomic data and the genus data together using the "Sample_ID" column to select only
# genus data for samples that have matching metabolomic data.
merged_all_data <- merge(metbabolomic_data, genus_data_t_merged, by = c("Sample_ID"))
#here, for Klagenfurt, only 7 samples remain

# Select the genera into a new data frame.
genus_table <- merged_all_data %>%
  select(Sample_ID,
         Bifidobacterium, Lactobacillus, Staphylococcus)
# Make the Sample_ID column into rownames.
genus_table <- column_to_rownames(genus_table, var = "Sample_ID")

# Select the top 10 genus in the dataframe.
# Order the data frame by the abundance of the genera.
genus_table_ordered <- genus_table[,order(colSums(genus_table),decreasing=TRUE)]
# Extract a list of top N Genera names.
N <- 3
taxa_list <- colnames(genus_table_ordered)[1:N]
view(taxa_list)

# Select only the top ten most abundant genus data.
genus_table_top <- data.frame(genus_table_ordered[,colnames(genus_table_ordered) %in% taxa_list])
names(genus_table_top)
# Transform the genus data to make the distribution normal for correlations.

#install.packages("bestNormalize")
library("bestNormalize")
help("bestNormalize")

hist(genus_table_top$Bifidobacterium)
xtrans <- bestNormalize(genus_table_top$Bifidobacterium, loo=TRUE)
hist(xtrans$x.t)
genus_table_top$Bifidobacterium <- xtrans$x.t

hist(genus_table_top$Lactobacillus)
xtrans <- bestNormalize(genus_table_top$Lactobacillus, loo=TRUE)
hist(xtrans$x.t)
genus_table_top$Lactobacillus <- xtrans$x.t

hist(genus_table_top$Staphylococcus)
xtrans <- bestNormalize(genus_table_top$Staphylococcus, loo=TRUE)
hist(xtrans$x.t)
genus_table_top$Staphylococcus <- xtrans$x.t


# Select the metabolites into a final dataframe to use in the correlation analysis.
metabolite_table <- merged_all_data %>%
  select(Sample_ID,
         Capric_acid:Methyl)

# Make the Sample_ID column into rownames.
metabolite_table <- column_to_rownames(metabolite_table, var = "Sample_ID")

# Order the dataframe by the abundance of the genera.
#metabolite_table_ordered <- metabolite_table[,order(colSums(metabolite_table),decreasing=TRUE)]
# Extract a list of top N Genera names.
#N <- 49
#metabolite_list <- colnames(metabolite_table_ordered)[1:N]
# Select only the top ten most abundant genus data.
#metabolite_table_top <- data.frame(metabolite_table_ordered[,colnames(metabolite_table_ordered) %in% metabolite_list])
metabolite_table_top <- metabolite_table


names <- names(metabolite_table_top)
names <- as.data.frame(names)
names

hist(metabolite_table_top$Capric_acid)
xtrans <- bestNormalize(metabolite_table_top$Capric_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Capric_acid <- xtrans$x.t

hist(metabolite_table_top$Butyric_acid)
xtrans <- bestNormalize(metabolite_table_top$Butyric_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Butyric_acid <- xtrans$x.t

hist(metabolite_table_top$L_Leucine)
xtrans <- bestNormalize(metabolite_table_top$L_Leucine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Leucine <- xtrans$x.t

hist(metabolite_table_top$L_Valine)
xtrans <- bestNormalize(metabolite_table_top$L_Valine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Valine <- xtrans$x.t

hist(metabolite_table_top$L_Isoleucine)
xtrans <- bestNormalize(metabolite_table_top$L_Isoleucine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Isoleucine <- xtrans$x.t

hist(metabolite_table_top$Propionic_acid)
xtrans <- bestNormalize(metabolite_table_top$Propionic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Propionic_acid <- xtrans$x.t

hist(metabolite_table_top$Propylene_glycol)
xtrans <- bestNormalize(metabolite_table_top$Propylene_glycol, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Propylene_glycol <- xtrans$x.t

hist(metabolite_table_top$Isopropyl_alcohol)
xtrans <- bestNormalize(metabolite_table_top$Isopropyl_alcohol, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Isopropyl_alcohol <- xtrans$x.t

hist(metabolite_table_top$Valeric_acid)
xtrans <- bestNormalize(metabolite_table_top$Valeric_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Valeric_acid <- xtrans$x.t

hist(metabolite_table_top$L_Lactic_acid)
xtrans <- bestNormalize(metabolite_table_top$L_Lactic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Lactic_acid <- xtrans$x.t

hist(metabolite_table_top$L_Alanine)
xtrans <- bestNormalize(metabolite_table_top$L_Alanine,loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Alanine <- xtrans$x.t

hist(metabolite_table_top$L_Lysine)
xtrans <- bestNormalize(metabolite_table_top$L_Lysine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Lysine <- xtrans$x.t

hist(metabolite_table_top$Acetic_acid)
xtrans <- bestNormalize(metabolite_table_top$Acetic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Acetic_acid <- xtrans$x.t

hist(metabolite_table_top$Aminopentanoic_acid)
xtrans <- bestNormalize(metabolite_table_top$Aminopentanoic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Aminopentanoic_acid <- xtrans$x.t

hist(metabolite_table_top$L_Glutamic_acid)
xtrans <- bestNormalize(metabolite_table_top$L_Glutamic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Glutamic_acid <- xtrans$x.t

hist(metabolite_table_top$Succinic_acid)
xtrans <- bestNormalize(metabolite_table_top$Succinic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Succinic_acid <- xtrans$x.t

hist(metabolite_table_top$beta_Alanine)
xtrans <- bestNormalize(metabolite_table_top$beta_Alanine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$beta_Alanine <- xtrans$x.t

hist(metabolite_table_top$Methylamine)
xtrans <- bestNormalize(metabolite_table_top$Methylamine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Methylamine <- xtrans$x.t

hist(metabolite_table_top$L_Aspartic_acid)
xtrans <- bestNormalize(metabolite_table_top$L_Aspartic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Aspartic_acid <- xtrans$x.t

hist(metabolite_table_top$Sarcosine)
xtrans <- bestNormalize(metabolite_table_top$Sarcosine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Sarcosine <- xtrans$x.t

hist(metabolite_table_top$Sarcosine)
xtrans <- bestNormalize(metabolite_table_top$Sarcosine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Sarcosine <- xtrans$x.t

hist(metabolite_table_top$Trimethylamine)
xtrans <- bestNormalize(metabolite_table_top$Trimethylamine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Trimethylamine <- xtrans$x.t

hist(metabolite_table_top$L_Asparagine)
xtrans <- bestNormalize(metabolite_table_top$L_Asparagine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Asparagine <- xtrans$x.t

hist(metabolite_table_top$Malonic_acid)
xtrans <- bestNormalize(metabolite_table_top$Malonic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Malonic_acid <- xtrans$x.t

hist(metabolite_table_top$Choline)
xtrans <- bestNormalize(metabolite_table_top$Choline, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Choline <- xtrans$x.t

hist(metabolite_table_top$Betaine)
xtrans <- bestNormalize(metabolite_table_top$Betaine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Betaine <- xtrans$x.t

hist(metabolite_table_top$Taurine)
xtrans <- bestNormalize(metabolite_table_top$Taurine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Taurine <- xtrans$x.t

hist(metabolite_table_top$Phenylacetic_acid)
xtrans <- bestNormalize(metabolite_table_top$Phenylacetic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Phenylacetic_acid <- xtrans$x.t

hist(metabolite_table_top$Glycine)
xtrans <- bestNormalize(metabolite_table_top$Glycine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Glycine <- xtrans$x.t

hist(metabolite_table_top$Glycerol)
xtrans <- bestNormalize(metabolite_table_top$Glycerol, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Glycerol <- xtrans$x.t

hist(metabolite_table_top$D_Fructose)
xtrans <- bestNormalize(metabolite_table_top$D_Fructose, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$D_Fructose <- xtrans$x.t

hist(metabolite_table_top$L_Proline)
xtrans <- bestNormalize(metabolite_table_top$L_Proline, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Proline <- xtrans$x.t

hist(metabolite_table_top$Lactulose)
xtrans <- bestNormalize(metabolite_table_top$Lactulose, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Lactulose <- xtrans$x.t

hist(metabolite_table_top$D_Glucose)
xtrans <- bestNormalize(metabolite_table_top$D_Glucose, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$D_Glucose <- xtrans$x.t

hist(metabolite_table_top$Galactose)
xtrans <- bestNormalize(metabolite_table_top$Galactose, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Galactose <- xtrans$x.t

hist(metabolite_table_top$Sucrose)
xtrans <- bestNormalize(metabolite_table_top$Sucrose, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Sucrose <- xtrans$x.t

hist(metabolite_table_top$Uracil)
xtrans <- bestNormalize(metabolite_table_top$Uracil, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Uracil <- xtrans$x.t

hist(metabolite_table_top$Fumaric_acid)
xtrans <- bestNormalize(metabolite_table_top$Fumaric_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Fumaric_acid <- xtrans$x.t

hist(metabolite_table_top$Deoxyuridine)
xtrans <- bestNormalize(metabolite_table_top$Deoxyuridine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Deoxyuridine <- xtrans$x.t

hist(metabolite_table_top$Dihydroxybenzeneacetic_acid)
xtrans <- bestNormalize(metabolite_table_top$Dihydroxybenzeneacetic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Dihydroxybenzeneacetic_acid <- xtrans$x.t

hist(metabolite_table_top$L_Tyrosine)
xtrans <- bestNormalize(metabolite_table_top$L_Tyrosine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Tyrosine <- xtrans$x.t

hist(metabolite_table_top$Gallic_acid)
xtrans <- bestNormalize(metabolite_table_top$Gallic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Gallic_acid <- xtrans$x.t

hist(metabolite_table_top$L_Histidine)
xtrans <- bestNormalize(metabolite_table_top$L_Histidine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Histidine <- xtrans$x.t

hist(metabolite_table_top$p_Cresol)
xtrans <- bestNormalize(metabolite_table_top$p_Cresol, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$p_Cresol <- xtrans$x.t

hist(metabolite_table_top$L_Phenylalanine)
xtrans <- bestNormalize(metabolite_table_top$L_Phenylalanine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Phenylalanine <- xtrans$x.t

hist(metabolite_table_top$L_Tryptophan)
xtrans <- bestNormalize(metabolite_table_top$L_Tryptophan, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$L_Tryptophan <- xtrans$x.t

hist(metabolite_table_top$Xanthine)
xtrans <- bestNormalize(metabolite_table_top$Xanthine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Xanthine <- xtrans$x.t

hist(metabolite_table_top$Hypoxanthine)
xtrans <- bestNormalize(metabolite_table_top$Hypoxanthine, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Hypoxanthine <- xtrans$x.t

hist(metabolite_table_top$Formic_acid)
xtrans <- bestNormalize(metabolite_table_top$Formic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Formic_acid <- xtrans$x.t

hist(metabolite_table_top$Nicotinic_acid)
xtrans <- bestNormalize(metabolite_table_top$Nicotinic_acid, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Nicotinic_acid <- xtrans$x.t

hist(metabolite_table_top$Two_FL)
xtrans <- bestNormalize(metabolite_table_top$Two_FL, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Two_FL <- xtrans$x.t

hist(metabolite_table_top$Three_SL)
xtrans <- bestNormalize(metabolite_table_top$Three_SL, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Three_SL <- xtrans$x.t

hist(metabolite_table_top$Six_SL)
xtrans <- bestNormalize(metabolite_table_top$Six_SL, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Six_SL <- xtrans$x.t

hist(metabolite_table_top$LDNDFH)
xtrans <- bestNormalize(metabolite_table_top$LDNDFH, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LDNDFH <- xtrans$x.t

hist(metabolite_table_top$LNFP1)
xtrans <- bestNormalize(metabolite_table_top$LNFP1, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LNFP1 <- xtrans$x.t

hist(metabolite_table_top$LNFP1)
xtrans <- bestNormalize(metabolite_table_top$LNFP1, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LNFP1 <- xtrans$x.t

hist(metabolite_table_top$LNFP3)
xtrans <- bestNormalize(metabolite_table_top$LNFP3, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LNFP3 <- xtrans$x.t

hist(metabolite_table_top$LNnT)
xtrans <- bestNormalize(metabolite_table_top$LNnT, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LNnT <- xtrans$x.t

hist(metabolite_table_top$LNT)
xtrans <- bestNormalize(metabolite_table_top$LNT, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LNT <- xtrans$x.t

hist(metabolite_table_top$Three_FL)
xtrans <- bestNormalize(metabolite_table_top$Three_FL, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Three_FL <- xtrans$x.t

hist(metabolite_table_top$LDFT)
xtrans <- bestNormalize(metabolite_table_top$LDFT, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LDFT <- xtrans$x.t

hist(metabolite_table_top$LSTa)
xtrans <- bestNormalize(metabolite_table_top$LSTa, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LSTa <- xtrans$x.t

hist(metabolite_table_top$LSTb)
xtrans <- bestNormalize(metabolite_table_top$LSTb, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LSTb <- xtrans$x.t

hist(metabolite_table_top$LSTc)
xtrans <- bestNormalize(metabolite_table_top$LSTc, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$LSTc <- xtrans$x.t

hist(metabolite_table_top$Methyl)
xtrans <- bestNormalize(metabolite_table_top$Methyl, loo=TRUE)
hist(xtrans$x.t)
metabolite_table_top$Methyl <- xtrans$x.t


# Put the two data frames into x and y files to put into the correlation script.
x <- genus_table_top
y <- metabolite_table_top

# Create single grouping variable.
# If not using a grouping variable just leave this set to a created variable of group 1.
study_no_sub <- dplyr::select(merged_all_data,
                              c(Sample_ID))
study_no_sub$group <- study_no_sub$Sample_ID
study_no_sub$group <- 1
# Add grouping data even if group is not used.
groups <- study_no_sub[,2]

#You can use kendall, spearman, or pearson below:
method <- "pearson"
colnames(x)
colnames(y)

#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp
      }
      else{
        df<-rbind(df,tmp)
      }
    }
  }
}
?cor
a
b
i

df <- data.frame(row.names=NULL,df)
df

colnames(df)<-c("Metabolite", "Taxa", "Correlation", "Pvalue", "Group")

df$Pvalue<-as.numeric(as.character(df$Pvalue))

df$AdjPvalue<-rep(0,dim(df)[1])

df$Correlation<-as.numeric(as.character(df$Correlation))

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label <- c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment <- 4

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Metabolite)){
    for(j in unique(df$Type)){
      sel<-df$Metabolite==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Metabolite)){
    sel<-df$Metabolite==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for significant values
df$Significance <- cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
df


ggplot(aes(x = Metabolite,
           y = Taxa,
           fill = Correlation),
       data = df) +
  geom_tile() +
  scale_fill_gradient2(low = "slateblue",
                       mid = "white",
                       high = "red") +
  labs(title="Graz Hospital (Pearson)",
       x="Genus",
       y = "Metabolites") +
  theme(legend.position = "right",
        plot.title = element_text(size=8),
        legend.title = element_text(size=8), #Define the legend title.
        legend.text = element_text(size=8), #Define the size of the legend text.
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), #Define the plot margin size.
        panel.spacing = unit(.25, "lines"),
        legend.background = element_blank(), #Remove the grey background.
        panel.background = element_blank(), #Remove the grey plot background.
        panel.border = element_blank(), #Remove the plot border.
        panel.grid.major = element_blank(), #Remove the major plot gridlines.
        panel.grid.minor = element_blank(), #Remove the minor plot gridlines.
        strip.text.x = element_text(size = 8, colour = "black", angle = 90, vjust = .5, hjust = 0.5),
        axis.title.x = element_text(size=8, colour = "black"), #Define x axis title.
        axis.text.x = element_blank(), #Define x axis labels.
        axis.title.y = element_text(size=8, colour = "black"), #Define y axis title text size.
        axis.text.y = element_text(size=8, colour = "black"), #Define the axis label text size.
        axis.ticks.x = element_blank(),
        aspect.ratio = 30
  ) +
  geom_text(aes(label=Significance),
            color="black",
            size=3)+
  labs(y=NULL, x=NULL, fill=method) +
  facet_grid(. ~ Metabolite,
             drop=TRUE
             #labeller=labeller(Env_labeller)
  )

ggsave("Klagenfurt_heatmap_12_22.svg",
       width = 20,
       height = 25,
       units = c("cm"),
       dpi = 300)

write.csv(df,"Klagenfurt_heatmap_12_22.csv")
