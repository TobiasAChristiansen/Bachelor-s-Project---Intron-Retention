###############################################################################
#                                Initialization
###############################################################################

#install.packages("tidyverse")
#install.packages("ggsignif")
#install.packages("patchwork")
library("tidyverse")
library("ggsignif")
library("patchwork")

#For accessibility reasons, the color palette will be changed to Okabe-Ito
#This is a color palette that is colorblind friendly
palette.colors(palette = "Okabe-Ito")























###############################################################################
#                                Loading data
###############################################################################

#Reading the data from the tsv file
data <- paste0(getwd(), "/Bachelorprojekt/Proteomics_results/result_HLAsignificantMatches_revised.txt") |>
    read_delim(col_names = TRUE, delim = "\t")


##Loading metadata from BREG22_ages.txt
#moreinfo <- paste0(getwd(), "/Bachelorprojekt/Proteomics_results/BREG22_ages.txt") |>
#    read_delim(col_names = TRUE, delim = "\t") |>
#    select(ID, Sex)
#
##Joining the to data frames, thereby adding Diagnosis, Age and Sex information to the data
#data <-
#    data |>
#    inner_join(moreinfo, join_by(ID == ID))


#Removing potential NA's
data <-
    data |>
    filter(!is.na(all()))

#Because of the charges in the proteomics data,
#some the observations in the same intron variants
#in the same brain regions of the same patients are duplicate
#Therefore, we're using "distinct" to remove these duplicates
data <-
    data |>
    distinct()

#Printing the 5 first lines of the data frame
data |>
    print(n = 5)






























###############################################################################
#                                Data wrangling
###############################################################################


#Counting 9-mers for every brain region for every patient
Counts_9mer_BREG <-
    data |>
    select(ID, Sex, Age, Brain_region, Diagnosis, Peptide_fragment) |>
    distinct() |>
    group_by(ID, Sex, Age, Brain_region, Diagnosis) |>
    summarise(count = n()) |>
    ungroup()


#Counting total 9-mers (Not per brain region)
Counts_9mer <-
    data |>
    select(ID, Sex, Age, Diagnosis, Peptide_fragment) |>
    distinct() |>
    group_by(ID, Sex, Age, Diagnosis) |>
    summarise(count = n()) |>
    ungroup()


#Counting intron retention variants (found in Gene_name column) for every brain region for every patient
Counts_IR_var_BREG <-
    data |>
    select(ID, Sex, Age, Brain_region, Diagnosis, Gene_name) |>
    distinct() |>
    group_by(ID, Sex, Age, Brain_region, Diagnosis) |>
    summarise(count = n()) |>
    ungroup()


#Counting total 9-mers (Not per brain region)
Counts_IR_var <-
    data |>
    select(ID, Sex, Age, Diagnosis, Gene_name) |>
    distinct() |>
    group_by(ID, Sex, Age, Diagnosis) |>
    summarise(count = n()) |>
    ungroup()


#Wide data containing 1s and 0s for absence or prescence of HLA-relevant IR variants
wide_IR_vars <- 
    data |>                                                             #Starting point
    distinct(ID, Diagnosis, Gene_name) |>                               #Getting unique observations
    mutate(Expression = 1) |>                                           #Changing the expression column to 1
    mutate(across(ID, as.character)) |>                                 #Converting IDs to strings
    pivot_wider(names_from = Gene_name, values_from = Expression) |>    #Pivoting wider
    mutate_if(is.double, ~replace(., is.na(.), 0))                      #Changing NA to 0


#Wide data containing 1s and 0s for absence or prescence of HLA-relevant 9-mers
wide_9mer <- 
    data |>                                                             #Starting point
    distinct(ID, Diagnosis, Peptide_fragment) |>                               #Getting unique observations
    mutate(Expression = 1) |>                                           #Changing the expression column to 1
    mutate(across(ID, as.character)) |>                                 #Converting IDs to strings
    pivot_wider(names_from = Peptide_fragment, values_from = Expression) |>    #Pivoting wider
    mutate_if(is.double, ~replace(., is.na(.), 0))                      #Changing NA to 0

























###############################################################################
#                                Statistics
###############################################################################


#Number of patients in each patient group
data |>
    select(ID, Diagnosis) |>
    distinct() |>
    group_by(Diagnosis) |>
    summarise(count = n())



















#Printing the counts of the IR_variants in each disease group in each brain region
data |>                                                 #Starting point
    select(ID, Diagnosis, Brain_region, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(Diagnosis, Brain_region) |>    #Grouping by Diagnosis and Brain_region
    summarise(count = n()) |>                           #Counting the number of observations
    pivot_wider(names_from = Brain_region, values_from = count) |> #Pivoting wider
    print()                                             #Printing the result


#Printing the average number of IR_variants in each disease group in each brain region
data |>                                                 #Starting point
    select(ID, Diagnosis, Brain_region, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Brain_region) |>    #Grouping by Diagnosis and Brain_region
    summarise(count = n()) |>                           #Counting the number of observations
    group_by(Diagnosis, Brain_region) |>
    summarise(mean = mean(count), sd = sd(count)) |>
    mutate(mean_and_sd = paste0(signif(mean, 3), "+-", signif(sd, 3))) |>
    select(Diagnosis, Brain_region, mean_and_sd) |>
    pivot_wider(names_from = Brain_region, values_from = mean_and_sd) |> #Pivoting wider
    print()                                             #Printing the result


#Printing the counts of the IR_variants in each disease group across all brain regions
data |>                                                                 #Starting point
    select(ID, Diagnosis, Gene_name) |>                                 #Selecting the columns we want to use
    distinct() |>                                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>                                          #Grouping by Diagnosis and ID
    summarise(count = n()) |>                                           #Counting the number of observations
    group_by(Diagnosis) |>                                              #Grouping only by disease group
    summarise(total = sum(count), mean = mean(count), sd = sd(count)) |>     #Calculating mean and standard deviation
    print()                                                             #Printing the result




#Mean and standard deviation of number of intron retention variants. Results are by Diagnosis and Sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, Gene_name) |>            #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>                     #Grouping by ID, Diagnosis, and Sex
    summarise(IR_var_count = n()) |>                     #Counting intron retention variants
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>                         #Grouping by Diagnosis and Sex
    summarise(K_mer_mean = mean(IR_var_count), K_mer_sd = sd(IR_var_count)) |> #Calculating the mean of the numeric columns
    mutate(mean_and_sd = paste0(signif(K_mer_mean, digits = 3), "+-", signif(K_mer_sd, digits = 3))) |>
    select(Diagnosis, Sex, mean_and_sd) |>
    pivot_wider(names_from = Sex, values_from = mean_and_sd) |>
    print()




















#Printing the counts of the K_mers in each disease group across all brain regions
data |>                                                                 #Starting point
    select(ID, Diagnosis, Peptide_fragment) |>                                 #Selecting the columns we want to use
#    distinct() |>                                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>                                          #Grouping by Diagnosis and ID
    summarise(count = n()) |>                                           #Counting the number of observations
    group_by(Diagnosis) |>                                              #Grouping only by disease group
    summarise(total = sum(count), mean = mean(count), sd = sd(count)) |>     #Calculating mean and standard deviation
    print()                                                             #Printing the result




#Mean and standard deviation of number of HLA-relevant 9mers. Results are by Diagnosis and Sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, Peptide_fragment) |>            #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>                     #Grouping by ID, Diagnosis, and Sex
    summarise(k_mer_count = n()) |>                     #Counting intron retention variants
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>                         #Grouping by Diagnosis and Sex
    summarise(K_mer_mean = mean(k_mer_count), K_mer_sd = sd(k_mer_count)) |> #Calculating the mean of the numeric columns
    mutate(mean_and_sd = paste0(signif(K_mer_mean, digits = 3), "+-", signif(K_mer_sd, digits = 3))) |>
    select(Diagnosis, Sex, mean_and_sd) |>
    pivot_wider(names_from = Sex, values_from = mean_and_sd) |>
    print()



#Mean and standard deviation of number of HLA-relevant 9mers. Results are by Diagnosis
data |>                                                 #Starting point
    select(ID, Diagnosis, Peptide_fragment) |>            #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>                     #Grouping by ID, Diagnosis, and Sex
    summarise(k_mer_count = n()) |>                     #Counting intron retention variants
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis) |>                         #Grouping by Diagnosis and Sex
    summarise(K_mer_mean = mean(k_mer_count), K_mer_sd = sd(k_mer_count)) |> #Calculating the mean of the numeric columns
    mutate(mean_and_sd = paste0(signif(K_mer_mean, digits = 3), "+-", signif(K_mer_sd, digits = 3))) |>
    select(Diagnosis, mean_and_sd) |>
    print()




























#Finding the frequencies of patients having intron retention variants across diagnoses
frequencies_IR_var <-
    wide_IR_vars |>                                                 #Starting point
    select(-ID) |>
    group_by(Diagnosis) |>    #Grouping by ID and Diagnosis
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_longer(!Diagnosis, names_to = "Intron_retention_variant", values_to = "Frequency")
frequencies_IR_var



#Finding the frequencies of patients having specific HLA-relevant 9-mers across diagnoses
frequencies_9mer <-
    wide_9mer |>                                                 #Starting point
    select(-ID) |>
    group_by(Diagnosis) |>    #Grouping by ID and Diagnosis
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_longer(!Diagnosis, names_to = "Peptide_fragment", values_to = "Frequency")
frequencies_9mer









###############################################################################
#                             Plotting - Boxplots
###############################################################################





#HLA-relevant 9-mers across diagnoses and brain regions
Counts_9mer_BREG |>
    ggplot(aes(y = count, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape=Sex), size = 5) +
    facet_wrap(~Brain_region) +
    labs(title = "HLA-relevant 9-mers in different brain regions",
         x = "Disease group",
         y = "Number of HLA-relevant 9-mers") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PD"),
                    c("CTRL", "PSP"),
                    c("MSA", "PD"),
                    c("MSA", "PSP"),
                    c("PD", "PSP")),
                step_increase = 0.1,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20))

Counts_9mer_BREG |>
    ggplot(aes(x = count, fill = Diagnosis)) +
    geom_histogram(bins = 8) +
    facet_wrap(~Diagnosis)













#HLA-relevant IR variants across diagnoses and brain regions
Counts_IR_var_BREG |>
    ggplot(aes(y = count, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape=Sex), size = 5) +
    facet_wrap(~Brain_region) +
    ylim(0,120) +
    labs(title = "HLA-relevant IR variants across brain regions",
         x = "Diagnosis",
         y = "HLA-relevant IR variants") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*3 <= 0.05, " *", "")), p*3),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PSP"),
                    c("MSA", "PSP")),
                step_increase = 0.1,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15), 
          axis.text.y=element_text(size=15),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20))




















#Distribution of the ages of the patient groups
Counts_IR_var |>
    select(ID, Sex, Age, Diagnosis) |>
    distinct() |>
    ggplot(aes(y = Age, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Age distribution",
         subtitle = "Boxplot - Stratified by diagnosis",
         x = "Diagnosis",
         y = "Age") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PSP"),
                    c("CTRL", "PD"),
                    c("MSA", "PSP"),
                    c("MSA", "PD"),
                    c("PSP", "PD")),
                step_increase = 0.07,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))





#HLA-relevant matches overall among patients separated by diagnosis
Counts_IR_var |>
    ggplot(aes(y = count, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Distribution of HLA-relevant IR variants",
         subtitle = "Boxplot - stratified by diagnosis",
         x = "Diagnosis",
         y = "Unique HLA-relevant IR variants") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PD"),
                    c("CTRL", "PSP"),
                    c("MSA", "PD"),
                    c("MSA", "PSP"),
                    c("PD", "PSP")),
                step_increase = 0.07,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1.5, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))












#HLA-relevant matches overall among patients separated by diagnosis
Counts_9mer |>
    ggplot(aes(y = count, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "HLA-relevant 9-mers - Unique counts",
         subtitle = "Boxplot - stratified by diagnosis",
         x = "Diagnosis",
         y = "HLA-relevant 9-mers") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PD"),
                    c("CTRL", "PSP"),
                    c("MSA", "PD"),
                    c("MSA", "PSP"),
                    c("PD", "PSP")),
                step_increase = 0.07,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1.5, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))

t.test()












#Counting total 9-mers (Not per brain region)
data |>
    select(ID, Sex, Age, Diagnosis, Gene_name) |>
    group_by(ID, Sex, Age, Diagnosis) |>
    summarise(count = n()) |>
    ungroup() |>
    ggplot(aes(y = count, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "HLA-relevant 9-mers - Non-unique counts",
         subtitle = "Boxplot - stratified by diagnosis",
         x = "Diagnosis",
         y = "Amount of HLA-relevant 9-mers") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PD"),
                    c("CTRL", "PSP"),
                    c("MSA", "PD"),
                    c("MSA", "PSP"),
                    c("PD", "PSP")),
                step_increase = 0.07,
                textsize = 7,
                test = "t.test") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1.5, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))








#Distribution 
































###############################################################################
#                          Plotting - Scatter plots
###############################################################################

#Scatter plot with fitted lines - number of HLA-relevant 9-mers vs. age
Counts_9mer |>
    ggplot(aes(x = Age, y = count, color = Diagnosis)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Number of HLA-relevant 9-mer matches vs. age",
         subtitle = "Linear models, stratified by disease group",
         x = "Age",
         y = "Amount of HLA-relevant 9-mer matches") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold")) +
    geom_smooth(data = Counts_9mer |>
                    filter(Diagnosis == "MSA"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_9mer |>
                    filter(Diagnosis == "PD"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_9mer |>
                    filter(Diagnosis == "PSP"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_9mer |>
                    filter(Diagnosis == "CTRL"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)


Counts_9mer |>
    filter(Diagnosis == "CTRL") |>
    lm(count ~ Age, data = _) |>
    summary()

#Summaries of the linear models
Counts_9mer |>
    filter(Diagnosis == "CTRL") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_9mer |>
    filter(Diagnosis == "MSA") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_9mer |>
    filter(Diagnosis == "PD") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_9mer |>
    filter(Diagnosis == "PSP") |>
    lm(count ~ Age, data = _) |>
    summary()




#Scatter plot with fitted lines - number of HLA-relevant intron retention variants vs. age
Counts_IR_var |>
    ggplot(aes(x = Age, y = count, color = Diagnosis)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "HLA-relevant IR variants vs. age",
         subtitle = "Linear models, stratified by disease group",
         x = "Age",
         y = "HLA-relevant IR variants") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold")) +
    geom_smooth(data = Counts_IR_var |>
                    filter(Diagnosis == "MSA"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_IR_var |>
                    filter(Diagnosis == "PD"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_IR_var |>
                    filter(Diagnosis == "PSP"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Counts_IR_var |>
                    filter(Diagnosis == "CTRL"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)



#Summaries of the linear models
Counts_IR_var |>
    filter(Diagnosis == "CTRL") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_IR_var |>
    filter(Diagnosis == "MSA") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_IR_var |>
    filter(Diagnosis == "PD") |>
    lm(count ~ Age, data = _) |>
    summary()

Counts_IR_var |>
    filter(Diagnosis == "PSP") |>
    lm(count ~ Age, data = _) |>
    summary()





#Number of intron retention variants by sex and patient group
Count_ir_vars_plot_by_age <-
    data |>
    select(ID, Sex, Age, Diagnosis, Gene_name) |>                                      
    distinct() |>
    group_by(ID, Sex, Age, Diagnosis) |>
    summarise(count = n())

Count_ir_vars_plot_by_age |>
    ggplot(aes(x = Age, y = count, color = Diagnosis, lty = Sex)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "HLA-relevant IR variants vs. age",
         subtitle = "Stratified by sex and facet wrapped on patient group",
         x = "Age",
         y = "HLA-relevant IR variants") +
    facet_wrap(~Diagnosis) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20)) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "MSA"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "PD"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "PSP"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "CTRL"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)



















###############################################################################
#                          Plotting - Heat maps
###############################################################################


#Heatmap of the different IR variants
frequencies_IR_var |>
    ggplot(aes(x = Diagnosis, y = Intron_retention_variant, fill = Frequency)) +
    geom_tile() +
    scale_fill_gradient(low = "white",
                         high = "darkred") +
    labs(title = "Patient frequency of specific HLA-relevant IR variants",
         subtitle = "Heatmap of frequencies- stratified by diagnosis",
         x = "Diagnosis",
         y = "HLA-relevant IR variants") +
    theme(axis.text.x=element_text(size=15), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))








#Heatmap of the different HLA-relevant 9-mers
frequencies_9mer |>
    ggplot(aes(x = Diagnosis, y = Peptide_fragment, fill = Frequency)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkred") +
    labs(title = "Patient frequency of specific HLA-relevant 9-mers",
         subtitle = "Heatmap of frequencies- stratified by diagnosis",
         x = "Diagnosis",
         y = "HLA-relevant 9-mers") +
    theme(axis.text.x=element_text(size=15), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))





















































#Creating a wide format data frame from the Diagnosis column
frequencies_IR_var_wide <-
    frequencies_IR_var |>
    pivot_wider(names_from = Diagnosis, values_from = Frequency)



#Finding the lines containing frequencies over 0.5 among disease groups and below 0.5 for controls
frequencies_IR_var_wide_MSA <-
    frequencies_IR_var_wide |>
    filter(MSA > 0.5 & CTRL < 0.5)
frequencies_IR_var_wide_PSP <-
    frequencies_IR_var_wide  |>
    filter(PSP > 0.5 & CTRL < 0.5)
frequencies_IR_var_wide_PD <-
    frequencies_IR_var_wide  |>
    filter(PD > 0.5 & CTRL < 0.5)


#Creating a data frame containing all of the above and sorting it
all_freqs <-
      frequencies_IR_var_wide_MSA |>
      union(frequencies_IR_var_wide_PSP) |>
      union(frequencies_IR_var_wide_PD)

#Sums of disease group frequencies
all_freqs <-
    all_freqs |>
    select(CTRL, MSA, PD, PSP) |>
    apply(1, sum) |>
    mutate(all_freqs, sum = _)

#Max disease group frequencies
all_freqs <-
    all_freqs |>
    select(MSA, PD, PSP) |>
    apply(1, max) |>
    mutate(all_freqs, max = _)

#Sorting by penalizing high sums and high controls
all_freqs |>
    mutate(score = max / sum) |>
    arrange(desc(score)) |>
    select(Intron_retention_variant, CTRL, MSA, PD, PSP, score) |>
    print(n = 50)
















#Creating a wide format data frame from the Diagnosis column
frequencies_kmers_wide <-
    frequencies_9mer |>
    pivot_wider(names_from = Diagnosis, values_from = Frequency)

#Finding the lines containing frequencies over 0.5 among disease groups and below 0.5 for controls
frequencies_kmers_vars_wide_MSA <-
    frequencies_kmers_wide |>
    filter(MSA > 0.5 & CTRL < 0.5)
frequencies_kmers_vars_wide_PSP <-
    frequencies_kmers_wide  |>
    filter(PSP > 0.5 & CTRL < 0.5)
frequencies_kmers_vars_wide_PD <-
    frequencies_kmers_wide  |>
    filter(PD > 0.5 & CTRL < 0.5)


#Creating a data frame containing all of the above and sorting it
all_kmer_freqs <-
      frequencies_kmers_vars_wide_MSA |>
      union(frequencies_kmers_vars_wide_PSP) |>
      union(frequencies_kmers_vars_wide_PD)

#Sums of disease group frequencies
all_kmer_freqs <-
    all_kmer_freqs |>
    select(CTRL, MSA, PD, PSP) |>
    apply(1, sum) |>
    mutate(all_kmer_freqs, sum = _)

#Max disease group frequencies
all_kmer_freqs <-
    all_kmer_freqs |>
    select(MSA, PD, PSP) |>
    apply(1, max) |>
    mutate(all_kmer_freqs, max = _) |>
    remove_missing()

#Sorting by penalizing high sums and high controls
all_kmer_freqs |>
    mutate(score = max / sum) |>
    arrange(desc(score)) |>
    distinct(CTRL, MSA, PD, PSP, sum, max, .keep_all = TRUE) |>
    select(Peptide_fragment, CTRL, MSA, PD, PSP, score) |>
    print(n = 50)























p3_1 <-
    data_peptides_2_1 |>
    filter(Sex == "Male") |>
    ggplot(aes(x = Age, y = count, color = Diagnosis)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Number of HLA-relevant matches vs. age - Males",
         subtitle = "Linear models",
         x = "Age",
         y = "Amount of HLA-relevant matches") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20)) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "MSA" & Sex == "Male"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PD" & Sex == "Male"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PSP" & Sex == "Male"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "CTRL" & Sex == "Male"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)


p3_2 <-
    data_peptides_2_1 |>
    filter(Sex == "Female") |>
    ggplot(aes(x = Age, y = count, color = Diagnosis)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Number of HLA-relevant matches vs. age - Females",
         subtitle = "Linear models",
         x = "Age",
         y = "Amount of HLA-relevant matches") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20)) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "MSA" & Sex == "Female"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PD" & Sex == "Female"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PSP" & Sex == "Female"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "CTRL" & Sex == "Female"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)

p3_1/p3_2



















#HLA-relevant matches vs age, stratified by disease and sex
data_peptides_2_1 |>
    ggplot(aes(x = Age, y = count, color = Diagnosis, lty = Sex)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Number of HLA-relevant matches vs. age - by disease and sex",
         subtitle = "Linear models",
         x = "Age",
         y = "Amount of HLA-relevant matches") +
    facet_wrap(~Diagnosis) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20)) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "MSA"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PD"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "PSP"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = data_peptides_2_1 |>
                    filter(Diagnosis == "CTRL"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)


















#Mean ir variants: non-unique IR variants
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, Gene_name) |>  #Selecting the columns we want to use                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(IR_var_count = n()) |>           #Counting the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_wider(names_from = Sex, values_from = IR_var_count) |>
    print() 





#Mean ir variants: unique IR variants
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(IR_var_count = n()) |>           #Calculating the mean number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_wider(names_from = Sex, values_from = IR_var_count) |>
    print()  



#Mean K_mers non-unique
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, K_mer) |>  #Selecting the columns we want to use
#    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>           #Calculating number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_wider(names_from = Sex, values_from = K_mer_count) |>
    print()  


#Mean K_mers unique
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, K_mer) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    pivot_wider(names_from = Sex, values_from = K_mer_count) |>
    print()  
















#Mean ir variants unique with standard deviation - Not by sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(IR_var_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis) |>
    summarise(IR_var_mean = mean(IR_var_count), IR_var_sd = sd(IR_var_count))
    #summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the of the numeric columns
    print()  




#Mean K_mers unique with standard deviation - Not by sex
data |>                                                 #Starting point
    select(ID, Diagnosis, K_mer) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis) |>
    summarise(K_mer_mean = mean(K_mer_count), K_mer_sd = sd(K_mer_count))
    #summarise(across(where(is.numeric), mean, na.rm = FALSE)) |> #Calculating the mean of the numeric columns
    print()  







#Mean ir variants unique with standard deviation - By sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(IR_var_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(IR_var_mean = mean(IR_var_count), IR_var_sd = sd(IR_var_count)) |> #Calculating the mean of the numeric columns
    mutate(value = paste0(signif(IR_var_mean, digits = 3), " +- ", signif(IR_var_sd, digits = 3))) |>
    select(Diagnosis, Sex, value) |>
    pivot_wider(names_from = Sex, values_from = value) |>
    print()  









#Mean K_mers unique with standard deviation - By sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, K_mer) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(K_mer_mean = mean(K_mer_count), K_mer_sd = sd(K_mer_count)) |> #Calculating the mean of the numeric columns
    mutate(value = paste0(signif(K_mer_mean, digits = 3), " +- ", signif(K_mer_sd, digits = 3))) |>
    select(Diagnosis, Sex, value) |>
    pivot_wider(names_from = Sex, values_from = value) |>
    print()  






























#Count_ir_vars_plot_by_age
Count_ir_vars_plot_by_age <-
    data |>
    select(ID, Sex, Age, Diagnosis, Gene_name) |>                                      
    distinct() |>
    group_by(ID, Sex, Age, Diagnosis) |>
    summarise(count = n())

Count_ir_vars_plot_by_age |>
    ggplot(aes(x = Age, y = count, color = Diagnosis, lty = Sex)) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Intron retention variant counts vs. age - by disease and sex",
         subtitle = "Linear models",
         x = "Age",
         y = "Amount of HLA-relevant matches") +
    facet_wrap(~Diagnosis) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20)) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "MSA"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "PD"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "PSP"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE) +
    geom_smooth(data = Count_ir_vars_plot_by_age |>
                    filter(Diagnosis == "CTRL"), aes(x = Age, y = count), method = "lm", formula = y ~ x, se = FALSE)



Counts_IR_var_BREG |>
    filter(Diagnosis == "CTRL") |>
    lm(count ~ Age, data = _) |>
    summary()











#Mean ir variants unique with standard deviation - By sex
sd_data_irvars_collected <-
    data |>                                                 #Starting point
    select(ID, Diagnosis, Gene_name) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(IR_var_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis) |>
    summarise(IR_var_mean = mean(IR_var_count), IR_var_sd = sd(IR_var_count)) |> #Calculating the mean of the numeric columns
    ungroup()

sd_data_irvars_collected |>
    ggplot() +
    geom_bar(aes(x = Diagnosis, y = IR_var_mean, fill = Diagnosis), stat = "identity") +
    geom_errorbar(aes(x = Diagnosis, ymin = IR_var_mean - IR_var_sd, ymax = IR_var_mean + IR_var_sd), width=0.4, alpha=0.9, size=1.3, color = "black")+
    labs(title = "HLA-relevant intron retention variant means by disease",
         subtitle = "Bar plot - with error bars",
         x = "Diagnosis",
         y = "HLA-relevant intron retention variants") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          strip.text = element_text(size=20))









#Mean K_mers unique with standard deviation - By sex
data |>                                                 #Starting point
    select(ID, Diagnosis, Sex, K_mer) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis, Sex) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>           #Calculating the number of IR_variants
    ungroup() |>
    select(-ID) |>                                      #Removing the ID column
    group_by(Diagnosis, Sex) |>
    summarise(K_mer_mean = mean(K_mer_count), K_mer_sd = sd(K_mer_count)) |> #Calculating the mean of the numeric columns
    mutate(value = paste0(signif(K_mer_mean, digits = 3), " +- ", signif(K_mer_sd, digits = 3))) |>
    select(Diagnosis, Sex, value) |>
    pivot_wider(names_from = Sex, values_from = value) |>
    print()  





data |>                                                 #Starting point
    select(ID, Diagnosis, K_mer) |>  #Selecting the columns we want to use
    distinct() |>                                       #Getting unique observations (Only 1 count per IR_variant per patient)
    group_by(ID, Diagnosis) |>    #Grouping by Diagnosis, Brain_region and IR_variant
    summarise(K_mer_count = n()) |>
    group_by(Diagnosis) |>
    summarise(Sum = sum(K_mer_count), Mean = mean(K_mer_count), sd = sd(K_mer_count))












#Unique IR variants in each patient group in the filtered data
data |>
    select(Diagnosis, Gene_name) |>
    distinct() |>
    group_by(Diagnosis) |>
    summarise(Unique_ir_vars_in_group = n())

#Total unique IR variants in the filtered data
data |>
    select(Gene_name) |>
    distinct() |>
    summarise(total_unique_ir_vars = n())










data |>
    select(ID, Diagnosis, Brain_region, Gene_name) |>
    distinct() |>
    group_by(ID, Diagnosis, Brain_region) |>
    summarise(count_IR_var = n()) |>
    group_by(Diagnosis) |>
    summarise(mean_IR_per_patient_BREG = mean(count_IR_var), sd = sd(count_IR_var))



data |>
    select(ID, Sex, Diagnosis, Brain_region, Peptide_fragment) |>
    distinct() |>
    group_by(ID, Sex, Diagnosis, Brain_region) |>
    summarise(count_IR_var = n()) |>
    ggplot(aes(y = count_IR_var, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "HLA-relevant IR variants per brain region",
         subtitle = "Boxplot - Stratified by diagnosis",
         x = "Diagnosis",
         y = "HLA-relevant IR variants") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PSP"),
                    c("CTRL", "PD"),
                    c("MSA", "PSP"),
                    c("MSA", "PD"),
                    c("PSP", "PD")),
                step_increase = 0.07,
                textsize = 7) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=15), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))







data |>
    select(ID, Sex, Diagnosis, Brain_region) |>
    distinct() |>
    group_by(ID, Sex, Diagnosis) |>
    summarise(num_of_BREG = n()) |>
    ggplot(aes(y = num_of_BREG, x = Diagnosis, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(aes(shape = Sex), size = 5) +
    labs(title = "Number of brain regions for patients in HLA-relevant data",
         subtitle = "Boxplot - Stratified by diagnosis",
         x = "Diagnosis",
         y = "Brain regions") +
    geom_signif(map_signif_level = function(p) sprintf(paste0("adj. p = %.2g", ifelse(p*6 <= 0.05, " *", "")), p*6),
                comparisons = list(
                    c("CTRL", "MSA"),
                    c("CTRL", "PSP"),
                    c("CTRL", "PD"),
                    c("MSA", "PSP"),
                    c("MSA", "PD"),
                    c("PSP", "PD")),
                step_increase = 0.07,
                textsize = 7) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20), 
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.text=element_text(size=20, face = "bold"),
          legend.title = element_text(size=20, face = "bold"),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size=25, face = "bold"),
          plot.subtitle = element_text(size=18),
          axis.title=element_text(size=20,face="bold"))



data |>
    select(Diagnosis, ID, Gene_name) |>
    distinct() |>
    group_by(Diagnosis, ID) |>
    summarise(counts_ir_vars = n()) |>
    group_by(Diagnosis) |>
    summarise(mean_ir_vars = mean(counts_ir_vars),
              sd_ir_vars = sd(counts_ir_vars))



data |>
    select(Diagnosis, Sex, ID, Gene_name) |>
    distinct() |>
    group_by(Diagnosis, Sex, ID) |>
    summarise(counts_ir_vars = n()) |>
    group_by(Diagnosis, Sex) |>
    summarise(mean_ir_vars = mean(counts_ir_vars),
              sd_ir_vars = sd(counts_ir_vars)) |>
    mutate(mean_sd = paste0(signif(mean_ir_vars, 3), "+/-", signif(sd_ir_vars, 3))) |>
    select(Diagnosis, Sex, mean_sd) |>
    pivot_wider(names_from = Sex, values_from = mean_sd)




data |>
    select(Diagnosis, ID, Gene_name) |>
    distinct() |>
    group_by(Diagnosis, ID) |>
    summarise(counts_ir_vars = n()) |>
    group_by(Diagnosis) |>
    summarise(mean_ir_vars = mean(counts_ir_vars)) #|>
    #group_by(Diagnosis) |>
    #summarise(mean_matches = mean(mean_BREG),
              #sd_matches = sd(mean_BREG))




data |>
    select(Diagnosis, ID, Age) |>
    distinct() |>
    ggplot(aes(x = Age, fill = Diagnosis)) +
    geom_histogram(bins = 5) +
    facet_wrap(~Diagnosis)






frequencies_9mer |>
    pivot_wider(names_from = Diagnosis, values_from = Frequency) |>
    filter(MSA >= 0.5 & CTRL >= 0.5 & PD >= 0.5 & PSP >= 0.5)
