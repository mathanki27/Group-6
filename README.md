![Deakin Image](https://github.com/mathanki27/Group-6/blob/e1e6959b57337d9f8c6478ec521e62b2a74955a3/deakin%20image.png?raw=true)

Group-6
====================================================

The **Group-6** repository is created by students from Deakin University as part of an academic assignment. This project aims to explore various aspects of data analysis and visualization in the field of bioinformatics. Through this repository, we aim to demonstrate our understanding of gene expression analysis, growth data evaluation, and the application of statistical methods. We hope that this project serves as a valuable resource for anyone interested in these topics.

Details
-------
First created through urgency and adrenaline by Aaron Quinlan Spring 2009. 
Maintained by the Quinlan Laboratory at the University of Virginia.

1. **Lead Contributors**:         Mathanki Mehra, Anakha Prasad, Kalyla Dsouza
2. **Repository**:                https://github.com/mathanki27/Group-6.git
3. **License**:                   Released under MIT license

Part 1: Importing files, data wrangling, mathematical operations, plots and saving code on GitHub
====================================================
In our bioinformatics assignment 4, we were assigned to perform an array of operations using R markdown, prepare a report, and collaborate with other group team members through the GitHub repository. The variety of tasks was divided into two parts

Description
--------------
Part 1 – This involved an extended analysis of the RNA seq_count data and circumference measurement which was done with the help of R markdown. The RNAseq dataset was downloaded and deeply analysed using different code chunks to estimate the genes with the highest mean expression, determine the number of genes with a mean less than 10, and make a histogram for the same. The other data which was the growth data was also downloaded and studied which included estimating the mean and standard deviation of tree circumferences, calculating the mean growth, generating a boxplot, and estimation of p-value with the help of t.test. Collaborating through Git Hub has been effective and helpful in engaging in productive teamwork and analyzing the datasets easily.

Downloading Files
--------------
Step 1: Download Files from GitHub
To download the required files from GitHub, use the following R code:
To download a file from GitHub to R, the download.file() function can be employed. The following steps outline the procedure:

1. **Find the Raw File URL**:
   - Navigate to the desired file within the GitHub repository.
   - Click on the file and then select the "Raw" button to obtain the direct link to the file.

2. **Use download.file() in R**:
   - Utilize the download.file() function, specifying the URL of the raw file as the first argument and the destination path on the local machine as the second argument.

```sh
# URL of the first raw file on GitHub
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"

# Path to save the file locally (in the current directory)
destination <- "gene_expression.tsv"

# Download the file
download.file(url, destination)

# URL of the second raw file on GitHub
url1 <- "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"

# Path to save the file locally (in the current directory)
destination <- "growth_data.csv"

# Download the file
download.file(url1, destination)

```
### Step 1
To read the gene_expression.tsv file into R and set the gene identifiers as the row names, the read.table() function can be used. Subsequently, the first six genes can be displayed using the head() function. The following code illustrates these steps:
```sh
# Read in the gene expression data and making gene identifiers as row names
gene_expression <- read.table("gene_expression.tsv", header = TRUE, row.names = 1, sep = "\t")

# Show a table of values for the first six genes
head(gene_expression)
```

**Code Explanation:**

- header = TRUE: This tells R that the first row of the file contains the column names (i.e., headers for each column). This allows R to correctly assign the names to the data frame's columns.

- row.names = 1: This argument specifies that the first column of the data file (index 1) should be treated as row names instead of a regular data column. This means that the values in the first column will not be included as part of the data frame; instead, they will serve as labels for the rows.

- sep = "\t": This indicates that the data is tab-separated. Since the file is a .tsv (tab-separated values), this ensures that R reads the file correctly by treating each tab as a separator between columns.

By using the row.names = 1 argument, the gene identifiers from the first column are set as the row names of the gene_expression data frame. As a result, when we display the data frame using head(gene_expression), we will see the gene identifiers as the row labels.

We can also see the structure of the data frame after reading it, using str() command.
This shows the column names and the row names, allowing us to confirm that the gene identifiers are indeed set as row names. 

### Step 2
The following code is used to calculate the mean of columns as well as add a new column "mean_count" to the existing data frame. 

```sh
# Read in the gene expression data
gene_expression <- read.table("gene_expression.tsv", header = TRUE, row.names = 1, sep = "\t")

# Calculate the mean of the counts across the columns for each gene
gene_expression$mean_count <- rowMeans(gene_expression)

# Show a table of values for the first six genes, including the new mean column
head(gene_expression)
```

**Code Explanation:**

- rowMeans(gene_expression): This function calculates the mean of the counts for each row (gene) across all the specified columns. The result is a vector containing the mean counts for each gene.

- gene_expression$mean_count <-: This assigns the calculated mean counts to a new column called mean_count in the gene_expression data frame.

- head(gene_expression): This displays the first six genes along with all their corresponding count data and the new mean column.

We can further use the View() command to visualize the data frame. By using this command we can further see the newly added column.

### Step 3
To identify the 10 genes with the highest mean expression from the gene_expression data frame, the order() function can be utilized to sort the data based on the mean_count column in descending order. Subsequently, the top 10 rows can be extracted. The following steps outline this process:
```sh
# Read in the gene expression data
gene_expression <- read.table("gene_expression.tsv", header = TRUE, row.names = 1, sep = "\t")

# Calculate the mean of the counts across the columns for each gene
gene_expression$mean_count <- rowMeans(gene_expression)

# Get the 10 genes with the highest mean expression
top_10_genes <- gene_expression[order(-gene_expression$mean_count), ][1:10, ]

# Display the top 10 genes
top_10_genes

```
**Code Explanation:**

- order(-gene_expression$mean_count): This orders the indices of gene_expression based on the mean_count column in descending order (the negative sign - indicates descending order).

- gene_expression[...]: This uses the ordered indices to subset the gene_expression data frame, sorting the rows based on mean expression.

- ([1:10, ]): This selects the first 10 rows of the sorted data frame, which correspond to the genes with the highest mean expression.

- top_10_genes: This stores the resulting data frame containing the top 10 genes, and the final line displays it.

### Step 4
The next step is to estimate the number of genes with a mean that is less than 10. This can be done with the help of the codes given below:

```sh
# Count the number of genes with a mean expression < 10
num_genes_below_10 <- sum(gene_expression$mean_count < 10)

# Display the result
num_genes_below_10
```
**Code Explanation:**
To determine the number of genes that have mean values less than 10, we see the sum total of the categories that is done with the help of the code given above, which is 35988.
- gene_expression$mean_count < 10: This creates a logical vector where each gene's mean expression is checked to see if it's less than 10.
- sum(...): This counts the number of TRUE values in the logical vector, which corresponds to the number of genes with a mean expression less than 10.

### Step 5
This step involves the construction of the histogram which is a bar graph representation of the mean expression values .


```sh
# Create a histogram of the mean expression values
hist(gene_expression$mean_count, 
     main = "Histogram of Mean Gene Expression Values", 
     xlab = "Mean Expression Value", 
     ylab = "Frequency", 
     col = "orange", 
     border = "black")
```

**Code Explanation:**
The code given above is used to generate a histogram of the mean expression values where in the code the function "hist" is used to produce a histogram for the specified data , that is, gene_expression. In the code the main heading for the graph and the x and y axis headings are also specified.
- hist(): This function generates the histogram.
- main =: Sets the title of the plot.
- xlab = and ylab =: Label the x-axis and y-axis.
- col = "orange": Sets the color of the bars.
- border = "black": Adds a black border to the bars for better visualization.

### Step 6
```sh
#Import the csv file
growth_data <- read.csv("growth_data.csv")
```

**Code Explanation:**

The code "read.csv" given above is used to import the csv file named growth_data.csv and converts it into R object , where it is saved as growth_data.

"<-" function is used to direct the output of the read.csv function to the growth_data.

```sh
#Display the column names
colnames(growth_data)
```

**Code Explanation:**

The code given above is used to list the column names and that is done with the help of the function "colnames" after which we also specify the name of the data frame ,that is growth_data.

### Step 7
This step calculates the mean and standard deviation of tree circumference at the start and end end of the study at both sites. 
```sh
# Install new packages

install.packages("dplyr")

install.packages("tidyr")

library("dplyr")

# Ensure the circumference columns are numeric
growth_data$Circumf_2005_cm <- as.numeric(growth_data$Circumf_2005_cm)
growth_data$Circumf_2020_cm <- as.numeric(growth_data$Circumf_2020_cm)
# Calculate mean and standard deviation for both years at each site
summary_stats <- growth_data %>%
  group_by(Site) %>%
  summarise(
    Mean_Circumf_2005 = mean(Circumf_2005_cm, na.rm = TRUE),
    SD_Circumf_2005 = sd(Circumf_2005_cm, na.rm = TRUE),
    Mean_Circumf_2020 = mean(Circumf_2020_cm, na.rm = TRUE),
    SD_Circumf_2020 = sd(Circumf_2020_cm, na.rm = TRUE)
  )
# Print the summary statistics
print(summary_stats)

```

**Code Explanation:**
This command installs the dplyr package in R.
dplyr is a popular R package used for data manipulation. It provides a range of functions to efficiently filter, select, mutate, arrange, and summarize data.
Installing the package ensures that its functions are available to be used in the R environment.

The next command installs the tidyr package in R.
tidyr is an R package designed to help tidy data, making it easier to work with. It offers functions to reshape messy data into a more structured format (e.g., wide or long formats), which is critical for analysis.
Installing tidyr ensures that tools like gather(), spread(), pivot_longer(), and pivot_wider() are available for use in your code.

- Loading the dplyr Library: The code first loads the dplyr library, which is necessary for data manipulation. dplyr provides functions like group_by(), summarise(), and the pipe operator (%>%), all of which are used in this code to handle the dataset.
- Ensuring Columns are Numeric: The columns Circumf_2005_cm and Circumf_2020_cm are converted to numeric using the as.numeric() function. This step ensures that these columns can be used for calculations such as the mean and standard deviation. If these columns were initially non-numeric (e.g., due to text or missing values), this conversion handles that issue.
- The code uses dplyr's group_by() and summarise() functions to calculate summary statistics for each site. It groups the data by the Site column and then calculates the following for each site:

1. The mean circumference for 2005.
2. The standard deviation of the circumferences for 2005.
3. The mean circumference for 2020.
4. The standard deviation of the circumferences for 2020.

- Finally, the calculated summary statistics (mean and standard deviation for both years) are printed. The result is a table that provides insight into the tree growth patterns for each site in both 2005 and 2020.

### Step 8
This step involves the making of a box plot of tree circumference at the start and end of the study. The following code is used to make the box plot:
```sh
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Reshape the data into a long format for easier plotting
long_growth_data <- growth_data %>%
  select(Site, Circumf_2005_cm, Circumf_2020_cm) %>%
  pivot_longer(cols = c(Circumf_2005_cm, Circumf_2020_cm), 
               names_to = "Year", 
               values_to = "Circumference") %>%
  mutate(Year = ifelse(Year == "Circumf_2005_cm", "2005", "2020"))

# Create the box plot
ggplot(long_growth_data, aes(x = Year, y = Circumference, fill = Site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference at Start (2005) and End (2020) of Study",
       x = "Year",
       y = "Tree Circumference (cm)") +
  theme_minimal()
```

**Code Explanation:**
The code first loads three libraries:

1. ggplot2: Used for creating visualizations in R, including box plots.
2. dplyr: A data manipulation library used for filtering, selecting, and reshaping data.
3. tidyr: A library used to reshape data, particularly for tasks like converting between wide and long formats.

The data is reshaped into a long format to make it easier to plot multiple years in a single graph. This involves:

- Selecting the Site, Circumf_2005_cm, and Circumf_2020_cm columns from the growth_data dataset.
- Using the pivot_longer() function to convert the data from wide format (separate columns for 2005 and 2020 circumferences) to long format (a single Year column indicating whether the data is from 2005 or 2020, and a Circumference column for the actual measurements).
- The mutate() function is used to change the Year values from "Circumf_2005_cm" and "Circumf_2020_cm" to "2005" and "2020", making the years easier to understand in the plot.

A box plot is generated using the ggplot2 library. The ggplot() function takes long_growth_data as input and specifies:

- aes(x = Year, y = Circumference, fill = Site): This sets the aesthetics (variables for the plot). The x-axis will represent the Year, the y-axis will show the Circumference, and the Site will be used to differentiate colors for different sites.
- geom_boxplot(): Adds the actual box plots for each year, separated by site.
- labs(): Sets the title of the plot as well as the labels for the x and y axes.
- theme_minimal(): Applies a minimal theme to the plot, giving it a clean appearance with less clutter.

This code reshapes the dataset and creates a box plot to compare tree circumferences in 2005 and 2020, grouped by site. The box plot visually represents the distribution of tree circumferences for each year, allowing for easy comparison between the start and end of the study across different sites.

### Step 9
This step involves the calculation of mean growth over the last 10 years at each site. The following code calculates how much the trees at each site have grown between 2010 and 2020, then computes the average growth for each site. The output is a summary of the mean growth per site over the specified 10-year period.

```sh
# Calculate the growth over the last 10 years (2020 - 2010)
growth_data$Growth_2010_2020 <- growth_data$Circumf_2020_cm - growth_data$Circumf_2010_cm

# Calculate the mean growth for each site
mean_growth <- growth_data %>%
  group_by(Site) %>%
  summarise(Mean_Growth_2010_2020 = mean(Growth_2010_2020, na.rm = TRUE))

# View the result
print(mean_growth)
```

**Code Explanation:**
The code calculates the tree growth over the last 10 years (from 2010 to 2020) for each tree. This is done by subtracting the circumference measured in 2010 (Circumf_2010_cm) from the circumference measured in 2020 (Circumf_2020_cm), creating a new column called Growth_2010_2020 in the growth_data dataset. This new column represents the growth (in cm) over that period.
The code then calculates the mean growth for each site by:

- Grouping the data by Site using the group_by() function. This ensures that the following calculations are done separately for each site.
- Using the summarise() function to calculate the mean of the Growth_2010_2020 column for each group (site). The na.rm = TRUE argument is used to exclude any missing values (NA) from the calculation.

This step creates a new data frame, mean_growth, where each row corresponds to a site, and the column Mean_Growth_2010_2020 contains the average growth for that site.
Finally, the print() function is used to display the mean_growth data frame, which contains the mean growth for each site over the period 2010 to 2020.

### Step 10
This step involves the estimation of p-value using t.test to show that the growth is different at the two sites. The code below is used to test if there is a statistically significant difference in the growth of trees between two different sites over the period from 2010 to 2020. The t-test provides insights into whether the observed difference in growth is likely due to random variation or reflects a true difference between the sites.

```sh
# Step 1: Perform a t-test comparing the growth between the two sites
t_test_result <- t.test(Growth_2010_2020 ~ Site, data = growth_data)

# Step 2: View the p-value and other test results
print(t_test_result)
```
**Code Explanation:**

The code performs a t-test to compare the growth from 2010 to 2020 (Growth_2010_2020) between two different sites in the growth_data dataset. The t-test is used to determine if there is a significant difference in growth between the two sites.

- t.test: This function in R performs a t-test, which is used to test if there is a significant difference between the means of two groups.
- Growth_2010_2020 ~ Site:
  1. Growth_2010_2020 is the numeric variable representing the growth in tree circumference over the last 10 years (2010 to 2020).
  2. Site is a factor variable that indicates which site each growth observation comes from. In this case, you have two sites: "northeast" and "southwest."
  3. The formula Growth_2010_2020 ~ Site tells R to compare the Growth_2010_2020 variable between the two levels of the Site factor (northeast vs southwest).
-	data = growth_data: This specifies that both Growth_2010_2020 and Site are columns in the growth_data data frame. It tells the t.test function to look for these variables in that data frame.

**Welch Two Sample t-test:**
Welch's t-test is used here because it does not assume equal variances between the two groups, which makes it more robust in cases where the variances might differ.

The output of the Welch Two Sample t-test shows the following key points:
-	t-value: **1.8882**
-	Degrees of freedom (df): **87.978**
-	p-value: **0.06229**
Interpretation:
-	P-value: The p-value of 0.06229 is slightly above the common significance threshold of 0.05. This means that there is not enough statistical evidence at the 5% level to reject the null hypothesis that the mean 10-year growth is the same at both sites. However, the p-value is relatively close to 0.05, suggesting a potential trend towards a difference in growth between the two sites, but not at a statistically significant level.

-	Mean growth:
 1. Northeast site: 42.94 cm
 2. Southwest site: 35.49 cm

The northeast site appears to have a slightly higher average growth over the 10 years compared to the southwest site, though the difference is not statistically significant at the 5% level.

Part 2: Examining biological sequence diversity
====================================================

Description
--------------
Part 2- the second part of the assignment focussed on comparing the DNA sequences of two different organisms namely Escherichia coli and Bifidobacteriaceace bacterium. The first step was aimed at downloading both the sequences and estimating the coding sequences which had to be presented in tabular form. Following that we obtained the total coding DNA which was presented in tabular form. We also estimated the length of the coding sequence that was used to generate the boxplot, following which we determined the mean and median values of the coding sequences. Subsequently, we calculated the frequency of bases and protein sequence which was used to generate the bar plots. Additionally, we generated the codon usage table and quantified the codon usage bias. Finally, we identified the over and under-represented k-mers in the protein sequences. This task has been immensely helpful in comparative genomic studies of both organisms.

### Step 1

To download whole set of E. coli. coding DNA sequences, we use the following code:

```sh
# Load the necessary package

library("R.utils")

# Define the URL for the Escherichia coli K-12 substrain MG1655 CDS file

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"

# Download the CDS file from the URL and save it locally as 'ecoli_cds.fa.gz'
download.file(URL,destfile="ecoli_cds.fa.gz")

# Decompress the gzipped file to extract the fasta file (.fa)
gunzip("ecoli_cds.fa.gz")

# List all files in the current working directory to confirm successful extraction
list.files()
```

To download whole set of Bifidobacteriaceae bacterium coding DNA sequences, we use the following code:

```sh
# Load the necessary package for file decompression
library("R.utils")

# URL for Bifidobacteriaceae bacterium CDS

URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_63_collection/bifidobacteriaceae_bacterium_wp012_gca_003585625/cds/Bifidobacteriaceae_bacterium_wp012_gca_003585625.ASM358562v1.cds.all.fa.gz"

# Download the CDS file
download.file(URL,destfile="bifido_cds.fa.gz")

# Decompress the file
gunzip("bifido_cds.fa.gz")

# List the files in the current directory
list.files()

```

To see the number of coding sequences in the organisms, we use the following code:

```sh
# Load necessary library
library("seqinr")

# Reading of the uncompressed files for both the bacteria
cds <- seqinr::read.fasta("ecoli_cds.fa")
cds1 <- seqinr::read.fasta("bifido_cds.fa")

# Shows the coding sequences in E.coli and Bifidobacteriaceae
length(cds)
length(cds1)
```

A more refined code can be used to display the number of sequences in the downloaded coding sequence of both the organisms.

```sh
# Load the necessary library
library(seqinr)

# Function to count the number of coding sequences in a fasta file using seqinr
count_CDS <- function(file) {
  # Read the fasta file
  fasta_content <- read.fasta(file, seqtype = "DNA")
  
  # The length of fasta_content gives the number of sequences (CDSs)
  num_cds <- length(fasta_content)
  
  return(num_cds)
}

# Count the number of CDS in E. coli 
ecoli_cds_count <- count_CDS("ecoli_cds.fa")
cat("Number of coding sequences in Escherichia coli:", ecoli_cds_count, "\n")

# Count the number of CDS in Bifidobacteriaceae bacterium
bifido_cds_count <- count_CDS("bifido_cds.fa")
cat("Number of coding sequences in Bifidobacteriaceae bacterium:", bifido_cds_count, "\n")

```
**Code Explanation:**

- The seqinr library provides tools for biological sequence analysis, such as reading, writing, and manipulating sequence data in various formats, including FASTA files.
- count_CDS: This function takes a FASTA file as an argument and returns the number of coding sequences (CDS) present in that file.
- read.fasta(file, seqtype = "DNA"): Reads the FASTA file specified by file. It is instructed to treat the sequences as DNA (seqtype = "DNA").
- fasta_content holds all the sequences from the file.
- length(fasta_content): This counts the number of sequences in the fasta_content list, which is equivalent to the number of CDS in the file.
- The function returns the number of sequences (CDS) found in the file.
- count_CDS("ecoli_cds.fa"): Calls the count_CDS function and passes the E. coli CDS FASTA file (ecoli_cds.fa) as the input. The result is stored in ecoli_cds_count.
- cat(): Prints a message to the console displaying the number of coding sequences in the E. coli genome.
- Similar to the previous step, this section reads the FASTA file for Bifidobacteriaceae bacterium (bifido_cds.fa), counts the CDS, and prints the result.

To create a table for the number of coding sequences, we will download a new package called "dplyr" and we use the following code:
```sh
# Download the package
install.packages("dplyr")

# Load necessary library for data manipulation
library(dplyr)

# Create a table to compare the number of CDSs in both organisms
cds_table <- data.frame(
  Organism = c("Escherichia coli", "Bifidobacteriaceae bacterium"),
  CDS_Count = c(ecoli_cds_count, bifido_cds_count)
)

# Display the table
print(cds_table)
```

**Code Explanation:**

- install.packages("dplyr"): This command installs the dplyr package, which is a powerful tool for data manipulation in R. You only need to run this once unless the package is not already installed or needs to be updated.
- library(dplyr): Loads the dplyr library into the R session so you can use its functions. It’s a popular package for data manipulation, offering tools like filtering, selecting, grouping, and summarizing data.
- data.frame(): This function creates a data frame, which is a two-dimensional, tabular data structure in R (similar to a spreadsheet or SQL table).
 
1. Organism = c("Escherichia coli", "Bifidobacteriaceae bacterium"): This creates the first column in the data frame called Organism, with the names of two organisms: Escherichia coli and Bifidobacteriaceae bacterium.
2. CDS_Count = c(ecoli_cds_count, bifido_cds_count): This creates the second column called CDS_Count, which stores the number of coding sequences (CDS) in each organism. The values come from the previously calculated variables ecoli_cds_count and bifido_cds_count, which hold the CDS counts for Escherichia coli and Bifidobacteriaceae bacterium, respectively.

- print(cds_table): Displays the contents of the cds_table data frame in the console. This allows you to see a side-by-side comparison of the number of CDSs in the two organisms.

The difference between the 2 organisms is clear in terms of number of sequences. E. coli. has more number of coding sequences than Bifidobacteriaceae bacterium. This can also be seen using the following code:
```sh
# Describe the difference between the two organisms
if (ecoli_cds_count > bifido_cds_count) {
  cat("Escherichia coli has more coding sequences (CDSs) than Bifidobacteriaceae bacterium.\n")


} else if (ecoli_cds_count < bifido_cds_count) {
  cat("Bifidobacteriaceae bacterium has more coding sequences (CDSs) than Escherichia coli.\n")
} else {
  cat("Both organisms have the same number of coding sequences (CDSs).\n")
}
```
### Step 2
Calculation of amount of coding DNA present in total for both organisms and presentation of the data in the form of a tabular column.

```sh
# Load the necessary library
library(seqinr)

# Function to calculate the total length of coding DNA in a fasta file
total_coding_DNA <- function(file) {
  # Read the fasta file
  fasta_content <- read.fasta(file, seqtype = "DNA")
  
  # Calculate the length of each sequence and sum them up
  total_dna <- sum(sapply(fasta_content, length))
  
  return(total_dna)
}

# Calculate the total coding DNA in Escherichia coli
ecoli_total_dna <- total_coding_DNA("ecoli_cds.fa")
cat("Total coding DNA in Escherichia coli:", ecoli_total_dna, "base pairs\n")

# Calculate the total coding DNA in Bifidobacteriaceae bacterium
bifido_total_dna <- total_coding_DNA("bifido_cds.fa")
cat("Total coding DNA in Bifidobacteriaceae bacterium:", bifido_total_dna, "base pairs\n")

```

**Code Explanation:**

The code begins by loading the seqinr library, which provides tools for handling biological sequences, such as reading and manipulating FASTA files (a standard format for representing nucleotide or protein sequences).

The function total_coding_DNA() is created to calculate the total length of coding DNA sequences in a given FASTA file.

Input: The function takes a file name (in this case, a FASTA file) as its argument.
- Step 1: It reads the contents of the FASTA file using read.fasta(), which loads the DNA sequences into R as a list, with each entry representing a sequence. The seqtype = "DNA" argument specifies that the sequences in the file are DNA sequences.
- Step 2: It calculates the length of each sequence in the FASTA file by applying the length() function to each sequence in the list using sapply(). This creates a list of lengths.
- Step 3: The function sums up the lengths of all sequences using sum(), yielding the total length of coding DNA in base pairs.
- Output: The total length of coding DNA is returned by the function.

The function total_coding_DNA() is called with the argument "ecoli_cds.fa", which represents the FASTA file containing the coding DNA sequences of Escherichia coli. The result, stored in the variable ecoli_total_dna, is the total length of coding DNA for E. coli. This result is printed using the cat() function, showing the total number of base pairs.

Similarly, the function total_coding_DNA() is called with the argument "bifido_cds.fa", representing the FASTA file for Bifidobacteriaceae bacterium. The total length of coding DNA for this organism is stored in bifido_total_dna and printed with cat().

The code below shows, how to create a table for the above obtained data.

```sh
# Load necessary library for data manipulation
library(dplyr)

# Create a table to compare the total coding DNA in both organisms
dna_table <- data.frame(
  Organism = c("Escherichia coli", "Bifidobacteriaceae bacterium"),
  Total_Coding_DNA = c(ecoli_total_dna, bifido_total_dna)
)

# Display the table
print(dna_table)

```

**Code Explanation:**
The code starts by loading the dplyr library, which is commonly used for data manipulation tasks such as creating and modifying tables.
The data.frame() function is used to create a table that compares the total coding DNA in two organisms: Escherichia coli and Bifidobacteriaceae bacterium. The table contains two columns:

- Organism: A vector that lists the names of the two organisms ("Escherichia coli" and "Bifidobacteriaceae bacterium").
- Total_Coding_DNA: A vector that holds the previously calculated total coding DNA lengths for these organisms (ecoli_total_dna for Escherichia coli and bifido_total_dna for Bifidobacteriaceae bacterium).

This creates a data frame (dna_table) with one row for each organism and a corresponding value for its total coding DNA.

The print() function is used to display the dna_table, showing the names of the organisms alongside their respective total coding DNA values.

The code creates a simple comparison table showing the total coding DNA for two different organisms and prints it to the console for easy visualization. It helps compare the lengths of coding DNA between Escherichia coli and Bifidobacteriaceae bacterium.

```sh
# Describe the difference in total coding DNA between the two organisms
if (ecoli_total_dna > bifido_total_dna) {
  cat("Escherichia coli has more coding DNA than Bifidobacteriaceae bacterium.\n")
} else if (ecoli_total_dna < bifido_total_dna) {
  cat("Bifidobacteriaceae bacterium has more coding DNA than Escherichia coli.\n")
} else {
  cat("Both organisms have the same amount of coding DNA.\n")
}
```

**Code Explanation:**
The code uses a conditional statement (an if-else structure) to compare the total coding DNA lengths of Escherichia coli and Bifidobacteriaceae bacterium and provide a descriptive output based on the comparison.

This code provides a clear and concise comparison of the total coding DNA lengths between two organisms. Depending on the values of ecoli_total_dna and bifido_total_dna, it outputs an appropriate message to describe which organism has more coding DNA or if they are equal. This helps in understanding the genetic differences between the two organisms based on their coding DNA content.

### Step 3
The first part of the question involves the estimation of the length of the coding sequence in E.coli.
```sh
#Loading the seqinr package
library(seqinr)

##Read the fasta file of E.coli
cds <- seqinr::read.fasta("ecoli_cds.fa")
```
**code explanation**
library(seqinr)- this code is used to load the seqinr package which is used in R as it contains various tools which helps in examining and extracting the genomic sequences of diverse organisms.

cds <- seqinr::read.fasta("ecoli_cds.fa")-this function is used to convert and organises the sequences in the fasta file into R and save it into the variable "cds" . read.fasta is a sequinr tool which is used for reading the sequences in fasta file to R.


```sh
#Estimate the length of coding sequence in bifidobacterium
library(seqinr)

#Read the fasta file of Bifidobacterium 
cds1 <- seqinr::read.fasta("bifido_cds.fa")

```
**code explanation**
The code chunk above is used to install the seqinr package which consists of various tools like read.fasta which is used to convert the fasta file into R so that the sequences stored in the file can be studied and analysed.

The code below defines a function to calculate the lengths of coding sequences (CDS) from FASTA files for two organisms: *Escherichia coli and Bifidobacteriaceae bacterium*.

```sh
# Function to get the lengths of all coding sequences in a fasta file
get_CDS_lengths <- function(file) {
  # Read the fasta file
  fasta_content <- read.fasta(file, seqtype = "DNA")
  
  # Calculate the length of each sequence
  cds_lengths <- sapply(fasta_content, length)
  
  return(cds_lengths)
}

# Get CDS lengths for Escherichia coli
ecoli_cds_lengths <- get_CDS_lengths("ecoli_cds.fa")
cat("Escherichia coli CDS lengths calculated.\n")

# Get CDS lengths for Bifidobacteriaceae bacterium
bifido_cds_lengths <- get_CDS_lengths("bifido_cds.fa")
cat("Bifidobacteriaceae bacterium CDS lengths calculated.\n")

```

**Code Explanation:**
- Function Definition:
The get_CDS_lengths function takes one parameter, file, which represents the path to a FASTA file.
read.fasta(file, seqtype = "DNA") reads the content of the specified FASTA file, treating the sequences as DNA.
sapply(fasta_content, length) applies the length function to each sequence in the FASTA file, returning a vector of lengths corresponding to each coding sequence.
The function returns the vector of CDS lengths.

- Calculating CDS Lengths for Escherichia coli:

ecoli_cds_lengths <- get_CDS_lengths("ecoli_cds.fa") calls the function with the FASTA file for E. coli and stores the resulting lengths in ecoli_cds_lengths.
cat("Escherichia coli CDS lengths calculated.\n") prints a message indicating that the lengths have been calculated.

- Calculating CDS Lengths for Bifidobacteriaceae bacterium:

bifido_cds_lengths <- get_CDS_lengths("bifido_cds.fa") calls the function with the FASTA file for Bifidobacteriaceae and stores the lengths in bifido_cds_lengths.
cat("Bifidobacteriaceae bacterium CDS lengths calculated.\n") prints a message indicating the completion of this calculation.


The next part of step three is to create boxplot for the coding sequence length in both the organisms.

```sh
# Combine the lengths into a single data frame for plotting
cds_length_data <- data.frame(
  Length = c(ecoli_cds_lengths, bifido_cds_lengths),
  Organism = c(rep("Escherichia coli", length(ecoli_cds_lengths)), rep("Bifidobacteriaceae bacterium", length(bifido_cds_lengths)))
)

# Load the necessary library for plotting
library(ggplot2)

# Create a boxplot of CDS lengths for both organisms
ggplot(cds_length_data, aes(x = Organism, y = Length)) +
  geom_boxplot(fill = c("skyblue", "lightgreen")) +
  labs(title = "Boxplot of Coding Sequence Lengths in Escherichia coli and Bifidobacteriaceae bacterium",
       x = "Organism", y = "CDS Length (base pairs)") +
  theme_minimal()

```
**Code Explanation:**
1. Combining Lengths into a Data Frame:

A data frame named cds_length_data is created with two columns:
- Length: This column contains all CDS lengths from both organisms, combining ecoli_cds_lengths and bifido_cds_lengths.
- Organism: This column specifies the organism for each length using the rep function to repeat the organism names for each length. It assigns "Escherichia coli" for lengths in ecoli_cds_lengths and "Bifidobacteriaceae bacterium" for lengths in bifido_cds_lengths.

2. Loading Required Libraries:

library(ggplot2) loads the ggplot2 library, which is used for creating the plot.
Creating the Boxplot:

- ggplot(cds_length_data, aes(x = Organism, y = Length)) initializes a ggplot object using cds_length_data, with the x-axis representing the organism and the y-axis representing the CDS lengths.
- geom_boxplot(fill = c("skyblue", "lightgreen")) adds boxplots to the graph, with each organism colored differently (sky blue for E. coli and light green for Bifidobacteriaceae).
labs(title = "Boxplot of Coding Sequence Lengths in Escherichia coli and Bifidobacteriaceae bacterium", x = "Organism", y = "CDS Length (base pairs)") sets the title and axis labels for the plot.
- theme_minimal() applies a minimal theme to the plot for a cleaner appearance.


```sh

#Summaries the first column as vector
len <- as.numeric(summary(cds)[,1])

#Find the mean of the coding sequence length for E.coli
mean(len)

#Find the median of the coding sequence length for E.coli
median(len)
```
**Code Explanation**
The code chunk used above is used to find the mean and median. First the fasta file is converted to R and then the first column is then summarised as a numerial vector by the name "len" which is then used to obtain the mean and median of the coding sequence length for E.coli, thus the functions mean and median are used respectively.

```sh
#Summarise the first column as vector
len <- as.numeric(summary(cds1)[,1])

#Find the mean of the coding sequence length for Bifidobacterium
mean(len)

#Find the median of the coding sequence for Bifidobacterium
median(len)
```
**Code Explanation**
The code chunk above is used to obtain the mean and median values for the coding sequence lengths of the bifidobacterium.

The next part of the question is to find differences between the two organisms:

```sh
# Describe the differences in mean and median lengths
if (ecoli_mean_length > bifido_mean_length) {
  cat("Escherichia coli has a higher mean coding sequence length than Bifidobacteriaceae bacterium.\n")
} else if (ecoli_mean_length < bifido_mean_length) {
  cat("Bifidobacteriaceae bacterium has a higher mean coding sequence length than Escherichia coli.\n")
} else {
  cat("Both organisms have the same mean coding sequence length.\n")
}

if (ecoli_median_length > bifido_median_length) {
  cat("Escherichia coli has a higher median coding sequence length than Bifidobacteriaceae bacterium.\n")
} else if (ecoli_median_length < bifido_median_length) {
  cat("Bifidobacteriaceae bacterium has a higher median coding sequence length than Escherichia coli.\n")
} else {
  cat("Both organisms have the same median coding sequence length.\n")
}

```
**Code Explanation:**
1. Mean Length Comparison:

The code first checks if the mean coding sequence length of Escherichia coli is greater than that of Bifidobacteriaceae bacterium.
- If true, it prints a message indicating that Escherichia coli has a higher mean length.
- If false, it checks if the mean length of Bifidobacteriaceae bacterium is greater and prints the corresponding message.
- If both means are equal, it indicates that they have the same mean length.

2. Median Length Comparison:

Next, the code performs a similar comparison for the median coding sequence lengths.
- It checks if the median length of Escherichia coli is greater than that of Bifidobacteriaceae bacterium and prints the appropriate message based on the result.
- It also checks if Bifidobacteriaceae bacterium has a greater median length, printing the corresponding message.
- If both medians are equal, it states that they have the same median length.

**The differences in coding sequence (CDS) length between Escherichia coli and Bifidobacteriaceae can be attributed to several factors:**

1. Genomic Organization:

E. coli has a more complex genomic architecture, which can lead to longer coding sequences. Its genome is organized in a manner that often contains larger genes or genes that are part of operons, allowing for longer CDS lengths.
Bifidobacteriaceae, on the other hand, tends to have a more streamlined genome, with a focus on essential metabolic pathways, potentially resulting in shorter coding sequences.

2. Functional Adaptation:

The functions and adaptations of each organism may also play a role. E. coli is a versatile bacterium capable of thriving in various environments, which may require more complex proteins and longer coding sequences to support its metabolic versatility.
In contrast, Bifidobacteriaceae, often found in the gut microbiota, might have evolved to rely on shorter coding sequences for specific functions tailored to their niche environment.

3. Evolutionary Pressure:

Different evolutionary pressures can influence genome size and gene length. E. coli is subject to a wide range of environmental conditions, necessitating a broader set of genes and potentially longer CDS to adapt.
Bifidobacteriaceae may experience different selective pressures, focusing more on efficiency and specialization.

**Importance of Determining Mean and Median Coding Sequences in Bacteria**

- Understanding Genome Structure: Analyzing mean and median coding sequence lengths can provide insights into the overall structure and organization of bacterial genomes, helping to identify evolutionary patterns and functional adaptations.

- Comparative Genomics: Mean and median values allow for the comparison of different bacterial species, shedding light on their evolutionary relationships and functional diversity.

- Functional Annotation: Knowing the distribution of coding sequence lengths can aid in the annotation of genomes, helping researchers predict gene function and understand metabolic capabilities.

- Bioinformatics Applications: In bioinformatics, the mean and median CDS lengths are useful for developing algorithms and models that rely on gene size for various analyses, such as gene prediction or functional annotation.

- Health and Disease Implications: Understanding the coding sequences of gut bacteria like *Bifidobacteriaceae* can have implications for human health, as their metabolic functions can influence gut health and microbiome composition.

In summary, the differences in coding sequence lengths between *E. coli and Bifidobacteriaceae* highlight their evolutionary adaptations and functional roles, while analyzing these metrics is crucial for various applications in genomics and bioinformatics.

### Step 4
Calculation of the frequency of DNA bases and protein in the total coding sequence and creation of a bar plot for the same.

```sh
# Load the required libraries
library(seqinr)
library(ggplot2)

```

**Code Explanation**
seqinr package is used for sequenece manipulation.The unlist function is used to combine the elements of bifido_cds.fa and Ecoli_cds.fa.The c function combines the vectors from both the sequences

```sh
# Function to calculate base frequencies
calculate_base_frequencies <- function(coding_sequences) {
  # Concatenate all coding sequences into a single string
  total_sequence <- paste(coding_sequences, collapse = "")
  
  # Calculate frequencies of each base
  base_counts <- table(unlist(strsplit(total_sequence, "")))
  
  # Convert counts to numeric and calculate frequencies
  base_frequencies <- as.data.frame(base_counts)
  base_frequencies$Frequency <- as.numeric(base_frequencies$Freq) / sum(base_frequencies$Freq)
  names(base_frequencies) <- c("Base", "Count", "Frequency")  # Renaming columns
  
  return(base_frequencies)
}

# Calculate frequencies for E. coli
base_frequencies_ecoli <- calculate_base_frequencies(cds[[2]])

# Calculate frequencies for Bifidobacteriaceae
base_frequencies_bifido <- calculate_base_frequencies(cds1[[2]])

```
**Code Explanation:**
1. Function Definition:
The function calculate_base_frequencies takes a list of coding sequences as input.
It concatenates all coding sequences into a single string using paste.

2. Base Frequency Calculation:
The combined string is split into individual bases (nucleotides) using strsplit.
The frequency of each base (A, T, C, G) is calculated using table, which counts the occurrences of each base.
These counts are converted to a data frame, and the frequency of each base is computed by dividing the count of each base by the total number of bases.

3. Column Renaming:
The columns of the resulting data frame are renamed to "Base," "Count," and "Frequency" for clarity.

4. Calculating Frequencies:
The function is then called to calculate base frequencies for Escherichia coli and Bifidobacteriaceae using the coding sequences stored in cds[[2]] and cds1[[2]], respectively.

This code effectively summarizes the composition of nucleotide bases within the provided coding sequences, allowing for further analysis of their genetic characteristics.

The protein sequences are extracted using the following code:

```sh
# Extract protein sequences from coding sequences for E. coli
prots <- translate(cds[[2]])

# Extract protein sequences from coding sequences for Bifidobacteriaceae
prots1 <- translate(cds1[[2]])

```
**Code Explanation**
The code translates cds and cds1 (contains the coding sequences for E. coli and Bifidobacteriaceace) and extracts protein sequences from these coding sequences using the translate function.
The resulting protein sequences are stored in the variable prots and prots1.

The following code calculates amino acid frequencies from protein sequences for Escherichia coli and Bifidobacteriaceae:
```sh
# Function to calculate amino acid frequencies
calculate_aa_frequencies <- function(protein_sequences) {
  # Concatenate all protein sequences into a single string
  total_protein <- paste(protein_sequences, collapse = "")
  
  # Calculate frequencies of each amino acid
  aa_counts <- table(unlist(strsplit(total_protein, "")))
  
  # Convert counts to numeric and calculate frequencies
  aa_frequencies <- as.data.frame(aa_counts)
  aa_frequencies$Frequency <- as.numeric(aa_frequencies$Freq) / sum(aa_frequencies$Freq)
  names(aa_frequencies) <- c("AminoAcid", "Count", "Frequency")  # Renaming columns
  
  return(aa_frequencies)
}

# Calculate frequencies for E. coli
aa_frequencies_ecoli <- calculate_aa_frequencies(prots)

# Calculate frequencies for Bifidobacteriaceae
aa_frequencies_bifido <- calculate_aa_frequencies(prots1)

```

**Code Explanation**
1. Amino Acid Frequency Calculation Function:

The function calculate_aa_frequencies takes a list of protein sequences as input.
It concatenates all protein sequences into a single string using paste().
The function then counts the occurrences of each amino acid using table() and splits the total protein string into individual characters.

2. Frequency Calculation:

The counts of amino acids are converted into a data frame called aa_frequencies, which includes:
AminoAcid: The amino acid type.
Count: The number of occurrences of each amino acid.
Frequency: The proportion of each amino acid in relation to the total number of amino acids.

3. Frequencies for Each Organism:

The function is called with the protein sequences for E. coli (stored in prots), and the resulting frequencies are stored in aa_frequencies_ecoli.
The same function is called with the protein sequences for Bifidobacteriaceae (stored in prots1), and the results are stored in aa_frequencies_bifido.

The following code effectively visualizes the frequency distribution of DNA bases between the two organisms, allowing for easy comparison.
```sh
# Combine base frequencies for plotting
base_frequencies_ecoli$Organism <- "E. coli"
base_frequencies_bifido$Organism <- "Bifidobacteriaceae"
combined_base_frequencies <- rbind(base_frequencies_ecoli, base_frequencies_bifido)

# Plot for DNA base frequencies
ggplot(combined_base_frequencies, aes(x = Base, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "DNA Base Frequency Comparison",
       x = "Base",
       y = "Frequency")

```
**Code Explanation**
1. Combining Base Frequencies:

Two new columns, Organism, are added to the base frequency data frames (base_frequencies_ecoli and base_frequencies_bifido) to label each organism.
The data frames are then combined into one using rbind(), resulting in combined_base_frequencies.

2. Plotting DNA Base Frequencies:

A bar plot is created using ggplot() with the combined data.
The x-axis represents the Base, and the y-axis shows the Frequency of each base.
The fill aesthetic is used to differentiate between the two organisms.
The plot employs a dodge position to show the bars for each organism side by side.
Additional formatting includes a minimal theme and axis labels, with a title "DNA Base Frequency Comparison."

The following code that combines amino acid frequencies from *Escherichia coli and Bifidobacteriaceae* for plotting:

```sh
# Combine amino acid frequencies for plotting
aa_frequencies_ecoli$Organism <- "E. coli"
aa_frequencies_bifido$Organism <- "Bifidobacteriaceae"
combined_aa_frequencies <- rbind(aa_frequencies_ecoli, aa_frequencies_bifido)

# Plot for amino acid frequencies
ggplot(combined_aa_frequencies, aes(x = AminoAcid, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Amino Acid Frequency Comparison",
       x = "Amino Acid",
       y = "Frequency")

```

**Code Explanation:**
1. Combining Amino Acid Frequencies:

New columns, Organism, are added to the amino acid frequency data frames (aa_frequencies_ecoli and aa_frequencies_bifido) to label each organism.
The data frames are combined into one using rbind(), resulting in combined_aa_frequencies.

2. Plotting Amino Acid Frequencies:

A bar plot is created using ggplot() with the combined data.
The x-axis represents the Amino Acid, and the y-axis shows the Frequency of each amino acid.
The fill aesthetic differentiates between the two organisms.
The plot uses a dodge position to show the bars for each organism side by side.
Additional formatting includes a minimal theme and axis labels, with a title "Amino Acid Frequency Comparison."

### Step 5
This questions demands to generate the codon usage table followed by estimating the codon usage bias for all coding sequences.this also requires to highlight the differences between the codon usage bias of both the organisms and proving the same with the help of plots.

First step involves generating the codon usage list.
```sh
# Load necessary packages
library(seqinr)
library(ggplot2)
```
**Code Explanation** 
This code snippet loads essential R packages for biological analysis and visualization:

- library(seqinr): Used for handling and analyzing DNA and protein sequences.
- library(ggplot2): A powerful tool for creating sophisticated visualizations, such as bar plots and scatter plots.

The codon usage for E.coli is calculated as follows:
```sh
uco(cds, index = "freq")
```

The codon usage for Bifidobacteriaceae is calculated as follows:
```sh
uco(cds1, index = "freq")
```
**Note**: Since both the outputs show zero value for codons, we calculate the codon usage of single nucleotides.

The following code shows the nucleotide usage in both the organisms:
```sh
calculate_codon_usage_all <- function(coding_sequences) {
  # Initialize an empty data frame to hold codon counts
  codon_usage <- data.frame(Codon = character(), Count = integer(), Frequency = numeric())
  
  # Loop through each coding sequence
  for (sequence in coding_sequences) {
    # Combine all coding sequences into a single character vector
    all_sequences <- paste(sequence, collapse = "")
    
    # Split the sequence into codons (3 bases each)
    codons <- unlist(strsplit(all_sequences, "(?<=.)(?=.{3})", perl = TRUE))
    
    # Count occurrences of each codon
    codon_counts <- table(codons)
    
    # Convert to data frame
    temp_codon_usage <- as.data.frame(codon_counts)
    names(temp_codon_usage) <- c("Codon", "Count")
    
    # Calculate the frequency of each codon
    temp_codon_usage$Frequency <- temp_codon_usage$Count / sum(temp_codon_usage$Count)
    
    # Combine with the overall codon usage
    codon_usage <- rbind(codon_usage, temp_codon_usage)
  }
  
  # Aggregate counts and frequencies by codon
  codon_usage <- aggregate(Count ~ Codon, data = codon_usage, sum)
  codon_usage$Frequency <- codon_usage$Count / sum(codon_usage$Count)
  
  return(codon_usage)
}

# Calculate codon usage for all genes in E. coli
codon_usage_ecoli <- calculate_codon_usage_all(cds[[2]])

# Calculate codon usage for all genes in Bifidobacteriaceae
codon_usage_bifido <- calculate_codon_usage_all(cds1[[2]])
```

**Code Explanation**
The function calculates codon usage for multiple coding sequences. It initializes an empty data frame to store codon counts. For each coding sequence, it concatenates the sequences, splits them into codons, and counts the occurrences of each codon. These counts are then converted into a data frame, where frequencies are calculated based on the total count of codons. After processing all sequences, the function aggregates the counts and frequencies by codon and returns the final codon usage data.

The codon usage is calculated separately for all genes in *E. coli* and *Bifidobacteriaceae*.
```sh
head(codon_usage_bifido)
head(codon_usage_ecoli)
```

**Code Explanation** 
The head() function displays the first six rows of the codon usage data for *Bifidobacteriaceae* and *E. coli*. This provides insights into the most frequently used codons in their respective coding sequences, allowing for a comparison of codon usage patterns between the two organisms.


```sh
# Add an organism column for plotting
codon_usage_ecoli$Organism <- "E. coli"
codon_usage_bifido$Organism <- "Bifidobacteriaceae"

# Combine both datasets for plotting
combined_codon_usage <- rbind(codon_usage_ecoli, codon_usage_bifido)

# Ensure the Codon column is treated as a factor
combined_codon_usage$Codon <- factor(combined_codon_usage$Codon, levels = unique(combined_codon_usage$Codon))
```

**Code Explanation** 
The provided code snippet adds an “Organism” column to the codon usage datasets for both *E. coli* and *Bifidobacteriaceae*. It combines the two datasets into a single data frame for easier plotting. Additionally, it ensures that the “Codon” column is treated as a factor, maintaining the order of unique codons as levels. This preparation is important for accurate visualization in plots.

```sh
# Plot codon usage
ggplot(combined_codon_usage, aes(x = Codon, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Codon Usage Comparison",
       x = "Codon",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

**Code Explanation**
The code snippet creates a bar plot to compare codon usage between *E. coli* and *Bifidobacteriaceae*. It uses ggplot2 to generate the plot, specifying:

- x-axis: Codon
- y-axis: Frequency
- fill: Organism (to differentiate between E. coli and Bifidobacteriaceae)
The geom_bar function with position = "dodge" displays the bars side by side for each codon. The theme_minimal() function is applied for a clean layout, and the x-axis text is rotated 45 degrees for better readability. The plot also includes labels for the title, x-axis, and y-axis.

It is observed from the graph that the occurrence of A and T is more in Bifidobacteriaceae than in E.coli and GC content is more in E.coli than in Bifidobacteriaceae

### Step 6
This step identifies 10 protein sequence for variable k-mer lengths which are most over and under represented in *E.coli and Bifidobacteriaceae*.

```sh
# Load necessary libraries
library(seqinr)

# Extract protein sequences from CDS for both organisms
prot_ecoli <- sapply(cds, function(x) translate(x))           # *E. coli*
prot_bifido <- sapply(cds1, function(x) translate(x))         # *Bifidobacteriaceae*

# Combine the protein sequences into single strings for easy k-mer counting
prots_ecoli <- unlist(prot_ecoli)
prots_bifido <- unlist(prot_bifido)

# Define the amino acid alphabet
aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

```
**Code Explanation:**
The code starts by loading the seqinr library, which provides tools for working with biological sequences, including DNA and protein translation.
- The sapply() function is applied to the cds object, which contains the coding DNA sequences (CDS) for Escherichia coli.
- The translate() function converts each CDS sequence (DNA) into its corresponding protein sequence (amino acids).
- The result is stored in prot_ecoli, which holds the translated protein sequences of E. coli.
- unlist() is used to flatten the list of protein sequences for both organisms into single sequences for easier handling.
- prots_ecoli stores the combined protein sequences for E. coli.
- prots_bifido stores the combined protein sequences for Bifidobacteriaceae.
- The vector aa defines the standard 20 amino acids, which will be used for further analysis, such as k-mer counting (analyzing short subsequences of amino acids).

The code below counts and calculates the frequencies of k-mers (subsequences of length 3, 4, and 5) in the protein sequences of Escherichia coli and Bifidobacteriaceae. It compares how often certain amino acid combinations appear in these organisms’ proteins, which can be useful in protein analysis, evolutionary studies, or identifying conserved patterns.

```sh
# Count k-mers for *E. coli*
mycounts_ecoli_3mer <- count(prots_ecoli, wordsize=3, alphabet=aa)
mycounts_ecoli_4mer <- count(prots_ecoli, wordsize=4, alphabet=aa)
mycounts_ecoli_5mer <- count(prots_ecoli, wordsize=5, alphabet=aa)

myfreq_ecoli_3mer <- count(prots_ecoli, wordsize=3, alphabet=aa, freq=TRUE)
myfreq_ecoli_4mer <- count(prots_ecoli, wordsize=4, alphabet=aa, freq=TRUE)
myfreq_ecoli_5mer <- count(prots_ecoli, wordsize=5, alphabet=aa, freq=TRUE)

# Count k-mers for *Bifidobacteriaceae*
mycounts_bifido_3mer <- count(prots_bifido, wordsize=3, alphabet=aa)
mycounts_bifido_4mer <- count(prots_bifido, wordsize=4, alphabet=aa)
mycounts_bifido_5mer <- count(prots_bifido, wordsize=5, alphabet=aa)

myfreq_bifido_3mer <- count(prots_bifido, wordsize=3, alphabet=aa, freq=TRUE)
myfreq_bifido_4mer <- count(prots_bifido, wordsize=4, alphabet=aa, freq=TRUE)
myfreq_bifido_5mer <- count(prots_bifido, wordsize=5, alphabet=aa, freq=TRUE)

```
**Code Explanation:**
1. Counting k-mers for Escherichia coli:
- What are k-mers?: k-mers are short subsequences of length k extracted from a sequence. For example, in a protein sequence, 3-mers would consist of subsequences of three amino acids.

**Counting k-mers for E. coli:**
- count(prots_ecoli, wordsize=3, alphabet=aa): This function counts all possible 3-mers (subsequences of length 3) in the combined protein sequence of Escherichia coli (prots_ecoli). 
- The alphabet=aa specifies that the sequences consist of amino acids (defined in the aa vector).
wordsize=4 and wordsize=5 are used similarly to count 4-mers and 5-mers, respectively, in the protein sequence.
- These counts are stored in mycounts_ecoli_3mer, mycounts_ecoli_4mer, and mycounts_ecoli_5mer for k-mers of size 3, 4, and 5.
- count(prots_ecoli, wordsize=3, alphabet=aa, freq=TRUE): This function works similarly to the one used for counting, but with the argument freq=TRUE, it calculates the frequency (relative proportion) of each 3-mer rather than just the raw count.
- Similar functions calculate the frequencies for 4-mers and 5-mers, stored in myfreq_ecoli_3mer, myfreq_ecoli_4mer, and myfreq_ecoli_5mer.

We can visualize the K-mers counts and frequency using the head command, as shown below:
 
```sh
# Visualizing K-mer count and frequency
head(mycounts_bifido_3mer)
head(myfreq_bifido_3mer)
```
The same above command can be used to visualize k-mers of different lengths for both the bacteria.
The following code can be utilized to calculate over and under represented k-mers of different lengths in *Bifidobacteriaceae*. The code identifies and visualizes the most and least frequent (over- and under-represented) amino acid sequences of length 3 (3-mers), 4 (4-mers), and 5 (5-mers) in the protein sequences of *Bifidobacteriaceae*.
```sh
# Over and under-represented 3mers for bifido
over_represented_bifido_3mer <- sort(myfreq_bifido_3mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 3-mers
under_represented_bifido_3mer <- sort(myfreq_bifido_3mer, decreasing = FALSE)[1:10] # Top 10 under-represented 3-mers

# Over and under-represented 4mers for bifido
over_represented_bifido_4mer <- sort(myfreq_bifido_4mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 4-mers
under_represented_bifido_4mer <- sort(myfreq_bifido_4mer, decreasing = FALSE)[1:10] # Top 10 under-represented 4-mers

# Over and under-represented 5mers for bifido
over_represented_bifido_5mer <- sort(myfreq_bifido_5mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 5-mers
under_represented_bifido_5mer <- sort(myfreq_bifido_5mer, decreasing = FALSE)[1:10] # Top 10 under-represented 5-mers

# Visualize the k-mers
View(over_represented_bifido_3mer)
View(under_represented_bifido_3mer)
View(over_represented_bifido_4mer)
View(under_represented_bifido_4mer)
View(over_represented_bifido_5mer)
View(under_represented_bifido_5mer)
```

**Code Explanation:**
**1. Finding Over- and Under-Represented 3-mers for Bifidobacteriaceae:**
- Over-represented 3-mers:
sort(myfreq_bifido_3mer, decreasing = TRUE)[1:10]: This line sorts the frequencies of all 3-mers in Bifidobacteriaceae's protein sequences in descending order (most frequent first). The top 10 most frequent 3-mers are selected and stored in over_represented_bifido_3mer.
- Under-represented 3-mers:
sort(myfreq_bifido_3mer, decreasing = FALSE)[1:10]: This line sorts the 3-mer frequencies in ascending order (least frequent first), and the top 10 least frequent 3-mers are stored in under_represented_bifido_3mer.

**2. Finding Over- and Under-Represented 4-mers for Bifidobacteriaceae:**
- Over-represented 4-mers:
sort(myfreq_bifido_4mer, decreasing = TRUE)[1:10]: This sorts the frequencies of all 4-mers in Bifidobacteriaceae in descending order and selects the top 10 most frequent 4-mers. These are stored in over_represented_bifido_4mer.
- Under-represented 4-mers:
sort(myfreq_bifido_4mer, decreasing = FALSE)[1:10]: This sorts the frequencies of all 4-mers in ascending order and selects the top 10 least frequent 4-mers, stored in under_represented_bifido_4mer.

**3. Finding Over- and Under-Represented 5-mers for Bifidobacteriaceae:**
- Over-represented 5-mers:
sort(myfreq_bifido_5mer, decreasing = TRUE)[1:10]: Similar to the previous steps, this line sorts the 5-mer frequencies in descending order and stores the top 10 most frequent 5-mers in over_represented_bifido_5mer.
- Under-represented 5-mers:
sort(myfreq_bifido_5mer, decreasing = FALSE)[1:10]: This sorts the 5-mers in ascending order and selects the top 10 least frequent 5-mers, stored in under_represented_bifido_5mer.

The View() function is used to open a new window in RStudio to visualize the top over- and under-represented 3-mers, 4-mers, and 5-mers. This allows for easy inspection of these patterns directly in the RStudio environment.

The following code can be utilized to calculate over and under represented k-mers of different lengths in *E.coli*. The code identifies the most and least frequent k-mers (3-mers, 4-mers, and 5-mers) in *E. coli* protein sequences, helping to detect any recurring or rare amino acid sequence patterns.

```sh
# Over and under-represented 3mers for bifido
over_represented_ecoli_3mer <- sort(myfreq_ecoli_3mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 3-mers
under_represented_ecoli_3mer <- sort(myfreq_ecoli_3mer, decreasing = FALSE)[1:10] # Top 10 under-represented 3-mers

# Over and under-represented 4mers for bifido
over_represented_ecoli_4mer <- sort(myfreq_ecoli_4mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 4-mers
under_represented_ecoli_4mer <- sort(myfreq_ecoli_4mer, decreasing = FALSE)[1:10] # Top 10 under-represented 4-mers

# Over and under-represented 5mers for bifido
over_represented_ecoli_5mer <- sort(myfreq_ecoli_5mer, decreasing = TRUE)[1:10]  # Top 10 over-represented 5-mers
under_represented_ecoli_5mer <- sort(myfreq_ecoli_5mer, decreasing = FALSE)[1:10] # Top 10 under-represented 5-mers

# Visualize the k-mers
View(over_represented_ecoli_3mer)
View(under_represented_ecoli_3mer)
View(over_represented_ecoli_4mer)
View(under_represented_ecoli_4mer)
View(over_represented_ecoli_5mer)
View(under_represented_ecoli_5mer)
```
**Code Explanation:**
This code identifies the top 10 most frequent (**over-represented**) and least frequent (**under-represented**) k-mers (3-mers, 4-mers, and 5-mers) in *E. coli* protein sequences. It uses the `sort()` function to rank the k-mers based on their frequencies—descending for over-represented and ascending for under-represented k-mers. The results are stored in separate variables for each k-mer length (3, 4, and 5), and the `View()` function allows for visual inspection of these results in RStudio. This helps in detecting significant amino acid sequence patterns in *E. coli*.


The following code creates data frames for the top over-represented k-mers (3-mers, 4-mers, and 5-mers) from both *E. coli* and *Bifidobacteriaceae*. Each data frame consists of two columns: kmer, which holds the names of the k-mers, and freq, which contains their corresponding frequencies converted to numeric values.
Instead of directly using over_represented_ecoli_3mer for the plots, data frames are created to avoid any column mismatch issues.

```sh
# Creating data frames for 3mers, 4mers and 5mers
plot_data_ecoli_3mer <- data.frame(
  kmer = names(over_represented_ecoli_3mer),
  freq = as.numeric(over_represented_ecoli_3mer)   # Convert to numeric values
)

plot_data_bifido_3mer <- data.frame(
  kmer = names(over_represented_bifido_3mer),
  freq = as.numeric(over_represented_bifido_3mer)  # Convert to numeric values
)

plot_data_ecoli_4mer <- data.frame(
  kmer = names(over_represented_ecoli_4mer),
  freq = as.numeric(over_represented_ecoli_4mer)
)

plot_data_bifido_4mer <- data.frame(
  kmer = names(over_represented_bifido_4mer),
  freq = as.numeric(over_represented_bifido_4mer)
)

plot_data_ecoli_5mer <- data.frame(
  kmer = names(over_represented_ecoli_5mer),
  freq = as.numeric(over_represented_ecoli_5mer)
)

plot_data_bifido_5mer <- data.frame(
  kmer = names(over_represented_bifido_5mer),
  freq = as.numeric(over_represented_bifido_5mer)
)

```

After creating the data frames, the head() function is used to inspect the first few rows of the 3-mer data frames for both organisms, ensuring that the data has been correctly formatted. This step can be repeated for the 4-mers and 5-mers to verify their correctness as well.

```sh
# Visualizing and Checking if data frame is correct
head(plot_data_ecoli_3mer)
head(plot_data_bifido_3mer)

```

The following code creates bar plots to visualize the top over-represented 3-mers in the protein sequences of E. coli and Bifidobacteriaceae.

```sh
# Plot for *E. coli* over-represented 3-mers
ggplot(plot_data_ecoli_3mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 3-mers in E. coli",
       x = "3-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

# Plot for *Bifidobacteriaceae* over-represented 3-mers
ggplot(plot_data_bifido_3mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 3-mers in Bifidobacteriaceae",
       x = "3-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

```
**Code Explanation:**
Plotting for E. coli:

- The ggplot() function initializes the plot using the plot_data_ecoli_3mer data frame, specifying the x-axis as kmer and the y-axis as freq.
- geom_bar(stat = "identity") creates a bar plot with the heights of the bars representing the frequency of each k-mer.
- theme_minimal() applies a clean and minimalistic theme to the plot.
labs() adds titles and labels to the axes, including the plot title "Top 10 Over-represented 3-mers in E. coli."
- theme(axis.text.x = element_text(angle = 45, hjust = 1)) rotates the x-axis labels by 45 degrees for better visibility.

Plotting for Bifidobacteriaceae:

- The same steps are repeated for the Bifidobacteriaceae k-mer data, using the plot_data_bifido_3mer data frame and updating the title to "Top 10 Over-represented 3-mers in Bifidobacteriaceae."

The following code creates bar plots to visualize the top over-represented 4-mers in the protein sequences of E. coli and Bifidobacteriaceae.

```sh
# Plot for *E. coli* over-represented 4-mers
ggplot(plot_data_ecoli_4mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 4-mers in E. coli",
       x = "4-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Plot for *Bifidobacteriaceae* over-represented 4-mers
ggplot(plot_data_bifido_4mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 4-mers in Bifidobacteriaceae",
       x = "4-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

```

The following code creates bar plots to visualize the top over-represented 5-mers in the protein sequences of E. coli and Bifidobacteriaceae.

```sh
# Plot for *E. coli* over-represented 5-mers
ggplot(plot_data_ecoli_5mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 5-mers in E. coli",
       x = "5-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Plot for *Bifidobacteriaceae* over-represented 5-mers
ggplot(plot_data_bifido_5mer, aes(x = kmer, y = freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Over-represented 5-mers in Bifidobacteriaceae",
       x = "5-mer",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

```

The presence of different levels of specific k-mers or codon usage in the genomes of *E. coli* and *Bifidobacteriaceae* can be attributed to several biological, ecological, and evolutionary factors:

1. Genetic Drift and Evolution
- Adaptive Evolution: Each organism has adapted to its specific environment, which can influence codon usage. For example, E. coli is often found in the gut of mammals and can quickly adapt to different nutrient sources, leading to the evolution of different codon preferences.
- Mutation Pressure: Different mutation rates can affect the frequency of particular codons over generations. This can lead to the accumulation of certain codons while others may diminish.

2. Functional Constraints
- Protein Function: The coding sequences may reflect the functional needs of the proteins produced by the organism. Different proteins may have distinct requirements for translation efficiency, folding, and stability, which can influence codon choice.
- Gene Expression Levels: Genes that are expressed at higher levels often show a bias towards codons that correspond to more abundant tRNAs in the organism. This can lead to overrepresentation of certain k-mers in highly expressed genes.

3. tRNA Availability
- tRNA Pool: The availability of different tRNA molecules in the cell can greatly influence codon usage. Organisms have varying pools of tRNA, which can lead to differences in how often certain codons are used based on which tRNAs are most abundant.
- Codon Bias: A bias towards specific codons often correlates with the availability of corresponding tRNA species. If certain tRNAs are more plentiful, their associated codons will likely be used more frequently, leading to their overrepresentation.

4. Ecological Niche and Lifestyle
- Environment: The ecological niche occupied by an organism affects its metabolic pathways and thus its protein-coding sequences. Bifidobacteriaceae, often associated with fermentation processes in the gut, may have evolved different codon preferences to optimize their survival in a competitive microbial environment.
- Symbiosis: The lifestyle of Bifidobacteriaceae, which includes symbiotic relationships with hosts, may select for different metabolic functions and thus different k-mer frequencies compared to the more versatile E. coli.

5. Horizontal Gene Transfer
- Gene Acquisition: E. coli is known for its ability to acquire genes from other bacteria through horizontal gene transfer. This can lead to variability in codon usage depending on the origin of the transferred genes and the codon preferences of donor organisms.

**Conclusion**

In summary, the differences in k-mer or codon usage between E. coli and Bifidobacteriaceae are influenced by evolutionary history, functional constraints, ecological niches, and the availability of tRNAs. Understanding these factors provides insight into how organisms have adapted at the genetic level to thrive in their specific environments.


Citation
--------
1.	Moeckel, C., et al., A survey of k-mer methods and applications in bioinformatics. Computational and Structural Biotechnology Journal, 2024.
2.	Baldi, P. and A.D. Long, A Bayesian framework for the analysis of microarray expression data: regularized t-test and statistical inferences of gene changes. Bioinformatics, 2001. 17(6): p. 509-519.
3.	Wang, H. and Y. Zhao, RBinds: a user-friendly server for RNA binding site prediction. Computational and Structural Biotechnology Journal, 2020. 18: p. 3762-3765.
4.	Foster, Z.S., S. Chamberlain, and N.J. Grünwald, Taxa: An R package implementing data standards and methods for taxonomic data. F1000Research, 2018. 7.
5.	Gaidatzis, D., et al., Supplementary Material: QuasR: Quantification and annotation of short reads in R.
6.	Jones, M.B., et al., The new bioinformatics: integrating ecological data from the gene to the biosphere. Annu. Rev. Ecol. Evol. Syst., 2006. 37(1): p. 519-544.
7.	Rehm, B., Bioinformatic tools for DNA/protein sequence analysis, functional assignment of genes and protein classification. Applied microbiology and biotechnology, 2001. 57: p. 579-592.
8.	Anastassiou, D., Frequency-domain analysis of biomolecular sequences. Bioinformatics, 2000. 16(12): p. 1073-1081.
9.	Vad, V., et al., Generalized box-plot for root growth ensembles. BMC bioinformatics, 2017. 18: p. 1-15.
10.	Chang, C.C.H., et al., Bioinformatics approaches for improved recombinant protein production in Escherichia coli: protein solubility prediction. Briefings in bioinformatics, 2014. 15(6): p. 953-962.
11.	Bottacini, F., et al., Comparative genomics of the genus Bifidobacterium. Microbiology, 2010. 156(11): p. 3243-3254.
12.	Carvajal-Rodríguez, A., Myriads: P-value-based multiple testing correction. Bioinformatics, 2018. 34(6): p. 1043-1045.
13.	Hui, L., Quantitative evaluations of variations using the population mean as a baseline for bioinformatics interpretation. PeerJ, 2023. 11: p. e14955.
14.	Costa, A.M., J.T. Machado, and M.D. Quelhas, Histogram-based DNA analysis for the visualization of chromosome, genome and species information. Bioinformatics, 2011. 27(9): p. 1207-1214.
15.	Ramırez, O.R., The importance and implementation of rich documentation in bioinformatics software projects. 2019.
16.	Zhao, S., et al., Bioinformatics for RNA-seq data analysis. Bioinformatics—Updated Features and Applications: InTech, 2016: p. 125-49.
17.	Breidenbach, J.D., et al., GeneToList: a web application to assist with gene identifiers for the non-bioinformatics-savvy scientist. Biology, 2022. 11(8): p. 1113.
18.	Hampson, S., D. Kibler, and P. Baldi, Distribution patterns of over-represented k-mers in non-coding yeast DNA. Bioinformatics, 2002. 18(4): p. 513-528.
19.	Kawulok, J. and S. Deorowicz, CoMeta: classification of metagenomes using k-mers. PloS one, 2015. 10(4): p. e0121453.
20.	Goldman, N. and Z. Yang, A codon-based model of nucleotide substitution for protein-coding DNA sequences. Molecular biology and evolution, 1994. 11(5): p. 725-736.



