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


Details
-------
First created through urgency and adrenaline by Aaron Quinlan Spring 2009. 
Maintained by the Quinlan Laboratory at the University of Virginia.

1. **Lead Contributors**:         Mathanki Mehra, Anakha Prasad, Kalyla Dsouza
2. **Repository**:                https://github.com/mathanki27/Group-6.git
3. **License**:                   Released under MIT license


Part 2: Examining biological sequence diversity
====================================================

Description
--------------
Part 2- the second part of the assignment focussed on comparing the DNA sequences of two different organisms namely Escherichia coli and Bifidobacteriaceace bacterium. The first step was aimed at downloading both the sequences and estimating the coding sequences which had to be presented in tabular form. Following that we obtained the total coding DNA which was presented in tabular form. We also estimated the length of the coding sequence that was used to generate the boxplot, following which we determined the mean and median values of the coding sequences. Subsequently, we calculated the frequency of bases and protein sequence which was used to generate the bar plots. Additionally, we generated the codon usage table and quantified the codon usage bias. Finally, we identified the over and under-represented k-mers in the protein sequences. This task has been immensely helpful in comparative genomic studies of both organisms.

Citation
--------


