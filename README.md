Part 1: Importing files, data wrangling, mathematical operations, plots and saving code on GitHub
====================================================

Description
--------------


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


Citation
--------


