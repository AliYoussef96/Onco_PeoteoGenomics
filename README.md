# Onco_PeoteoGenomics

## The Mutation script

**The [mutations script](https://github.com/AliYoussef96/Onco_PeoteoGenomics/tree/main/mutations) consists of three files;**

**Mutations.functions.r** contains the necessary functions to run the second script **Run.mutation.r.** The final file named **Mutations.r** is only to run the main function, import the packages, and source the rest of the functions.

**How does it work?**

The input needed from the user is only the path to the mutations table file, and the code will run automatically;

1. Read the mutations table from the path given by the user.

2. Using protein IDs, the protein sequence will be downloaded from the Ensemble database.

3. Using the mutation IDs (for e.g. L100G), the downloaded protein sequence will be checked using the normal amino acid position (for e.g. L at 100 or not).

4. If the correct protein was downloaded, the mutation is induced using the mutation ID (for e.g. L will be converted to G in the 100 aa positions).

5. A final check is done to see if the protein is mutated at the correct position or not.

6. If yes, the wild protein sequence and the mutated one will be written in different fasta files with the same header.

```R
library(httr)
library(jsonlite)
library(xml2)
library(stringr)
library(seqinr)

source("Mutations.functions.r")
source("Run.mutation.r")
# Run 

Run.mutation("COSMIC-haematopoietic-point-mutations.csv")

```
