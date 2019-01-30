

### YAGA is a tool that performs ABBA/BABA analysis for Orthofinder results. ###

From an Orthofinder_Results folder it will conduct an automated filtering of single copy orthologs for the species provided and then, paired with a KEGG, or GO file (or both) will associate functional enrichment for the different species and characteristics. 

While ABBA/BABA was not relavent for my master's thesis because a lack of introgression (fungi are not expected to hybridize and horizontal gene transfer is quite rare), I release this tool in the hope it will be useful to others. Instead, I was looking for identity between genes of organisms that have desired traits but within the structure of the ABBA/BABA four-taxa tree. 

### Requirements ###

Needs [R](http://archive.linux.duke.edu/cran/), the statistical programming language. Also, if you want to use the "abba" option, then please install "evobiR" in R:

```r
install.packages("evobiR")
```

Please also install [mafft](https://mafft.cbrc.jp/alignment/software/source.html) and [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) with the "abba" option. If you have MacOS and [Homebrew](https://brew.sh/) just use:

```bash
brew install mafft bedtools
```
If on linux apt-get install should work as well.

The main functions are from python. Please have version 2.7+ 
Check the version with:

```bash
python -v
```
Python (2.7+) packages - ete2, argparse, and Biopython. They can be installed by pip. If you have not installed pip please follow insctructions from [here](https://www.makeuseof.com/tag/install-pip-for-python/). 

Then 
```bash
pip install biopython ete2 argparse
```

### Functionality ###

The algorithm for the original ABBA/BABA-like analysis of YAGA is as follows:

First, the file paths of pertinent files are deduced from the Orthfinder_Results folder as of Orthofinder version 2.2.1. It will also work for up to version 2.2.6 but I ignored those because of this bug. [put bug issue here] We create two new directories "YAGA" in the main folder, and "Single_Copy_Gene_Trees" in the "Orthologues" directory. The "YAGA" folder contains output files directly related to the analysis, while "Single_Copy_Gene_Trees" conatins the Gene trees which only have 1 gene per species tested, and all species must be present, essentially analogous to BUSCOs but for this dataset only. 

From the Orthogroups.csv file, we can gather the Orthogroups, the members of that orthogroup and the order in which they are specified in the file. 

Next, the "target.json" is parsed to glean the info for which species are what. Typically ABBA looks like this:

                 /\
                /  \
               /    \
              /\     \
             /  \     \
            /    \     \
           /\     \     \
          /  \     \     \
         /    \     \     \
        /      \     \     \
       /        \     \     \
      /          \     \     \
     Target  Nearest  Trait  Outgroup

The target and the nearest neighbor (by expected phylogeny) are the two closest for most of the trees. This is because they are specifically chosen to be close, and the gene sequence identity should follow this pattern. For genes that are necessary for a specific trait, a rock-inhabiting fungi for example, the "rock" genes of the "target" and "trait" species would be more similar indicating the tree had rearranged to BABA as shown below:

                 /\
                /  \
               /    \
              /\     \
             /  \     \
            /    \     \
           /\     \     \
          /  \     \     \
         /    \     \     \
        /      \     \     \
       /        \     \     \
      /          \     \     \
    Target    Trait  Nearest  Outgroup

If the signal was strong enough, for instance, out of 3000 single copy orthologs identifed, 30 genes were annotated for GO:0033961 out of 32 of that very same GO:0033961 present in the whole genome overall, then we might be able to say it was enriched, as the fisher test would look something like:

                   | Genes in X | Genes NOT in X |
    ---------------|------------|----------------|
    Genes in A     |    30      |         2      |
    Genes NOT in A |    500     |       2468     |

From the species, when we have the "target" we go through all the trees and find the "sisters" to the target species. From the sisters (assuming the tree is resolved, MSA usually needs to be run with Orthofinder), we can then infer which orthogroups are closer if the sister appears. 

So to reiterate, the GO term GO:0033961 was found 32 times in the whole genome and 30 of those instances were in the "enriched" set. 530 genes in total were in the "enriched" set. Notice I say genes, when I refer to a GO term, this is because if the interproscan file is included, YAGA gets an array of all the GO terms in the file (even duplicates) and then filters for unique per gene, so in effect it counts gene occurences linked to the GO term. Then it is correct to say genes, since the GO terms are deduplicated. 

The OUTPUT file is then fed into R, through the Rscript and the fisher.test is run for each line like so:

```R
> x <- matrix(c(30, 2, 100, 2868), byrow = TRUE, 2, 2)
> fisher.test(x)

Fisher's Exact Test for Count Data

data:  x
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  18.65206 641.52001
sample estimates:
odds ratio 
  73.90576
```

This would be an example of enrichment.


### True ABBA/BABA analysis ###

It might be useful to use the pipeline demonstrated in [this tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/). This does true ABBA/BABA, albeit with a custom ".geno" file. It also involves some preprocessing.

True ABBA/BABA would take SNPs and look at the D-Statistics. It would be easist to accept a GFF/GTF file and extract the DNA sequence from there. But, to make YAGA as accessible as possible, one can optionally use only a protein fasta and a refrence genome in fasta format. Exonerate can take the proteins and create alignments. (In progress)

Regardless of how the DNA is come by, the sequences can then be compared to see certain polymorphisms and evidence of introgression. 


### Running YAGA ###

```bash
python yaga.py abba --target /Users/daniel/Desktop/YAGA/target.json
```

The only required arugment for the "abba" option is "--target" or "-t". This file contains the names of taxa in the analysis, or at least a substring, so the taxa can be recognized in other files. It also needs the gff tracks with structural annotation (read: CDS locations). These are then extracted with bedtools "getfasta" and then joined to a mature DNA which is then aligned with MAFFT to a combination of sequences from the other 3 members of the four-taxa orthogroup. 

The aligned fastas are then read into R and processed by "evobiR", and the D-statistics returned for each single combination and the entire orthogroup average. Of course, this depends entirely on the combination 

Likewise the "yaga" option needs the "--target" option and the GO file (interproscan output) and an optional KEGG annotation. 

```bash
python yaga.py yaga --target /Users/daniel/Desktop/YAGA/target.json -g /Users/daniel/Desktop/protein_fasta.csv -k /Users/daniel/Desktop/kegg.csv
```

### In Progress ###

* Be able to align MAUVE co-linear blocks
* Get duplicates and missing/expanded ("yaga" option)
* Implement "geno" file and genome wide analysis
* Find average D-stats by gene present, not just overall combination

### Tips ###

Orthofinder names the "species" according to the protein fasta file names, so the "target.json" will need to have values similar to the filenames. If not, just run Orthofinder again with renamed protein fasta files. 

Tested on Orthofinder 2.2.1, other Orthofinder versions should work.




