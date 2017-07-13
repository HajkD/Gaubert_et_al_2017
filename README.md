# Reproducible Scripts for the Publication
Gaubert H, Sanchez DH, Drost H-G and Paszkowski J. __Developmental restriction of retrotransposition activated in Arabidopsis by environmental stress__. (2017).

```r
# Please install the following packages if you haven't installed them yet
# source("http://bioconductor.org/biocLite.R")
# biocLite('biomartr')
# install.packages("dplyr")
# install.packages("readxl")
# install.packages("ggplot2")

source("FUNCS.R")

# download A. thaliana genome from ENSEMBLGENOMES
Ath_genome <- biomartr::getGenome(organism = "Arabidopsis thaliana", db = "ensemblgenomes")
# import A. thaliana genome as Biostrings object
Ath_genome_seq <- biomartr::read_genome(Ath_genome)
# rename chromosomes for consistency
Ath_genome_seq_chr <- Ath_genome_seq[1:5]
names(Ath_genome_seq_chr) <- paste0("Chr", 1:5)
# import ONSEN insertion excel file
ONSEN_insertions <- readxl::read_xlsx("List_new_insertions_all_lines.xlsx", sheet = 1)
# rename chromosomes for consistency
ONSEN_insertions <- dplyr::mutate(ONSEN_insertions, chromosome = paste0("Chr", ONSEN_insertions$chromosome))
# compute length of for all insertions
ONSEN_insertions <- dplyr::mutate(ONSEN_insertions, width = end - start + 1)
# store ONSEN insertion data as *.bed file for later use with BEDTOOLS
write.table(
        dplyr::select(ONSEN_insertions, chromosome:gene),
        file = "ONSEN_insertions.bed",
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
)

# partition A. thaliana genomes into sliding windows of 1Mbps and store as *.bed file
# for later use with BEDTOOLS
partition2bed(data = partition_genome(Ath_genome_seq_chr, interval = 1000000), "Ath_partition_1Mb.bed")

# sample random loci from the same chromosomes and the same lengths as 
# new ONSEN insertions and save results as *.bed file for later use with BEDTOOLS
random_ONSEN_inserts_bed <- generate_random_loci(Ath_genome_seq_chr, ONSEN_insertions)
write.table(
        random_ONSEN_inserts_bed,
        file = "Ath_random_insertions.bed",
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
)

# Generate euchromatin and heterochromatin intersections
system(
        "bedtools intersect -a Ath_partition_1Mb.bed -b Ath_euchromatin.bed -c > Ath_partition_1Mb_euchromatin.bed"
)

system(
        "bedtools intersect -a Ath_partition_1Mb.bed -b Ath_heterochromatin.bed -c > Ath_partition_1Mb_heterochromatin.bed"
)

# Euchromatic regions
Ath_partition_1Mb_euchromatin <-
        readr::read_tsv("Ath_partition_1Mb_euchromatin.bed", col_names = FALSE)
names(Ath_partition_1Mb_euchromatin) <-
        c("chr", "start", "end", "count")

Ath_partition_1Mb_euchromatin <-
        dplyr::filter(Ath_partition_1Mb_euchromatin, count == 1)
partition2bed(
        dplyr::select(Ath_partition_1Mb_euchromatin, 1:3),
        "Ath_partition_1Mb_euchromatin_filtered.bed"
)

SlidingWindow_count_ONSEN_insertions_euchromatin <-
        count_overlap("Ath_partition_1Mb_euchromatin_filtered.bed",
                      "ONSEN_insertions.bed")

SlidingWindow_count_ONSEN_insertions_euchromatin <-
        dplyr::mutate(
                SlidingWindow_count_ONSEN_insertions_euchromatin,
                ID = paste0(
                        SlidingWindow_count_ONSEN_insertions_euchromatin$chr,
                        "_",
                        SlidingWindow_count_ONSEN_insertions_euchromatin$start,
                        "_",
                        SlidingWindow_count_ONSEN_insertions_euchromatin$end
                )
        )

SlidingWindow_count_ONSEN_insertions_euchromatin <-
        dplyr::filter(SlidingWindow_count_ONSEN_insertions_euchromatin,
                      ID != "Chr4_2000001_3000000")

# visualize new ONSEN insertion distribution
p1 <-
        plotOverlapCounts(
                SlidingWindow_count_ONSEN_insertions_euchromatin,
                ymax = 15,
                main = "ONSEN insertion counts in euchromatic loci (332 out of 338 insertions)",
                with_facet = FALSE
        )

SlidingWindow_count_random_ONSEN_insertions_euchromatin <-
        count_overlap("Ath_partition_1Mb_euchromatin_filtered.bed",
                      "Ath_random_insertions.bed")

# visualize random insertion distribution
p2 <-
        plotOverlapCounts(
                SlidingWindow_count_random_ONSEN_insertions_euchromatin,
                ylab = "Random insertion counts",
                ymax = 15,
                main = "Randomly sampled 'artificial' insertions in euchromatic loci",
                with_facet = FALSE
        )


Random_vars_1000_euchromatin <- random_insertion_vars(
        geneome          = Ath_genome_seq_chr,
        insertion_data   = ONSEN_insertions,
        genome_partition = "Ath_partition_1Mb_euchromatin_filtered.bed",
        iterations       = 1000,
        chr_separately = FALSE
)

true_vars_euchromatin <-
        tibble::tibble(vars = var(
                SlidingWindow_count_random_ONSEN_insertions_euchromatin$count
        ))

# visualise
plotRandomVars(Random_vars_1000_euchromatin,
               true_vars_euchromatin,
               chr_separately = FALSE)

gamma_p_vals(Random_vars_1000_euchromatin$vars,
             true_vars_euchromatin$vars)
# [1] 0.2203609



# Heterochromatic regions
Ath_partition_1Mb_heterochromatin <-
        readr::read_tsv("Ath_partition_1Mb_heterochromatin.bed", col_names = FALSE)
names(Ath_partition_1Mb_heterochromatin) <-
        c("chr", "start", "end", "count")

Ath_partition_1Mb_heterochromatin <-
        dplyr::filter(Ath_partition_1Mb_heterochromatin, count == 1)
partition2bed(
        dplyr::select(Ath_partition_1Mb_heterochromatin, 1:3),
        "Ath_partition_1Mb_heterochromatin_filtered.bed"
)


SlidingWindow_count_ONSEN_insertions_heterochromatin <-
        count_overlap("Ath_partition_1Mb_heterochromatin_filtered.bed",
                      "ONSEN_insertions.bed")

SlidingWindow_count_ONSEN_insertions_heterochromatin <-
        dplyr::mutate(
                SlidingWindow_count_ONSEN_insertions_heterochromatin,
                ID = paste0(
                        SlidingWindow_count_ONSEN_insertions_heterochromatin$chr,
                        "_",
                        SlidingWindow_count_ONSEN_insertions_heterochromatin$start,
                        "_",
                        SlidingWindow_count_ONSEN_insertions_heterochromatin$end
                )
        )

SlidingWindow_count_ONSEN_insertions_heterochromatin <-
        dplyr::filter(
                SlidingWindow_count_ONSEN_insertions_heterochromatin,
                ID != "Chr3_15000001_16000000"
        )

# visualize new ONSEN insertion distribution
p3 <-
        plotOverlapCounts(
                SlidingWindow_count_ONSEN_insertions_heterochromatin,
                ymax = 15,
                main = "ONSEN insertion counts in heterochromatic loci (6 out of 338 insertions)",
                with_facet = FALSE
        )


SlidingWindow_count_random_ONSEN_insertions_heterochromatin <-
        count_overlap("Ath_partition_1Mb_heterochromatin_filtered.bed",
                      "Ath_random_insertions.bed")

# visualize random insertion distribution
p4 <-
        plotOverlapCounts(
                SlidingWindow_count_random_ONSEN_insertions_heterochromatin,
                ylab = "Random insertion counts",
                ymax = 15,
                main = "Randomly sampled 'artificial' insertions in heterochromatic loci",
                with_facet = FALSE
        )



Random_vars_1000_heterochromatin <- random_insertion_vars(
        geneome          = Ath_genome_seq_chr,
        insertion_data   = ONSEN_insertions,
        genome_partition = "Ath_partition_1Mb_heterochromatin_filtered.bed",
        iterations       = 1000,
        chr_separately = FALSE
)

# visualise

true_vars_heterochromatin <-
        tibble::tibble(vars = var(
                SlidingWindow_count_random_ONSEN_insertions_heterochromatin$count
        ))

# visualise
plotRandomVars(Random_vars_1000_heterochromatin,
               true_vars_heterochromatin,
               chr_separately = FALSE)

gamma_p_vals(
        Random_vars_1000_heterochromatin$vars,
        true_vars_heterochromatin$vars,
        lower.tail = TRUE
)


p5 <- gridExtra::grid.arrange(p1, p2, p3, p4)


```





