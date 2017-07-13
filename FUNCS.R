#' @title Partition Genome into sliding windows
#' @description A given genome assemby as Biostrings object is
#' parititioned into intervals of length \code{n}.
#' @param genome Biostrings object storing a genome assembly as DNA strings.
#' @param interval Interval length of partitioned sliding windows.
#' @details This function partitions a genome assembly file into sliding windows
#' of a specified \code{interval} length and returns a \code{tibble} in
#' \code{*.bed} file format which can then be stored in a \code{*.bed}
#' file using the function \code{partition2bed}.
#' @author Hajk-Georg Drost
partition_genome <- function(genome, interval = 1000000) {
        chr_n <- length(genome)
        chr_lengths <- genome@ranges@width
        chr_names <- genome@ranges@NAMES
        chr_df <- vector("list", chr_n)
        
        for (i in seq_len(chr_n)) {
                split_intervals <- seq(1, chr_lengths[i], by = interval)
                chr_df[i] <-
                        list(tibble::tibble(
                                chr   = rep(chr_names[i], length(split_intervals)),
                                start = as.integer(split_intervals),
                                end   = as.integer(c(
                                        split_intervals[2:length(split_intervals)] - 1, chr_lengths[i]
                                ))
                        ))
                
        }
        
        res <- dplyr::bind_rows(chr_df)
        return(res)
}

#' @title Store \code{partition_genome} output as \code{*.bed} file.
#' @description Writes \code{partition_genome} intervals into \code{*.bed} file.
#' @param data a tibble in *.bed file format that is returned by \code{partition_genome}.
#' @param file file name (file path) of the *.bed file that shall be generated.
#' @author Hajk-Georg Drost
partition2bed <- function(data, file) {
        write.table(
                data,
                file = file,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE
        )
}

#' @title Generate random loci that are analogous to new ONSEN insertions
#' @description 
#' @param genome Biostrings object storing a genome assembly as DNA strings.
#' @param insertion_data ONSEN insertion data as \code{tibble}.
#' @details 
#' @author Hajk-Georg Drost
generate_random_loci <- function(genome, insertion_data) {
        chr_n <- length(genome)
        chr_lengths <- genome@ranges@width
        chr_names <- genome@ranges@NAMES
        chr_df <- vector("list", chr_n)
        
        for (i in seq_len(chr_n)) {
                chr_positions <- vector("numeric", chr_lengths[i])
                chr_positions <- seq(1, chr_lengths[i], by = 1)
                insertion_data_chr <-
                        dplyr::filter(insertion_data, chromosome == chr_names[i])
                sample_insertions <-
                        sample(chr_positions, nrow(insertion_data_chr))
                
                
                chr_df[i] <- list(tibble::tibble(
                        chr   = rep(chr_names[i], nrow(insertion_data_chr)),
                        start = as.integer(sample_insertions),
                        end   = as.integer(sample_insertions + insertion_data_chr$width)
                ))
                
        }
        res <- dplyr::bind_rows(chr_df)
        return(res)
}


#' @title Count for each sliding window how many ONSEN insertions were found
#' @description Using bedtools to count for each sliding window how many ONSEN insertions were found.
#' @param x partitioned genome file in *.bed format.
#' @param y ONSEN insertions file in *.bed format.
#' @author Hajk-Georg Drost
count_overlap <- function(x, y) {
        output_file <- file.path(tempdir(), "x_y_count.bed")
        
        # using bedtools to count for each feature in the “A” file,
        # the number of overlapping features in the “B” file.
        system(paste0("bedtools intersect -a ", x, " -b ", y, " -c > ", output_file))
        # import bedtools intersect output
        res <- readr::read_tsv(
                output_file,
                col_names = FALSE,
                col_types = readr::cols(
                        X1 = readr::col_character(),
                        X2 = readr::col_integer(),
                        X3 = readr::col_integer(),
                        X4 = readr::col_integer()
                )
        )
        
        names(res) <- c("chr", "start", "end", "count")
        return(res)
}

#' @title Visualize ONSEN insertion counts. 
#' @description Visualize ONSEN insertion count for each sliding window and chromosome. 
#' @param data a \code{tibble} returned by \code{count_overlap}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main title text.
#' @param ymax maximum value on y-axis (= ylim max).
#' @param with_facet shall counts be separately shown for each chromosome?
#' @author Hajk-Georg Drost
plotOverlapCounts <-
        function(data,
                 xlab = "Sliding Window",
                 ylab = "ONSEN Insertion Counts",
                 main = "",
                 ymax = 5, 
                 with_facet = TRUE) {
                
                data <-
                        dplyr::mutate(data, ID = paste0(data$chr, "_", data$start, "_", data$end))
                p <-
                        ggplot2::ggplot(data, ggplot2::aes(x = ID, y = count))
                
                if (with_facet)
                p <- p + ggplot2::facet_grid(chr ~ .)
                
                if (!with_facet)
                p <- p +
                        ggplot2::geom_point(size = 6, ggplot2::aes(colour = chr))
                
                if (with_facet)
                        p <- p +
                        ggplot2::geom_point(size = 6)
                
                p <- p +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(x = xlab,
                                      y = ylab,
                                      title = main) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 22, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 22, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 10, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.y     = ggplot2::element_text(
                                        size           = 22,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::expand_limits(y = c(-1, ymax)) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(
                                angle = 90,
                                vjust = 1,
                                hjust = 1
                        )) +
                        ggplot2::theme(panel.spacing = ggplot2::unit(2, "lines"))
                
                return(p)
                
        }


#' @title Boxplot visualising the \code{count_overlap} distributions between chromosomes.
#' @description The boxplot illustrates the differences between 
#' \code{count_overlap} distributions between chromosomes and internally performes 
#' a Kruskal-Wallis Rank Sum Test to assess the statistical difference between chromosome
#' distributions.
#' @param data a \code{tibble} returned by \code{count_overlap}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param ymax maximum value on y-axis (= ylim max).
#' @author Hajk-Georg Drost
plotChrDiff <-
        function(data,
                 xlab = "Arabidopsis Chromosome",
                 ylab = "Relative ONSEN Insertion Counts Per Sliding Window",
                 ymax = 15) {
                
                data2 <- dplyr::summarise(dplyr::group_by(data, chr), n = n())
                rel <- vector("list", nrow(data2))
                
                for (i in seq_len(nrow(data2))) {
                        
                        indiv_chr <- dplyr::filter(data, chr == data2$chr[i])
                        rel[i] <- list(indiv_chr$count / data2$n[i])
                }
                
                data <- dplyr::mutate(data, rel = unlist(rel))
                
                p <-
                        ggplot2::ggplot(data, ggplot2::aes(
                                x = chr,
                                y = rel,
                                colour = chr
                        )) +
                        ggplot2::geom_boxplot(size = 3) +
                        ggplot2::geom_point(size = 5) +
                        ggplot2::geom_jitter(position = ggplot2::position_jitter(0.1),
                                             size = 5) +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(
                                x = xlab,
                                y = ylab,
                                title = paste0("p-value = ", round(
                                        kruskal.test(data$rel, g = factor(data$chr))$p.value, 3
                                ))
                        ) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 22, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 22, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 22, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 22,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::expand_limits(y = c(0, ymax))
                
                message("Performing Kruskal-Wallis Rank Sum Test...")
                return(p)
        }


#' @title Generate random loci that are analogous to ONSEN insertions
#' @description Statistical assessment of \code{count_overlap} variance of ONSEN
#' insertions versus count variances of \code{n} times randomly sampled loci. 
#' @param genome Biostrings object storing a genome assembly as DNA strings.
#' @param insertion_data ONSEN insertion data as \code{tibble}.
#' @param genome_partition a \code{tibble} generated with \code{partition_genome}.
#' @param iterations Number of independent sampling runs.
#' @param chr_separately shall statsitics be performed separately for each chromosome 
#' (\code{chr_separately = TRUE}) or as a concatenation of all chromosomes ((\code{chr_separately = FALSE})). 
#' @details For 338 ONSEN insertions analogously 338 randomly sampled loci coming from the
#' same chromosomes and having the same lengths as 338 ONSEN insertions are generated
#' to test if randomly sampled loci can generate a similar sliding window count pattern
#' (see \code{plotOverlapCounts}) as the 338 ONSEN insertions.
#' @author Hajk-Georg Drost
random_insertion_vars <-
        function(geneome,
                 insertion_data,
                 genome_partition,
                 iterations = 10,
                 chr_separately = TRUE) {
                
                vars_random <- vector("list", iterations)
                
                for (i in seq_len(iterations)) {
                        message("Iteration ", i, " ...")
                        random_ONSEN_inserts_bed <-
                                generate_random_loci(geneome, insertion_data)
                        write.table(
                                random_ONSEN_inserts_bed,
                                file = file.path(tempdir(), "Ath_random_insertions.bed"),
                                row.names = FALSE,
                                col.names = FALSE,
                                sep = "\t",
                                quote = FALSE
                        )
                        sliding_window_random_insertion_count <-
                                count_overlap(
                                        genome_partition,
                                        file.path(tempdir(), "Ath_random_insertions.bed")
                                )
                        
                        # compute count variance of randomly sampled insertions per chromosome
                        if (chr_separately) {
                                vars_random[i] <- list(dplyr::summarise(
                                        dplyr::group_by(sliding_window_random_insertion_count, chr),
                                        vars = var(count)
                                ))
                        }
                        
                        if (!chr_separately)
                                vars_random[i] <- var(sliding_window_random_insertion_count$count)
                }
                
                if (chr_separately) {
                        res <- dplyr::bind_rows(vars_random)
                        return(res)
                } else {
                        return(tibble::tibble(vars = unlist(vars_random)))
                }
        }

#' @title Visualize random variances versus ONSEN count variance.
#' @description Visualize Statistical assessment of \code{count_overlap} variance of ONSEN
#' insertions versus count variances of \code{n} times randomly sampled loci.
#' @param data a \code{tibble} returned by \code{random_insertion_vars}.
#' @param true_vars a \code{tibble} storing the true count variances of ONSEN insertions.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param chr_separately shall statsitics be performed separately for each chromosome 
#' (\code{chr_separately = TRUE}) or as a concatenation of all chromosomes ((\code{chr_separately = FALSE})). 
#' @author Hajk-Georg Drost

plotRandomVars <-
        function(data,
                 true_vars,
                 xlab = "Variances of randomly sampled insertions",
                 ylab = "Frequency",
                 chr_separately = TRUE) {
                
                p <- ggplot2::ggplot(data, ggplot2::aes(x = vars))
                if (chr_separately)
                        p <- p + ggplot2::facet_grid(chr ~ .)
                
                p <- p +
                        ggplot2::geom_histogram(bins = 40) +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(x = xlab, y = ylab) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 22, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 22, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 22, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 22, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.y     = ggplot2::element_text(
                                        size           = 22,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::geom_vline(
                                data = true_vars,
                                ggplot2::aes(xintercept = vars, colour = "red"),
                                size = 2.5,
                                show.legend = FALSE
                        ) +
                        ggplot2::theme(panel.spacing = ggplot2::unit(0.5, "lines"))
                return(p)
        }

#' @title Fitting gamma distribution
#' @description Fitting gamma distribution to input data and
#' using moment matching estimation to estimate 
#' the shape and scale parameters of the gamma distribution. 
#' @param x a numeric vector for which shape and scale parameters of 
#' the gamma distribution shall be estimated.
#' @param real_val a numeric value representing the true variance values of ONSEN count variances
#' @author Hajk-Georg Drost

gamma_p_vals <- function(x, real_val, lower.tail = FALSE) {
        gamma_MME <- fitdistrplus::fitdist(x, "gamma", method = "mme")
        ### estimate shape
        shape <- gamma_MME$estimate[1]
        ### estimate the rate
        rate <- gamma_MME$estimate[2]
        
        return(stats::pgamma(
                q = real_val,
                shape = shape,
                rate = rate,
                lower.tail = lower.tail
        ))
        
}

#' @title Fitting gamma distributions for all chromosomes individually
#' @description 
#' @param x a \code{tibble} storing the chromosome information in the first column
#' and the randomly sampled count variance for each chromosome in the second column.
#' @param true_vars a \code{tibble} storing the chromosome information in the first column
#' and the ONSEN count variance for each chromosome in the second column.
#' @author Hajk-Georg Drost

gamma_p_vals_all_chr <- function(x, true_vars) {
        chr_names <- names(table(x$chr))
        
        if (length(chr_names) != nrow(true_vars))
                stop(
                        "Please provide the same number of chromosomes for true variance values and random variance values..."
                )
        
        res <- vector("list", nrow(true_vars))
        for (i in seq_len(length(chr_names))) {
                Chr_random <- dplyr::filter(x, chr == chr_names[i])$vars
                res[i] <-
                        list(gamma_p_vals(Chr_random, true_vars$vars[i]))
        }
        
        res <- unlist(res)
        names(res) <- paste0("p-value: ", chr_names)
        return(res)
}
