# ...
# Spp & Speciation module
# Hons

# Ruan van Mazijk
# 2017-08-08

# Setup ------------------------------------------------------------------------

rm(list = ls())
if (!require(pacman)) install.packages("pacman")
p_load(tidyverse, magrittr, here, stringr, glue, rjson)


# Specificities of the setup in this branch ------------------------------------

set_here(here("Results"))
do_plots <- FALSE  # used below to turn off all the plot devices for source()

# Tidy up the MSW3 data --------------------------------------------------------

msw <- as_tibble(read.csv(here("Results", "msw3-all.csv")))
names(msw)
# I hate factors
unfactor_to_chr <- function(x) as.character(x)
unfactor_to_num <- function(x) as.numeric(unfactor_to_chr(x))
msw[, c(18, 19)] %<>% map(unfactor_to_num) 
msw[, -c(1, 18, 19)] %<>% map(unfactor_to_chr)
# But sometimes ggplot forces me to embrace them:
msw$TaxonLevel %<>% as.factor()
msw$Order %<>%
    str_to_title() %>%
    as.factor()


# New functions for plots ------------------------------------------------------

plot_new_taxa <- function(x, rank, from = 1975, ylab = "taxa", ymax = 22,
                          by = "Order", which_orders = NULL, which_by = NULL,
                          drop_fill = FALSE, drop_facet = FALSE) {
    if (not(is.null(which_by))) {
        x <- x[x[, by][[1]] %in% which_by, ]
    }
    x %>%
        filter(TaxonLevel == rank,
               Date >= from,
               Order %in% if (is.null(which_orders)) x$Order
                          else which_orders) %>%
        ggplot() +
        geom_histogram(aes_string(x = "Date", fill = by),
                       bins = 2005 - from) +
        scale_fill_discrete(drop = drop_fill) +
        xlim(from, 2005) +
        ylim(0, ymax) +
        xlab("Year") +
        ylab(glue("No. new {ylab}")) +
        ggtitle(label = if (rank != "SPECIES") rank else NULL) +
        facet_wrap(by, drop = drop_facet) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90),
              panel.grid = element_blank())
}
save_new_taxa <- function(x, ranks, suffix = "", width = 8, height = 6) {
    pmap(list(x = x, y = ranks), function(x, y) {
        ggsave(width = width,
               height = height,
               plot = x,
               filename = here(glue("new_taxa_msw_{y}{suffix}.pdf")))
    })
}
get_the_coolest <- function(x, rank, by = "SPECIES", cutoff = 10) {
    x[(x$Date >= 1975) &
      (x$TaxonLevel == by), ] %>%
    count_(vars = rank) %>%
    filter(n >= cutoff) %>%
    `[`(rank) %>%
    `[[`(1)
}
search_all_columns <- function(x, query) {
    out <- list()
    for (i in 1:nrow(x)) {
        if (any(str_detect(x[i, ], query))) {
            out[[i]] <- x[i, ]
        }
    }
    match_maybe <- map(out, function(x) not(is.null(x)))
    match_maybe <- out[unlist(match_maybe)]
    n_match_maybe <- length(match_maybe)
    out_tidy <- data.frame()
    for (i in 1:n_match_maybe) {
        out_tidy %<>% rbind(match_maybe[[i]])
    }
    return(out_tidy)
}


# Plot the no. new taxa ~ time -------------------------------------------------

levels(msw$TaxonLevel)
levels(msw$Order)

# The ranks I want to see ~ time
all_ranks <- c("FAMILY", "SUBFAMILY", "TRIBE",
               "GENUS", "SUBGENUS",
               "SPECIES", "SUBSPECIES")

if (do_plots) {

    # Make & save all those plots
    all_ranks %>%
        map(plot_new_taxa, x = msw) %T>%
        save_new_taxa(all_ranks)
    
    # Plot & save the HIGHLIGHTED no. new taxa ~ time
    all_ranks %>%
        map(plot_new_taxa, x = msw, drop_facet = TRUE) %T>%
        save_new_taxa(all_ranks, suffix = "2")
    
    # Now let's REALLY highlight the new SPECIES each year
    plot_new_taxa(msw, "SPECIES",
        by = "Order",
        which_orders = get_the_coolest(msw, "Order"),
        ylab = "spp.",
        drop_fill = TRUE,
        drop_facet = TRUE) %>%
        list() %>%
        save_new_taxa("SPECIES", suffix = "3", width = 5, height = 4)
    # Total no. new spp since 1975:
    msw %>%
        filter(TaxonLevel == "SPECIES", Date >= 1975) %>%
        nrow()  # 564
    # And for the "cool" orders:
    msw %>%
        filter(TaxonLevel == "SPECIES",
               Date >= 1975,
               Order %in% get_the_coolest(msw, "Order")) %>%
        nrow()  # 533
    
    
    # And some more
    plot_new_taxa(msw, "SPECIES",
        by = "Family",
        which_by = get_the_coolest(msw, "Family", cutoff = 7),  # aot cutoff=10
        drop_fill = TRUE,
        drop_facet = TRUE)
    
    plot_new_taxa(x = msw, rank = "GENUS",
        by = "Family",
        which_by = get_the_coolest(msw, "Family", by = "GENUS", cutoff = 3),
        which_orders = get_the_coolest(msw, "Order", by = "GENUS"),
        drop_fill = TRUE,
        drop_facet = TRUE)
    
}

    
# Plot no. taxa ~ time GIVEN diversity in taxa ---------------------------------

richnesses <- msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order")) %>% 
    count(Order)
new_spp_counts <- msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order")) %>% 
    select(Date, Order) %>%
    group_by(Order) %>%
    count(Date)
rich_stdised_n <- vector(length = nrow(new_spp_counts))
for (i in 1:nrow(new_spp_counts)) {
    rich_stdised_n[i] <-
        new_spp_counts$n[i] %>%
        divide_by(richnesses$n[richnesses$Order == new_spp_counts$Order[i]])
}
new_spp_counts <- tibble(Order = new_spp_counts$Order,
                         Date = new_spp_counts$Date,
                         n = new_spp_counts$n,
                         rich_stdised_n)

if (do_plots) {
    ggplot(new_spp_counts) +
        geom_col(aes(x = Date, y = rich_stdised_n, fill = Order)) +
        xlim(1975, 2005) +
        xlab("Year") +
        ylab(glue("No. new spp. ÷ richness")) +
        facet_wrap(~ Order) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90),
              panel.grid = element_blank())
}

# Check rich_stdised_n for all orders ------------------------------------------
# (except those with ONLY 1 species discovery in the 1975:2005 period)

richnesses <- msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order", cutoff = 2)) %>%
    count(Order)
new_spp_counts <- msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order", cutoff = 2)) %>% 
    select(Date, Order) %>%
    group_by(Order) %>%
    count(Date)
rich_stdised_n <- vector(length = nrow(new_spp_counts))
for (i in 1:nrow(new_spp_counts)) {
    rich_stdised_n[i] <-
        new_spp_counts$n[i] %>%
        divide_by(richnesses$n[richnesses$Order == new_spp_counts$Order[i]])
}
new_spp_counts <- tibble(Order = new_spp_counts$Order,
                         Date = new_spp_counts$Date,
                         n = new_spp_counts$n,
                         rich_stdised_n)

if (do_plots) {
    ggplot(new_spp_counts) +
        geom_col(aes(x = Date, y = rich_stdised_n, fill = Order)) +
        xlim(1975, 2005) +
        xlab("Year") +
        ylab(glue("No. new spp. ÷ richness")) +
        facet_wrap(~ Order) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90),
              panel.grid = element_blank())
}

# What about log(count(Date) + 1)? ---------------------------------------------

new_spp_counts <- msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order", cutoff = 2)) %>% 
    select(Date, Order) %>%
    group_by(Order) %>%
    count(Date) %>%
    mutate(log_n = log(n + 1))

if (do_plots) {
    new_spp_counts %>%
        filter(Order %in% get_the_coolest(msw, "Order")) %>%
        ggplot() +
        geom_col(aes(x = Date, y = log_n, fill = Order)) +
        xlim(1975, 2005) +
        xlab("Year") +
        ylab(expression(paste(italic("log"["e"]), "[ No. new spp. + 1 ]"))) +
        facet_wrap(~ Order) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90),
              panel.grid = element_blank())
}


# Get the most common journals -------------------------------------------------

if (do_plots) {
    
    # For any taxon record in the whole MSW3 database
    msw$CitationName %>%
        as.factor() %>%
        summary() %T>%
        write.csv("msw3_CitationName_frequencies.csv")
    
    # For sub-sp., species OR genus records since 1975
    msw %>%
        filter(TaxonLevel %in% c("SUBSPECIES", "SPECIES", "GENUS"),
               Date >= 1975) %>%
        count(CitationName) %>%
        arrange(desc(n)) %T>%
        write.csv("msw3_CitationName_frequencies_subspORspORgen_since1975.csv")
    
    # For species OR genus records since 1975
    msw %>%
        filter(TaxonLevel %in% c("SPECIES", "GENUS"), Date >= 1975) %>%
        count(CitationName) %>%
        arrange(desc(n)) %T>%
        write.csv("msw3_CitationName_frequencies_spORgen_since1975.csv")
    
    # For species records since 1975
    msw %>%
        filter(TaxonLevel == "SPECIES", Date >= 1975) %>%
        count(CitationName) %>%
        arrange(desc(n)) %T>%
        write.csv("msw3_CitationName_frequencies_sp_since1975.csv")
    
    # For the coolest_orders' records since 1975
    map(get_the_coolest(msw, "Order"), function(x) {
        msw %>%
            filter(TaxonLevel == "SPECIES",
                   Date >= 1975,
                   Order %in% x) %>%
            count(CitationName) %>%
            arrange(desc(n)) %T>%
            write.csv(glue("msw3_CitationName_frequencies_{x}_sp_since1975.csv"))
    })

}


# Region information -----------------------------------------------------------

# Cape region
cape_msw <- filter(msw,
                   TaxonLevel == "SPECIES",
                   Date >= 1975,
                   str_detect(Distribution, "Cape"))
cape_msw %>% select(Distribution)
cape_msw %>% count(CitationName)
cape_msw %>%
    select(Order, Genus, Species, CitationName) %>%
    mutate(binom = paste(Genus, Species))

# Cape "all_column" search
cape_sl_msw <- msw %>%
    search_all_columns("[Ww]estern [Cc]ape") %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order")) %>%
    select(ID, Order, Genus, Species, CitationName, Distribution)
paste(cape_sl_msw$Genus, cape_sl_msw$Species)

# South Africa
SA_sl_msw <- msw %>%
    search_all_columns("[Ss]outh [Aa]frica") %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975) %>%
           #Order %in% get_the_coolest(msw, "Order")) %>%
    select(ID, Order, Genus, Species, CitationName, Distribution)
paste(SA_sl_msw$Genus, SA_sl_msw$Species)
SA_sl_msw

# Africa
afr_msw <- filter(msw,
                  TaxonLevel == "SPECIES",
                  Date >= 1975,
                  str_detect(Distribution, "Africa"))
afr_msw %>% select(Distribution)
afr_msw %>% count(CitationName)
afr_msw %>%
    select(Order, Genus, Species, CitationName) %>%
    mutate(binom = paste(Genus, Species))

# Appearing in Ann. Transvaal Mus. after 1975
transv_msw <- filter(msw,
                     TaxonLevel == "SPECIES",
                     Date >= 1975,
                     str_detect(CitationName, "Transvaal"))
transv_msw %>%
    select(Order, Genus, Species, CitationName) %>%
    mutate(binom = paste(Genus, Species))

# Search "Afr" in all_columns
afr_sl_msw <- search_all_columns(msw, "[Aa]fr")
if (do_plots) {
    plot_new_taxa(afr_sl_msw, "SPECIES",
        from = 1975,
        by = "Order",
        which_orders = get_the_coolest(msw, "Order"),
        ylab = "spp.",
        ymax = 3,
        drop_fill = TRUE,
        drop_facet = TRUE) %>%
        list() %>%
        save_new_taxa("SPECIES", suffix = "_AFRICA", width = 6, height = 2)
    plot_new_taxa(afr_sl_msw, "SPECIES",
        from = 1915,
        by = "Order",
        which_orders = get_the_coolest(msw, "Order"),
        ylab = "spp.",
        ymax = 5,
        drop_fill = TRUE,
        drop_facet = TRUE) %>%
        list() %>%
        save_new_taxa("SPECIES", suffix = "_AFRICA2", width = 4, height = 4)
}
# No. new spp. total?
afr_sl_msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1975,
           Order %in% get_the_coolest(msw, "Order")) %>% 
    nrow()  # 523
afr_sl_msw %>%
    filter(TaxonLevel == "SPECIES",
           Date >= 1915,
           Order %in% get_the_coolest(msw, "Order")) %>% 
    nrow()  # 93


# Literature searching ---------------------------------------------------------

lit_search <- fromJSON(file = "lit_search.json")
journals <- map(lit_search, function(x) {
    x[1][[1]] %>%
        str_replace_all(" ", "_") %>%
        str_replace_all("\\.", "")
})
names(lit_search) <- journals
listviewer::jsonedit(lit_search)

map(lit_search$J_Mammal$results, function(x) x$taxa)


# JUNK -------------------------------------------------------------------------

## Frame plot 3 in the context of all orders
#
#plot4 <- msw %>%
#    filter(
#        TaxonLevel == "SPECIES",
#        Date >= 1975
#    ) %>%
#    ggplot() +
#    geom_histogram(
#        aes_string(x = "Date", fill = "Order"),
#        bins = 2005 - 1975
#    ) +
#    #scale_fill_discrete(drop = FALSE) +  # TODO?
#    xlim(1975, 2005) +
#    ylim(0, 22) +
#    xlab("Year") +
#    ylab("No. new spp.") +
#    facet_wrap("Order", drop = FALSE) +
#    theme_bw() +
#    theme(
#        legend.position = "none",
#        axis.text.x = element_text(angle = 90),
#        panel.grid = element_blank()
#    )
#
#ggsave(
#    width = 8, height = 6,
#    plot = plot4,
#    filename = here(paste0("new_taxa_msw_SPECIES4.pdf"))
#)
#
#
## Visualise ordinal-contribution to new taxa in each year
## FIXME
#plot_new_taxa3 <- function(x,
#                           which_rank = c(1, 2),
#                           rank = c("SPECIES", "SUBSPECIES"),
#                           from = 1975,
#                           by = "Order") {
#    x %>% filter(TaxonLevel %in% rank, Date >= from)
#    x$TaxonLevel %<>% as.factor()
#    x %>%
#        filter(TaxonLevel == rank[which_rank]) %>%
#        ggplot() +
#        geom_dotplot(
#            aes_string(x = "Date", fill = by, colour = by),
#            binwidth = 1,
#            stackgroups = TRUE,
#            method = "histodot",
#            binpositions = "all"
#        ) +
#        scale_fill_discrete(drop = FALSE) +
#        scale_colour_discrete(drop = FALSE) +
#        xlim(from, 2005) +
#        ylim(0, 40) +
#        xlab("Year") +
#        ylab("No. new taxa") +
#        ggtitle(rank[which_rank]) +
#        theme_bw() +
#        theme(axis.text.x = element_text(angle = 90))
#}
#
#plot_new_taxa3(msw, 1)
#plot_new_taxa3(msw, 2)


# </> --------------------------------------------------------------------------

