sites_order <- c("tamasag", # rist
                 "legacy",
                 "lincoln",
                 "timberline",
                 "prospect",
                 "boxelder", # elc
                 "archery",
                 "river bluffs")

if(network == "virridy"){ # this will be the new default
  # establish non-tributary sites
  tributary_sites <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                     "tamasag","legacy", "lincoln","timberline",
                                     "timberline virridy","prospect",
                                     "prospect virridy","boxelder","archery",
                                     "archery virridy","riverbluffs"))
  # establish order for all the non-tributary sites
  sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                    "legacy","lincoln","timberline","prospect","boxelder",
                    "archery","riverbluffs")
  # establish order for all the tributary sites
  trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                        NA, "penn", "sfm", "lbea")
}

# determining the sites relevant to the site of interest.
if (site_arg %in% tributary_sites$site) {

  plot_filter <- tributary_sites

  site_parent <- strsplit(site_arg, " ")[[1]][1]
  site_index <- which(sites_order == site_parent)
  site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

  plot_filter <- plot_filter %>%
    filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
           site != site_arg) %>%
    pull(site)

} else {

  plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                 "springcreek", "prospect", "prospect virridy",
                                 "penn", "sfm", "lbea"))

  site_index <- which(trib_sites_order == site_arg)
  site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

  plot_filter <- plot_filter %>%
    filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
           site != site_arg) %>%
    pull(site)

}
