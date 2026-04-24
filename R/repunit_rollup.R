#' Repunit Rollup for Ms.GSI
#'
#' @param mdl_out Optional. Ms.GSI output object name for combining group proportions and harvest of a single mixture.
#' @param new_pop_info Population information for the new grouping. A tibble with columns `collection`, `repunit` and `new_repunit`. `repunit` is the names of the original reporting groups. Can include a column for `collection` if reorganizing using collections.
#'
#' @return A tibble of proportions and harvest numbers by reporting group for combined mixtures/strata.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' new_groups = msgsi_dat$comb_groups %>%
#'  mutate(new_repunit = case_when(
#'    repunit=="Lower Yukon"~ "Yukon",
#'    repunit=="Middle Yukon"~ "Yukon",
#'    repunit=="Upper Yukon"~ "Yukon",
#'    .default = repunit
#'  ))
#'
#' repunit_rollup_msgsi(mdl_out = msgsi_out, new_pop_info = new_groups)
#' }
#'
repunit_rollup_msgsi <- function(mdl_out = NULL, new_pop_info = NULL,CIs = c(0.05,0.95)) {

  if (is.null(mdl_out)) {
    stop("`mdl_out not specified, please provide a Ms.GSI output.")
  }

  if (is.null(new_pop_info)) {
    stop("`new_pop_info not specified, please provide tibble of collections to roll up.")
  }

  #Check to see if all values of new_pop_info exist

    if (!"collection" %in% names(new_pop_info)) {
      stop("Please provide collection in new_pop_info.")
    } else if (!"repunit" %in% names(new_pop_info)) {
      stop("Please provide collection in new_pop_info.")
      } else if (!"new_repunit" %in% names(new_pop_info)) {
        stop("Please provide new_repunit in new_pop_info.")
    }

  # group information ----
  #probably a good idea to check the users new groupings and make sure
  #what was analyzed is there, or throw an error

  grp_info <- mdl_out$comb_groups

  if (any(!(new_pop_info$collection %in% grp_info$collection))) {
    stop("Collections in new_pop_info do not match collections in mdl_out")
  }

  # calculations ----

  nburn <- as.numeric(mdl_out$specs["nburn"])
  harv <- mdl_out$sstc_trace_t1 %>%
    dplyr::slice_min(itr) %>%
    dplyr::slice_min(ch) %>%
    dplyr::select(-c(itr,ch)) %>%
    rowSums(.)


  #Summarize Stock Proportions
  MCMC_chains <- mdl_out$trace_comb %>%
    dplyr::filter(itr > nburn) %>%
    dplyr::group_split(ch) %>%
    lapply(function(df) {
      df %>%
        reshape2::melt(id.vars = c('itr', 'ch')) %>%
        dplyr::rename(collection = variable) %>%
        dplyr::left_join(new_pop_info, by = 'collection') %>%
        dplyr::group_by(new_repunit, itr) %>%
        dplyr::summarize(sum_rg = sum(value), .groups = "drop") %>%
        reshape2::dcast(formula = itr ~ new_repunit, value.var = 'sum_rg') %>%
        dplyr::select(-itr) %>%
        coda::mcmc()
    }) %>%
    coda::as.mcmc.list()

  MCMCsum <- coda:::summary.mcmc.list(MCMC_chains,quantiles=c(CIs[1],0.5,CIs[2]))

  MCMCconvg <- coda::gelman.diag(MCMC_chains,autoburnin = F, multivariate = F)

  MCMCsampsz <- coda::effectiveSize(MCMC_chains)

   P0 <- MCMC_chains %>%
     do.call(rbind,.) %>%
     reshape2::melt()%>%
     dplyr::group_by(Var2)%>%
     dplyr::summarise(P0 = mean(value < 5e-7)) %>%
     dplyr::rename('group' = 'Var2' )

   Res <- tibble::tibble(group = rownames(MCMCsum$statistics),
                         Mean = MCMCsum$statistics[,1],
                         SD = MCMCsum$statistics[,2],
                         tibble::as_tibble(MCMCsum$quantiles),
                         Eff_n = MCMCsampsz,
                         P0 = P0$P0,
                         GR_PtEst=MCMCconvg$psrf[,1],
                         GR_UpCI=MCMCconvg$psrf[,2])

  #summarize stock specific harvest
 # but only if the harvest in the sstc is positive - suggesting the user supplied it at the time of analysis
  if(harv > 0 ){

  coll_names_t1 <- names(mdl_out$sstc_trace_t1)[!names(mdl_out$sstc_trace_t1) %in% names(mdl_out$sstc_trace_t2)]

  sstc_ru <-  mdl_out$sstc_trace_t1[, c(coll_names_t1, "itr", "ch")] %>%
    dplyr::left_join(mdl_out$sstc_trace_t2, by = dplyr::join_by(itr, ch)) %>%
    tidyr::pivot_longer(-c(itr, ch), names_to = "collection") %>%
    dplyr::filter(collection %in% grp_info$collection,itr > nburn) %>%  #remove collections not in combined group
    dplyr::left_join(new_pop_info, by = "collection") %>%
    dplyr::summarise(sstc = sum(value), .by = c(itr, ch, new_repunit)) %>%
    dplyr::summarise(mean_sstc = mean(sstc),
                   sd_sstc = stats::sd(sstc),
                   median_sstc = stats::median(sstc),
                   cilo_sstc = stats::quantile(sstc, CIs[1]),
                   cihi_sstc = stats::quantile(sstc, CIs[2]),
                   `P=0` = mean(sstc < 0.5),
                   .by = c(new_repunit))%>%
    dplyr::rename( 'group' = 'new_repunit')


Res <- Res %>%
  dplyr::select(-P0) %>%
  dplyr::left_join(.,sstc_ru,by='group')
}

 Res


}
utils::globalVariables(c("nCritters", "variable", "harv_p", "Var2"))
