# New functions for GoF analyses, to be sourced as I work on the code in
# Goodness_of_fit.Rmd, and then copied into sizeSpectra once done.

##' Combine bins together (sequentially) as necessary to ensure that each has a
##'  count of 5 or more
##'
##' Required for the G-test. No combining is needed if all bins have counts of 5
##' or more. If, say, bin 7 has a count of 4 then it is combined with bin 8 (and
##' higher bins if necessary to get a count of 5 or more). If, after this
##' sequential calculation, the last bin has a
##' count of <5 then it is combined with the penultimate one (which by
##' definition will have a count of 5 or more).
##'
##' @param binBreaks vector of breakpoints of bins
##' @param binCounts vector of counts in bins (of length `length(binBreaks) - 1`)
##' @param minCounts minimum number of Counts required in each bin (default 5). TODO:
##'   generalise, 5 is hardwired
##' @return list containing objects:
##' \describe{
##'   \item{combined_breaks}{vector of combined bin breaks}
##'   \item{combined_counts}{vector of combined bin counts}
##' @export
##' @author Andrew Edwards
##' @examples
##' @dontrun{
##' test_that("combine_bins works correctly", {
##' Simple, no combining needed:
##' expect_equal(combine_bins(binBreaks = 1:7,
##'                           binCounts = c(5, 9, 10, 12, 8, 6)),
##'              list(combined_breaks = 1:7,
##'                   combined_counts = c(5, 9, 10, 12, 8, 6)))
##' # Some combining needed, including final bin being <5
##' expect_equal(combine_bins(binBreaks = 1:7,
##'                           binCounts = c(5, 4, 0, 2, 3, 4)),
##'              list(combined_breaks = c(1, 2, 5, 7),
##'                   combined_counts = c(5, 6, 7)))
##' # Some combining needed, including final bin being <10, minCounts = 10
##' expect_equal(combine_bins(binBreaks = 1:7,
##'                           binCounts = c(5, 4, 0, 2, 3, 4),
##'                           minCounts = 10),
##'              list(combined_breaks = c(1, 7),
##'                   combined_counts = 18))
##' # })
##' @}
combine_bins <- function(binBreaks,
                         binCounts,
                         minCounts = 5){
  stopifnot(length(binBreaks) == length(binCounts) + 1)
  stopifnot(sum(binCounts) >= minCounts)

  # Check each in turn, combine with subsequent one if binCounts[i] is < minCounts:
  combined_breaks <- vector()
  combined_counts <- vector()

  combined_counts_i <- 1      # counter for combined_counts vector
  combined_breaks[1] <- binBreaks[1]

  i <- 1                      # counter for binCounts
  while(i <= length(binCounts)){
    if(binCounts[i] >= minCounts){
      combined_counts[combined_counts_i] <- binCounts[i]
      combined_breaks[combined_counts_i + 1] <- binBreaks[i + 1]
      combined_counts_i <- combined_counts_i + 1
      i <- i + 1
    } else {
      # Combine the next bins to get >=5 total count
      cumul <- cumsum(binCounts[i:length(binCounts)])
      if(max(cumul) >= minCounts){
        first_to_five <- min(which(cumul >= minCounts))  # first index of remaining counts
                                                 #  to have cumulative count of
                                                 #  >= minCounts (default 5,
                                                 #  hence first_to_five),
                                                 #  will be index >=2

      } else {
        first_to_five <- length(cumul)           # just take them all and fix afterwards
      }

      combined_counts[combined_counts_i] <- cumul[first_to_five]
      combined_breaks[combined_counts_i + 1] <- binBreaks[i + first_to_five]  # TODO CHECK THAT, may
                                        # be off by 1
      combined_counts_i <- combined_counts_i + 1
      i <- i + first_to_five
    }
  }

  # But count in final bin may be <5, if so then combine with penultimate bin
  # (which should be >=5 by definition)
  M = length(combined_counts)
  if(combined_counts[M] < minCounts){
     combined_counts[M-1] <- combined_counts[M-1] + combined_counts[M]
     combined_counts <- combined_counts[-M]
     combined_breaks <- combined_breaks[-M]
  }

  stopifnot(sum(combined_counts) == sum(binCounts))

  return(list(combined_breaks = combined_breaks,
              combined_counts = combined_counts))
}

##' Goodness-of-fit test (G-test) on binned data for a given b, and (optionally)
##'  xmin and xmax
##'
##' Functionalising the G-test (also called likelihood ratio test) with
##' Willams's correction, as desribed in Sokal and Rohlf (1995;Biometry: the
##' principles and practice of statistics in biological research, Third
##' edition), and used in
##' Edwards et al. (2007)\url{https://www.nature.com/articles/nature06199} and
##' expanded upon in
##' Edwards (2011)\url{https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/10-1182.1}.
##' For raw (unbinned) data see `GoF_PLB_unbinned` which first calculates bins
##'  and then calls this function.
##' @param bin_breaks breakpoints of bins
##' @param bin_counts counts in each bin, where bin_counts[i] is the count in
##'   the bin spanning bin_breaks[i] to bin_breaks[i+1]. First and last
##'   bin_counts must be empty.
##' @param b MLE estimate of size-spectrum exponent b
##' @param K number of parameters that are estimated from the data; default of 1
##'  is recommended based on simulations (to be written up)
##' @param xmin minimum possible value of the data (if NULL then the minimum of
##'   first bin)
##' @param xmax maximum possible value of the data (if NULL then the maximum of
##'   final bin)
##' @param min_in_each_bin minimum count to have in each combined (merged) bin for
##'   goodness-of-fit test; bins are combined to ensure this condition is
##'   met. Default of 5 is from Sokal and Rohlf (1995).
##' @return list containing objects:
##' \describe{
##'   \item{combined_breaks}{breakpoints of combined bins such that each count is >=5}
##'   \itemcombined_counts}{counts in each bin}
##'   \item{n}{total number of counts}
##'   \item{expected_counts expected counts in each bin}
##'   \item{dof}{degrees of freedom used in test}
##'   \item{GWilliams}{G-test statistic, with Williams's correction}
##'   \item{crit95}{critical value; if the  = crit95,
##'   \item{Pvalue}{resulting P value}
##'   \item{consistent}{TRUE if the data are consistent with the PLB
##'   distribution (not significantly different from it at the 0.05 level;
##'   i.e. GWilliams < crit95); FALSE if data are not consistent with the PLB distribution}
##' }
##'
##' @export
##' @author Andrew Edwards
##' @examples
##' @dontrun{
##' @
##' @}
GoF_PLB <- function(bin_breaks,
                    bin_counts,
                    b,
                    K = 1,
                    xmin = NULL,
                    xmax = NULL,
                    min_in_each_bin = 5){
  stopifnot(bin_counts[1] > 0)
  stopifnot(bin_counts[length(bin_counts)] > 0)

  if(is.null(xmin)) xmin = min(bin_breaks)
  if(is.null(xmax)) xmax = max(bin_breaks)

  combined <- combine_bins(bin_breaks,
                           bin_counts,
                           minCounts = min_in_each_bin)


  # combined$combined_breaks
  # combined$combined_counts
  n <- sum(combined$combined_counts)
  expected_probs <- pPLB(x = combined$combined_breaks,
                         b = b,
                         xmin = min(combined$combined_breaks),
                         xmax = max(combined$combined_breaks)) %>%
    diff()
  expected_counts <- expected_probs * n

  dof <- length(combined$combined_counts) - K - 1   # degrees of freedom

  G <- 2 * sum( combined$combined_counts * log( combined$combined_counts / expected_counts) )    # p692 of S+R,
  qWilliams <-  1 + ( length(combined$combined_counts)^2 - 1 ) / (6 * n * dof)          # p698 of S+R
  GWilliams <- G / qWilliams

  crit95 <- qchisq(0.95, dof)    # pchisq(crit95, dof)=0.95
  Pvalue <- 1 - pchisq(GWilliams, dof)

  consistent <- (GWilliams < crit95)    # data are consistent (not sig
                                        # different) from the PLB

  return(list(combined_breaks = combined$combined_breaks,
              combined_counts = combined$combined_counts,
              n = n,
              expected_counts = expected_counts,
              dof = dof,
              GWilliams = GWilliams,
              crit95 = crit95,
              Pvalue = Pvalue,
              consistent = consistent))
}


##' Goodness-of-fit test (G-test) on raw unbinned data for a given b, and (optionally)
##'  xmin and xmax.
##'
##' Creates the bins for a given unbinned set of data and fitted `b`, and then
##' calls `GoF_PLB()`.
##'
##' @param x vector of raw unbinned values (e.g. body masses)
##' @param b.use value of `b`, already fitted
##' @param K number of parameters that are estimated from the data.
##' @param xmin minimum possible value of the data - if NULL then `min(x)`
##' @param xmax maximum possible value of the data - if NULL then `max(x) * (1 +
##'   epsilon)` to avoid rounding errors that say max of max bin is larger than `max(x).
##' @param increment.for.GoF.bins
##' @param epsilon tolerance for making the maximum bin `1+epsilon` times `max(x)`
##'   to avoid rounding errors that give error in `binData()` (as happened with second
##'   simulated data set using default settings and "unbinned" in
##'   `MLEbin.GoF.simulate()`).
##' @param ... arguments passed to `GoF_PLB()`, namely `min_in_each_bin`
##' @return list containing objects, see `GoF_PLB()`
##' @export
##' @author Andrew Edwards
##' @examples
##' @dontrun{
##'
##' }
GoF_PLB_unbinned <- function(x,
                             b.use,
                             K = 1,
                             xmin = NULL,
                             xmax = NULL,
                             increment.for.GoF.bins = 0.1,
                             epsilon = 1e-08,
                             ...
                             ){
  if(is.null(xmin)) xmin = min(x)
  if(is.null(xmax)) xmax = max(x)

  # Create the bins based on the fitted distribution
  binBreaks = qPLB(seq(0,
                       1,
                       by = increment.for.GoF.bins),
                   b = b.use,
                   xmin = min(x),
                   xmax = max(x))

  binBreaks[length(binBreaks)] = binBreaks[length(binBreaks)] * (1 + epsilon)

  bins.list = binData(x,
                      binBreaks = binBreaks)
  binCounts = bins.list$binVals[,"binCount"]$binCount

  GoF_res = GoF_PLB(bin_breaks = binBreaks,
                    bin_counts = binCounts,
                    b = b.use,
                    K = K,
                    ...)

  return(GoF_res)
}


##' Simulate, bin and fit data using four different binning methods and MLEbin
##'  and test goodness of fit with number of estimated paramters (K) set to 1,
##'  2, and 3. Also do MLE without binning and then bin into 10 bins of expected
##'  equal frequency for the GoF test.
##'
##' Simulate multiple data sets from a known individual size distribution (the
##' PLB distribution), bin them using linear bins of width 1, 5 and 10, and
##' using bins that progressively double in width, and then fit each data set
##' using the MLEbin likelihood method (like in Figures 4 and 5 of MEPS paper,
##' simulations might not be exactly the same), and then test goodness of fit
##' based on the bins, but using the number of estimated parameters (K) set to
##' (default of) 1, 2, and 3 as it's not obvious which is most appropriate. Running
##' simulation in `Goodness_of_fit.Rmd` vignette.
##' All simulated data sets have the same parameters for PLB and the same sample
##' size `n`. Individual data sets are not saved as they quickly take up a lot
##' of memory (would be `num.reps` \eqn{\times} `n` random numbers, which for the
##' default values is 10^7).
##'
##' @param n sample size of each simulated data set (numeric)
##' @param b.known known fixed value of b for all simulations
##' @param xmin.known known fixed value of xmin (minimum allowable x value); currently needs to
##'   be a power of two (since makes it simpler to define the bin widths that
##'   double in size).
##' @param xmax.known known fixed value of xmax (maximum allowable x value)
##' @param K.vec vector of values of `K` (number of parameters assumed to be
##'   estimated from the data in the goodness-of-fit test.
##' @param num.reps number of random samples to draw, where each sample is a set
##'   of `n` random numbers (like throwing `n` PLB dice `num.reps` times)
##' @param seed seed for random number generator (default is the same as for MEE paper)
##' @param binType list containing numeric values for linear bin widths and/or
##'   "2k" for bins that double in size or "unbinned" to get equal expected
##'   frequency in each bin for the goodness-of-fit test when unbinned data are
##'   available (so using MLE not MLEbin). Values
##'   other than the defaults have not yet been tested but should work.
##' @param vecDiffVal value to go into `profLike()` to compute confidence intervals.
##' @param cut.off cut-off value - data are only sampled \eqn{\geq} `cut.off`, for
##'   Figure S.37 and S.38 and Table S.5 in MEPS paper. Each resulting sample still has
##'   size `n`.
##' @param full.mult multiplier to generate desired sample size when using a
##'   `cut.off` value.
##' @param b.type either "MLE" to use the MLE, or "confMin" or "confMax" to use
##'   the minimum or maximum of the confidence interal of `b`, for the
##'   goodness-of-fit test.
##' @param epsilon tolerance for making the maximum bin `1+epsilon` times `max(x)`
##'   for simulated data that are fit using `binType = "unbinned"`, to avoid
##'   rounding errors that give error in `binData()` (as happened with second
##'   simulated data set using default settings and "unbinned").
##' @param increment.for.GoF.bins increment for dividing full range into
##'   quantiles (use 0.1 to create deciles, giving 10 groups) for the
##'   goodness-of-fit test applied to data that are not binned
##'
##' @return list containing:
##'
##'   * MLE.array.parameters: list containing values of
##'     `n`,
##'     `b.known`,
##'     `xmin.known`,
##'     `xmax.known`,
##'     `num.reps`,
##'     `binType`,
##'     `binTypes`,
##'     `binType.name`,
##'     `K.vec`,
##'     `b.type`
##'
##'   * MLE.array: two-dimensional array with element `[i, j]` representing
##'     the estimate of *b* obtained from random sample `i`, bin type `j`.
##'     Size is `num.reps` \eqn{\times} `length(binType)`. (Was
##'     three-dimensional for `MLEbin.simulate` because we had two MLE methods).
##'
##'   * MLEconf.array: three-dimensional array with vector
##'     `MLEconf.array[i, j, ]` being the confidence interval
##'     `c(confMin, confMax)` for random sample `i` and bin type `j`.
##'
##'   * GoF.list: list of length `n` with `GoF.list[[i]][[j]][[k]]` being the list result
##'   of `GoF_PLB()` for random sample `i` and bin type `j` and `K` value `K.vec[k]`.
##'
##' @export
##' @author Andrew Edwards
MLEbin.GoF.simulate = function(n = 1000,
                           b.known = -2,
                           xmin.known = 1,
                           xmax.known = 1000,
                           K.vec = c(1, 2, 3),
                           num.reps = 10000,
                           seed = 42,
                           binType = list(1, 5, 10, "2k", "unbinned"),
                           vecDiffVal = 0.5,
                           cut.off = NA,
                           full.mult = 1.5,
                           b.type = "MLE",
                           epsilon = 1e-08,
                           increment.for.GoF.bins = 0.1
                           )
{
# For debugging paste this:
# i=1; j=1; n = 1000;                           b.known = -2;                           xmin.known = 1;                           xmax.known = 1000;                           K.vec = c(1, 2, 3);                           num.reps = 10000;                           seed = 42;                           binType = "unbinned";                           vecDiffVal = 0.5;                           cut.off = NA;                           full.mult = 1.5;                           b.type = "MLE";
  if(!is.wholenumber(log2(xmin.known)) | !is.wholenumber(xmin.known))
    { stop("If want xmin.known to not be an integer and not be an integer
             power of 2 then need to edit binData; may just need to think
             about whether can just set startInteger = FALSE, but will
             need to edit binData to be able to define the binWidth for
             2k method. Will need to think about this for real data; maybe
             best to just remove the first bin (and a fit a range that is
             encompassed by the binned data).")
}

  stopifnot(b.type %in% c("MLE", "confMin", "confMax"))

  if(!is.na(cut.off))
    {
      full.sample.size = full.mult * n / (1 - pPLB(cut.off,
                                                   b = b.known,
                                                   xmin=xmin.known,
                                                   xmax=xmax.known))
                            # full.sample.size to draw from to expect to get n
                            # values >cut.off.
                            # Originally did some simulations this and then used the
                            #  values >= cut.off, which gave a distribution of
                            #  realised sample sizes, minimum of which was 858
                            #  (x.min=1, cut.off=16, b=-2). So 1000/858 = 1.165
                            #  suggests multipling full.sample.size by 1.2 should
                            #  get 1000 for each sample, so use 1.5 just to be
                            #  safe (this isn't the time consuming part of the
                            #  code). This number will change for different
                            #  parameter values so may need to be calculated
                            #  analytically, or just increase the value of full.mult.
  }

  set.seed(seed)
  binTypes = length(binType)
  binType.name = binType
  binType.name[!binType %in% c("2k", "unbinned")] =
    paste0("Linear ",
           binType[which(!binType
                         %in% c("2k", "unbinned"))])

  # Do 2-dimensional array for MLEs and then an array for confMin and confMax
  MLE.array = array(NA,
                   dim=c(num.reps, binTypes),
                   dimnames=list(1:num.reps,
                                 unlist(binType.name)))
  # No need to name the rows, just index by simulation number
  # MLE.array[i,j] is random sample i, bin type j

  # Record the confidence intervals
  MLEconf.array = array(NA,
                        dim=c(num.reps, binTypes, 2),
                        dimnames=list(1:num.reps,
                                      unlist(binType.name),
                                      c("confMin", "confMax")))
  # MLEconf.array[i,j, ] is confidence interval [c(confMin, confMax)] for
  #   random sample i, bin type j

  # Record the GoF results, just initialising the length of the list here, add
  #  each element within the loop.
#  GoF.list = vector(mode = "list",
#                    length = num.reps)
  GoF.list = list()

  # Main loop for doing the fitting num.reps times
  for(i in 1:num.reps){

    GoF.list[[i]] = list()
    if(num.reps > 1000)
      {
      if(i %in% seq(1000, num.reps, 1000)) {print(paste("i = ", i))}
                                # to show progress
      }


    if(is.na(cut.off))
      {
      x = rPLB(n,
               b = b.known,
               xmin = xmin.known,
               xmax = xmax.known)
      } else {
      x.full = rPLB(full.sample.size,
                    b = b.known,
                    xmin = xmin.known,
                    xmax = xmax.known)

      x.cut = x.full[x.full >= cut.off]

      if(length(x.cut) < n)
          {
             stop("Need to increase full.mult in MLEbin.simulate().")
          }
      x = x.cut[1:n]       # sample has desired length n
    }

    for(j in 1:binTypes)                    # Loop over binning type
      {
        GoF.list[[i]][[j]] = list()    # To initialise

        if(binType[j] != "unbinned"){     # Using MLEbin to fit binned data
          bins.list = binData(x,
                              binWidth=binType[[j]])  # 1, 2, 5 or "2k"
          num.bins = dim(bins.list$binVals)[1]

          binBreaks = bins.list$binVals[,"binMin"]$binMin   # Loses the column names
          maxOfMaxBin = bins.list$binVals[num.bins, "binMax"]$binMax
          binBreaks = c(binBreaks, maxOfMaxBin)             # Append endpoint of final bin

          binCounts = bins.list$binVals[,"binCount"]$binCount
          binMids = bins.list$binVals[,"binMid"]$binMid     # Midpoints of bins

          if(sum(!is.wholenumber(binCounts)) > 0)
          { stop("Need to adapt code for noninteger counts")
          }
          # May be needed by someone, but not yet. Though I think
          #  negLL.PLB.counts can handle it.

          sumCntLogMids = sum(binCounts * log(binMids))

          # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
          #  as a starting point for nlm for MLE of b for PLB model.
          PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts) ) - 1

          # MLEbin (maximum likelihood on binned data) calculations.

          MLEbin.res = calcLike(negLL.fn = negLL.PLB.binned,
                                p = PL.bMLE,
                                w = binBreaks,
                                d = binCounts,
                                J = length(binCounts),
                                vecDiff = vecDiffVal)

          MLE.array[i, j] = MLEbin.res$MLE
          MLEconf.array[i, j, ] = MLEbin.res$conf

          # Then do (later) the goodness-of-fit test for each value of K.vec, with the
          #  appropriate value of b
          if(b.type == "MLE")     b.use = MLEbin.res$MLE
          if(b.type == "confMin") b.use = MLEbin.res$conf[1]
          if(b.type == "confMax") b.use = MLEbin.res$conf[2]
        } else {
          stopifnot(binType[j] == "unbinned")
          # binType[j] == "unbinned"   so we want do fit unbinned using MLE
          #  then bin the data to give equal expected frequencies in each of,
          #  say, 10 bins, (though I saw this idea in S&R but couldn't find it
          #  again, is mentioned in
          #  https://stats.stackexchange.com/questions/16921/how-to-understand-degrees-of-freedom?noredirect=1&lq=1 ) based on the MLE of b. That website says:
          #  "Equal-probability binning assures the chi-squared distribution
          #   really is a good approximation to the true distribution of the
          #   chi-squared statistic about to be described.)". Be good to find a
          #   better reference, but this seems a sensible approach.

          PL.bMLE = 1/( log(min(x)) - sum(log(x))/length(x)) - 1

          MLE.res = calcLike(negLL.fn = negLL.PLB,
                             p = PL.bMLE,
                             x = x,
                             n = length(x),
                             xmin = min(x),
                             xmax = max(x),
                             sumlogx = sum(log(x)))

          MLE.array[i, j] = MLE.res$MLE
          MLEconf.array[i, j, ] = MLE.res$conf

          # Now do the goodness-of-fit test for each value of K.vec, with the
          #  appropriate value of b
          if(b.type == "MLE")     b.use = MLE.res$MLE
          if(b.type == "confMin") b.use = MLE.res$conf[1]
          if(b.type == "confMax") b.use = MLE.res$conf[2]

          # Create the bins based on the fitted distribution, as discussed above
          #  (this is now wrapped up into GoF_PLB_unbinned(), but keep like this
          #  here since doing multiple values of K later)
          binBreaks = qPLB(seq(0, 1, by = increment.for.GoF.bins),
                               b = b.use,
                               xmin = min(x),
                               xmax = max(x)
                               )
          binBreaks[length(binBreaks)] =
            binBreaks[length(binBreaks)] * (1 + epsilon)  # ensure, by
                                        # avoiding rounding errors, that
                                        # max of max bin is larger than xmax
          bins.list = binData(x,
                              binBreaks = binBreaks)  # to give same output
                                        # as above loop

          binCounts = bins.list$binVals[,"binCount"]$binCount
        }

        # Goodness-of-fit tests
        for(k in 1:length(K.vec)){
          K = K.vec[k]
          GoF_res = GoF_PLB(bin_breaks = binBreaks,
                   bin_counts = binCounts,
                   b = b.use,
                   K = K)
          GoF.list[[i]][[j]][[k]] = GoF_res
        } # End of for(k in 1:length(K.vec)) loop
      }   # End of for(j in 1:binTypes) loop
  }       # End of for(i in 1:num.reps) loop

  # Save the used parameters also:
  MLE.array.parameters = list("n" = n,
                              "b.known" = b.known,
                              "xmin.known" = xmin.known,
                              "xmax.known" = xmax.known,
                              "num.reps" = num.reps,
                              "binType" = binType,
                              "binTypes" = binTypes,
                              "binType.name" = binType.name,
                              "K.vec" = K.vec,
                              "b.type" = b.type)
  return(list("MLE.array" = MLE.array,
              "MLEconf.array" = MLEconf.array,
              "MLE.array.parameters" = MLE.array.parameters,
              "GoF.list" = GoF.list))
}

##' One figure with histogram of P-values for MLEbin, with each value of K and four binning types
##'
##' From simulation results, plot a panel of figures similar to Figure 4 of MEPS
##' paper, but showing P-value of each of 10,000 simulated data sets, binned using four
##' methods (rows), fitted using  MLEbin, and goodness of fit tested with different
##' values of K (columns). Assume K = 1, 2 and 3 for now, and four binning types.
##'
##' @param results.list output list from `MLEbin.GoF.simulate()`; see
##'   `?MLEbin.GoF.simulate` for details. TODO:
##' @param vertCol colour for vertical line at true value of `b`
##' @param vertThick thickness of vertical line at true value of `b`
##' @param xrange range of x values (should change to xlimA for consistency)
##' @param xbigticks.by increment between big tick marks on x-axis (all labelled)
##' @param xsmallticks.by increment between small tick marks on x-axis
##' @param ylimA range of y values
##' @param yBigTickLab.by increment between labelled big tick marks on y-axis
##' @param yBigTickNoLab.by increment between all big tick marks on y-axis
##' @param ySmallTick.by increment between small tick marks on y-axis
##' @param binwidth bin width for histograms; if NA then gets calculated such
##'   that `b.known` is a midpoint of a bin
##' @param cexAxis font size for axis labels
##' @param legLab character vector for labelling each panel
##' @param omi,mfrow,mai,xaxs,yaxs,mgp,cex,inset standard options for `par()`,
##'   defaults are for Figure 4, though mfrow now calculated automatically
##' @return figure with eight histograms, one for each combination of binning
##'   type and fitting method.
##' @export
##' @author Andrew Edwards
MLEbin.Pvalues.hists = function(results.list,
                              vertCol = "red",
                              vertThick = 1,
                              xrange = c(0, 1),
                              xbigticks.by = 0.2,
                              xsmallticks.by = 0.1,
                              ylimA = c(0, 1400),
                              yBigTickLab.by = 600,
                              yBigTickNoLab.by = 200,
                              ySmallTick.by = 100,
                              binwidth = NA,
                              omi = c(0.12, 0.05, 0.2, 0.0),
                              mfrow = NA,
                              mai = c(0.5, 0.5, 0.0, 0.3),
                              xaxs="i",
                              yaxs="i",
                              mgp = c(2.0, 0.5, 0),
                              cex = 0.8,
                              cexAxis = 0.9,
                              inset = c(-0.01, -0.04),
                              inset2 = c(-0.02, 0.07),
                              legLab = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                         "(f)", "(g)", "(h)", "(i)", "(j)",
                                         "(k)", "(l)", "(m)", "(n)", "(o)"))
{
  # Extract required components
#  MLE.array = results.list$MLE.array
  num.reps = results.list$MLE.array.parameters$num.reps
#  b.known = results.list$MLE.array.parameters$b.known
  binTypes = results.list$MLE.array.parameters$binTypes
  binType.name = results.list$MLE.array.parameters$binType.name
  K.vec = results.list$MLE.array.parameters$K.vec

  if(is.na(mfrow)){
    mfrow = c(binTypes, length(K.vec))
  }

  xbigticks = seq(xrange[1], xrange[2], by = xbigticks.by)
  xsmallticks = seq(xrange[1], xrange[2], by = xsmallticks.by)

  yBigTickLab = seq(0, num.reps, yBigTickLab.by)
  yBigTickNoLab = seq(0, num.reps, yBigTickNoLab.by)
  ySmallTick = seq(0, num.reps, ySmallTick.by)

  breakshist =  seq(0, 1, 0.05)

  par(omi = omi,
      mfrow = mfrow,
      mai = mai,
      xaxs = xaxs,
      yaxs = yaxs,
      mgp = mgp,
      cex = cex)

  for(binTypeInd in 1:binTypes){
    for(k in 1:length(K.vec)){
      Pvalues = vector(length = num.reps)

      for(i in 1:num.reps){
        Pvalues[i] = results.list$GoF.list[[i]][[binTypeInd]][[k]]$Pvalue
      }

      prop.above.05 = sum(Pvalues > 0.05) / length(Pvalues)

      hist(Pvalues,
           xlim=xrange,
           breaks=breakshist,
           col = c("black", rep("red", length(breakshist) - 2)),
           xlab="",
           ylab="",
           main="",
           axes=FALSE,
           ylim = ylimA)
      histAxes(yBigTickLab = yBigTickLab,
               yBigTickNoLab = yBigTickNoLab,
               ySmallTick = ySmallTick,
               cexAxis = cexAxis,
               xbigticks = xbigticks,
               xsmallticks = xsmallticks,
               vertCol = vertCol,
               vertThick = vertThick)
      legend("topright",
             paste(legLab[length(K.vec)*(binTypeInd - 1) + k], binType.name[binTypeInd]),
             bty="n",
             inset=inset)

      legend("topright", paste(round(prop.above.05*100,
                                     digits = 0),
                               "%",
                               sep=""),
             bty="n",
             inset=inset2,
             text.col = "red")

    }
  }

  mtext(expression(paste(italic(P), " value")),
        side=1,
        outer=TRUE,
        line=-1)
  mtext("Frequency",
        side=2,
        outer=TRUE,
        line=-1)

  mtext("    K=1                                      K=2                                      K=3",
        side=3,
        outer=TRUE,
        line=0)
}
