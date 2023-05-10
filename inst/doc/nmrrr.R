## ----rmd-setup, include = FALSE-----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  error = FALSE,
  comment = "#>"
)

## ----packages-----------------------------------------------------------------
library(nmrrr)

library(ggplot2)
theme_set(theme_bw()) # set the default ggplot theme

## ----input-directories--------------------------------------------------------

SPECTRA_FILES <- system.file("extdata", "kfp_hysteresis", "spectra_mnova", package = "nmrrr")
PEAKS_FILES <- system.file("extdata", "kfp_hysteresis", "peaks_mnova_multiple", package = "nmrrr")


## ----spectra-import-----------------------------------------------------------
spectra_df <- nmr_import_spectra(path = SPECTRA_FILES,
                                 method = "mnova") 

str(spectra_df)

## ----spectra-process----------------------------------------------------------
spectra_df <- subset(spectra_df, ppm >= 0 & ppm <= 10)

str(spectra_df)

## ----spectra-plot, fig.height=5, fig.width=8----------------------------------
nmr_plot_spectra(dat = spectra_df,
                 binset = bins_Clemente2012,
                 label_position = 5,
                 mapping = aes(x = ppm, 
                               y = intensity, 
                               group = sampleID, 
                               color = sampleID),
                 stagger = 0.5) +
  # OPTIONAL PARAMETERS/LAYERS
  geom_rect(aes(xmin = 2, xmax = 4, ymin = 0, ymax = 5.5), 
            fill = "white", color = NA, alpha = 0.8)+
  labs(subtitle = "binset: Clemente et al. 2012")+
  ylim(0, 5.5)

## -----------------------------------------------------------------------------
spectra_bins <- nmr_assign_bins(dat = spectra_df,
                                binset = bins_Clemente2012)

## -----------------------------------------------------------------------------
spectra_bins <- subset(spectra_bins, group != "oalkyl")

## -----------------------------------------------------------------------------
relabund_integration <- nmr_relabund(dat = spectra_bins,
                                     method = "AUC")

## -----------------------------------------------------------------------------
peaks_df <- nmr_import_peaks(path = PEAKS_FILES, 
                             method = "multiple columns")

str(peaks_df)

## -----------------------------------------------------------------------------
peaks_df <- subset(peaks_df, Type == "Compound")

## -----------------------------------------------------------------------------
peaks_bins <- nmr_assign_bins(dat = peaks_df,
                              binset = bins_Clemente2012)

## -----------------------------------------------------------------------------
peaks_bins <- subset(peaks_bins, group != "oalkyl")

## -----------------------------------------------------------------------------
relabund_peaks <- nmr_relabund(dat = peaks_bins,
                               method = "peaks")

## ---- fig.height=5, fig.width=8-----------------------------------------------

ggplot(relabund_integration,
       aes(x = sampleID, y = relabund, fill = group))+
  geom_bar(stat = "identity")+
  labs(title = "Relative abundance by AUC",
       subtitle = "binset: Clemente et al. 2012")

## ----example_2, eval=FALSE----------------------------------------------------
#  library(tidyverse)
#  
#  SS_FILE <- system.file("extdata", "amp_burnseverity", "spectra_wide.csv", package = "nmrrr")
#  ss_data <- read.csv(SS_FILE)
#  
#  ## Make long form and do additional cleaning if needed.
#  ss_data_long =
#    ss_data %>%
#    pivot_longer(-ppm,
#                 names_to = "sampleID",
#                 values_to = "intensity") %>%
#    arrange(sampleID, ppm)
#  
#  ss_data_long = subset(ss_data_long, ppm >= 0 & ppm <= 250)
#  ss_data_long = subset(ss_data_long, intensity >= 0)
#  
#  ## Assign bins
#  data_long_bins = nmr_assign_bins(dat = ss_data_long,
#                                   binset = bins_ss_Clemente2012)
#  
#  ## Plot spectra
#  nmr_plot_spectra(dat = data_long_bins,
#                   binset = bins_ss_Clemente2012,
#                   mapping = aes(x = ppm, y = intensity,
#                                 group = sampleID,
#                                 color = sampleID),
#                   stagger = 15,
#                   label_position = 70)+
#    theme(axis.text.y = element_blank())+
#    xlim(210, 0)
#  
#  ## Calculate relative abundance
#  data_relabund = nmr_relabund(dat = data_long_bins,
#                               method = "AUC")
#  
#  ggplot(data = data_relabund,
#         aes(x =  sampleID,
#             y = relabund,
#             fill = group))+
#    geom_bar(stat = "identity")
#  

## -----------------------------------------------------------------------------

bins_Clemente2012


