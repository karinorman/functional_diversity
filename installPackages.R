# Install pacman if it isn't already installed

if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")


# Install analysis packages using pacman

pacman::p_load(dplyr, ggplot2, tidyr, stringr, sf, 
               FD, spData, sf, tmap)