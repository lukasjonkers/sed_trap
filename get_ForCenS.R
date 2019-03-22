# prepare ForCenS data
# Lukas Jonkers, last edit 22 Mar 2019

# Data citation:
# Siccha, M; Kucera, M (2017): ForCenS, a curated database of planktonic foraminifera census counts in marine surface sediment samples.
# Scientific Data, 4, 170109, https://doi.org/10.1038/sdata.2017.109

# Data link:
# https://doi.pangaea.de/10.1594/PANGAEA.873570

# this script removes species groups, morphospecies and unIDs
# also removes duplicates and simplifies the species names
# and replaces NAs with 0
# the second part splits the dataset by ocean basin as per Kucera et al., (2005): Reconstruction of sea-surface temperatures
# from assemblages of planktonic foraminifera: multi-technique approach based on geographically constrained calibration
# data sets and its application to glacial Atlantic and Pacific Oceans, Quaternary Science Reviews, 24, 951-998, 2005.
# New data from the Red Sea are added to the Indian Ocean domain

library(sp)

# munge ForCenS data ####
# download ForCenS from https://doi.pangaea.de/10.1594/PANGAEA.873570
# load data, make sure to delete 2nd row manually before import
ForCenS_raw <- read.delim('ForCenS.txt', na.strings = 'N/A')

# remove rows with species groups
rm.rows <- unique(c(	
  which(!is.na(ForCenS_raw$Globorotalia_menardii_._Globorotalia_tumida)),
  which(!is.na(ForCenS_raw$Globigerinoides_ruber_._Globigerinoides_white)),
  which(!is.na(ForCenS_raw$Turborotalita_humilis_._Berggrenia_pumilio))))

# remove columns with morphospecies and unIDs
rm.cols <- c(
  which(names(ForCenS_raw) == 'Globorotalia_menardii_._Globorotalia_tumida'),
  which(names(ForCenS_raw) == 'Globigerinoides_ruber_._Globigerinoides_white'),
  which(names(ForCenS_raw) == 'Turborotalita_humilis_._Berggrenia_pumilio'),
  which(names(ForCenS_raw) == 'Globorotalia_truncatulinoides_dextral_coiling'),
  which(names(ForCenS_raw) == 'Globorotalia_truncatulinoides_sinistral_coiling'),
  which(names(ForCenS_raw) == 'Trilobatus_sacculifer_w_sac_chamber'),
  which(names(ForCenS_raw) == 'Trilobatus_sacculifer_wo_sac_chamber'),
  which(names(ForCenS_raw) == 'Turborotalita_quinqueloba_dextral_coiling'),
  which(names(ForCenS_raw) == 'Turborotalita_quinqueloba_sinistral_coiling'),
  which(names(ForCenS_raw) == 'unidentified'))

dat <- ForCenS_raw[-rm.rows,-rm.cols]

# remove replicates by randomly selecting one sample from replicates (determined only by location)
dat <- dat[sample(nrow(dat)), ]
vars <- c('Latitude', 'Longitude')
dat <- dat[-which(duplicated(dat[,vars])), ]

# check
anyDuplicated(dat[,vars])

# split
meta <- dat[, 1:which(names(dat) == 'Count')]
species <- dat[, -(1:which(names(dat) == 'Count'))]
# normalise species
species <- sweep(species, 1, rowSums(species, na.rm = TRUE), FUN = '/')
# check if all samples have an ocean descriptor, you may have to go back to the text file to enter for a single sample
which(is.na(meta$Ocean))

# simplify species names
original_name <- names(species)
short_name <- strsplit(original_name, '_')
short_name <- make.unique(sapply(short_name, function(x){
  paste0(substr(x[1], 1, 1), '_', substr(x[2], 1, 3))
}), sep = '')
names(species) <- short_name

species_names <- cbind.data.frame(original_name, short_name)

# replace NAs in species df with 0s
species <- replace(species, is.na(species), 0)

# make a list
ForCenS_trim <- list(meta = meta, species = species, species_names = species_names)

saveRDS(ForCenS_trim, file = 'forcens_trimmed_compare.RDS')

# split by ocean/region ####

# indices for each region
# North Atlantic and Arctic
NAT <- which(ForCenS_trim$meta$Ocean == 7 | ForCenS_trim$meta$Ocean == 513)
# South Atlantic including samples in North Atlantic up to 10 N and Atlantic Sector of Southern Ocean
SAT <- c(which(ForCenS_trim$meta$Ocean == 11),
         which(ForCenS_trim$meta$Ocean == 7 & ForCenS_trim$meta$Latitude < 10),
         which(ForCenS_trim$meta$Ocean == 257 & ForCenS_trim$meta$Longitude > -60 & ForCenS_trim$meta$Longitude < 20))
# Mediterranean including NAtl (following Hayes et al., 2005) for LGM reconstructions
# North Atlantic extent of Med training set (Hayes et al., 2005)
MED_NAT_poly <- list(
  x = c(-13, -23, -30, -30, 5, 5, -0, 0, -13, -13),
  y = c(33, 25, 25, 70, 70, 63, 59, 42, 42, 33)				
)
# within polygon and N. pachyderma <0.7
MED_NAT_indx <- which(point.in.polygon(ForCenS_trim$meta$Longitude, ForCenS_trim$meta$Latitude, pol.x = MED_NAT_poly$x, pol.y = MED_NAT_poly$y) > 0 & ForCenS_trim$species$N_pac <0.7) 
MDX <- c(which(ForCenS_trim$meta$Ocean == 1025), MED_NAT_indx)
# Indian Ocean, including Southern Ocean and (new) Red Sea
IND <- c(which(ForCenS_trim$meta$Ocean == 129 | ForCenS_trim$meta$Ocean == 2049),
         which(ForCenS_trim$meta$Ocean == 257 & ForCenS_trim$meta$Longitude < 150 & ForCenS_trim$meta$Longitude >20))
# Pacific including off west Australia
PAC <- c(which(ForCenS_trim$meta$Ocean == 49 | ForCenS_trim$meta$Ocean == 81),
         which(ForCenS_trim$meta$Ocean == 129 & ForCenS_trim$meta$Longitude > 105),
         which(ForCenS_trim$meta$Ocean == 257 & ForCenS_trim$meta$Longitude > 105),
         which(ForCenS_trim$meta$Ocean == 257 & ForCenS_trim$meta$Longitude < -62))

domain.indeces <- list(
  NAT = NAT,
  SAT = SAT,
  MDX = MDX,
  IND = IND,
  PAC = PAC
)

saveRDS(domain.indeces, file = 'domain_indeces_compare.RDS')

species_domains <- list(
  NAT = list(meta = ForCenS_trim$meta[NAT,], species = ForCenS_trim$species[NAT, ]),
  SAT = list(meta = ForCenS_trim$meta[SAT,], species = ForCenS_trim$species[SAT, ]),
  MDX = list(meta = ForCenS_trim$meta[MDX,], species = ForCenS_trim$species[MDX, ]),
  IND = list(meta = ForCenS_trim$meta[IND,], species = ForCenS_trim$species[IND, ]),
  PAC =  list(meta = ForCenS_trim$meta[PAC,], species = ForCenS_trim$species[PAC, ]),
  species_names = ForCenS_trim$species_names
)

# remove species that do not occur in a region
spec.zero <- lapply(species_domains[1:5], function(x) which(colSums(x$species) <= 1e-08))

for(i in 1:5){
  if(length(spec.zero[[i]]) >0){
    species_domains[[i]]$species <- species_domains[[i]]$species[ ,-spec.zero[[i]]]
  }
}

saveRDS(species_domains, 'species_domains_compare.RDS')
rm(list = ls())
