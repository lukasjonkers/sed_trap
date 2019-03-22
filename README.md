# sed_trap
Data and code for "Global change drives modern plankton communities away from pre-industrial state" by Lukas Jonkers, Helmut Hillebrandt and Michal Kucera

Compare planktonic foraminifera species assemblages from sediments and sediment traps.

Scripts written by Lukas Jonkers

DATA SOURCES
* HadISST: Rayner, N. A. et al. Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century. Journal of Geophysical Research: Atmospheres 108, doi:10.1029/2002JD002670 (2003).
* ERSST v5: Huang, B. et al. NOAA Extended Reconstructed Sea Surface Temperature (ERSST), Version 5. Monthly mean. NOAA National Centers for Environmental Information. doi:10.7289/V5T72FNM. Access date: 14 Sep 2018.  (2017).
* sediment assemblages: Siccha, M. & Kucera, M. ForCenS, a curated database of planktonic foraminifera census counts in marine surface sediment samples. Scientific Data 4, 170109, doi:10.1038/sdata.2017.109 (2017).
* sediment traps: see manuscript, extended data table 1

DATA
1. Planktonic foraminifera shell flux time series\
1.1. all data: dat_sel.RDS\
1.2. time series with >125 and >150 micron data: dat_small.RDS

2. ForCenS core top sediment assemblages\
2.1. all data (excluding duplicates and samples with incomplete taxonomy): forcens_trimmed_compare.RDS\
2.2. all data, split by region: species_domains_compare.RDS\
2.3. indices of samples: domain_indeces_compare.RDS

3. SST\
3.1. average SST for each sample in ForCenS for 1870-1899 period based on HadISST: forcens_HadSST_1870-1899.RDS\
3.2. average SST for each sample in ForCenS for 1854-1883 period based on ERSST v5: forcens_ERSST_1854-1883.RDS\
3.3. average SST for each sediment trap site for 1870-1899 period based on HadISST: traps_HadSST_1870-1899.RDS\
3.4. average SST for each sediment trap site for 1854-1883 period based on ERSST v5: traps_ERSST_1854-1883.RDS\
3.5. average SST for each sediment trap site for deployment period based on HadISST: traps_HadSST_period.RDS\
3.6. average SST for each sediment trap site for deployment period based on ERSST v5: traps_ERSST_period.RDS\
3.7. linear SST trend between 1870 and 2015 based on HadISST: hadisst_trend_1870-2015.RDS\
3.8. average SST for period of sediment trap observations: hadisst_mean_1978-2013.RDS

CODE
1. get_ForCenS.R: selection of ForCenS data. Used to generate data 2.1-2.3
2. make_polygons.R: make circles with 100 km radius around trap and core tope sites used in 4 and 5
3. make_polygon_function.R: used by 2
4. extract_HadSST.R: extraction of data 3.1, 3.3, 3.5, 3.7, 3.8
5. extract_ERRSTv5.R: extraction of data 3.2, 3.4, 3.6
6. make_annual_assemblages.R: process shell flux time series (data 1.1 and 1.2)
7. make_annual_fluxes_function.R: used by 6
8. compare_trap_sed_publish.R: code to compare sediment trap and core top assemblages
9. compare_figs.R: code to create figures
10. maps_robinson.R: code to create maps
