# make figures for trap forcens comparison
# Lukas Jonkers, last edit: March 21st, 2019
library(ggplot2)
library(egg)
library(gridExtra)
library(weights)
library(reshape2)

# read plot data for results from integrated time series (LT) ####
plot.LT <- readRDS('plot_LT_HadSST_1870-1899.RDS')

# Fig 4 needs output from compare_trap_sed_publish!!

xlab1 <- expression('temperature\nchange ['*degree*C*']')
xlab2 <- expression('latitudinal\n displacement ['*10^{3}* 'km]')
ylab1 <- c('dissimilarity from\nnearest sample')
ylab2 <- c('dissimilarity from\nclosest analogue')

# Figure 2 ####
# dissimilarity to nearest and most similar sediment sample vs. historical temperature change
Fig2A <- ggplot(plot.LT, aes(x = abs(real.dT), y = sqcd.nearest, label = trap)) +
          stat_smooth(method = 'lm', mapping = aes(weight = duration), colour = 'firebrick3') +
          geom_point(aes(size = duration), alpha = 0.7) +
          coord_cartesian(ylim = c(0, 2), xlim = c(0, 1.7)) +
          ylab(ylab1) +
          xlab(xlab1) +
          scale_size_continuous(range = c(0.5, 2), guide = 'none') +
          theme_bw() +
          theme(text = element_text(size = 7),
                panel.grid = element_blank(),
                axis.title.x = element_text(vjust = 0.5))

cor.test(abs(plot.LT$real.dT), plot.LT$sqcd.nearest)
wtd.cor(abs(plot.LT$real.dT), plot.LT$sqcd.nearest, weight = plot.LT$duration)


Fig2B <- ggplot(plot.LT, aes(lat.displacement/1000)) +
          geom_histogram(breaks = c(seq(0, 2.0, by = 0.5), 2.6), fill = 'grey80', colour = 'black', lwd = 0.3) +
          ylim(c(0, 13)) +
          theme_bw() +
          xlab(xlab2) +
          theme(text = element_text(size = 7),
                panel.grid = element_blank(),
                axis.title.x = element_text(hjust = 0.5)
                )

Fig2C <- ggplot(plot.LT, aes(x = abs(real.dT), y = sqcd.min)) +
          stat_smooth(method = 'lm', mapping = aes(weight = duration), colour = 'firebrick3') +
          geom_point(aes(size = duration), alpha = 0.7) +
          coord_cartesian(ylim = c(0, 2), xlim = c(0, 1.7)) +
          ylab(ylab2) +
          xlab(xlab1) +
          scale_size_continuous(range = c(0.5, 2)) +
          theme_bw() +
          theme(text = element_text(size = 7),
                panel.grid = element_blank(),
                legend.key.height=unit(0.3, 'cm'),
                legend.box.margin = margin(0, 0, 0, -10)
                )

cor.test(abs(plot.LT$real.dT), plot.LT$sqcd.min)
wtd.cor(abs(plot.LT$real.dT), plot.LT$sqcd.min, weight = plot.LT$duration)

Fig2 <- egg::ggarrange(plots = list(Fig2A, Fig2C, Fig2B), ncol=3)

# Figure 3 ####
# histogram of counts per category
Fig3B <- ggplot(data = plot.LT) +
          geom_histogram(aes(group, colour = real.trend, fill = trap.trend), stat = 'count', lwd = 0.5, alpha = 0.8) +
          scale_fill_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
          scale_colour_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
          ylim(c(0, 30)) +
          theme_bw() +
          theme(axis.title.y = element_blank()) +
          coord_flip() +
          theme(text = element_text(size = 7),
            panel.grid = element_blank(),
            legend.position = 'right',
            rect = element_rect(fill = 'transparent'),
            legend.key.height=unit(0.3, 'cm'))

binom.test(sum(plot.LT$consistent), length(plot.LT$consistent), alternative = 'two.sided')

# Figure 4 ####
# interannual variability and sensitivity tests
# %% of years that shows change with same sign as the multiyear timeseries
Fig4A <- ggplot(interAnnual.no.impute, aes(x = N, y = consistent)) +
  geom_jitter(size = 1, height = 1, width = 0.05, colour = 'grey30', alpha = 0.5) +
  stat_summary(fun.y = "mean", geom = "point", size = 2, color = "firebrick3") +
  ylab('interannual\nconsistency [%]') +
  xlab('minimum time series length [years]') +
  coord_cartesian(ylim = c(40, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        legend.position = 'top',
        legend.box.margin = margin(c(0, 0, -15, 0)),
        legend.key.width = unit(1, 'line'),
        legend.spacing.x = unit(0.2, 'cm'),
        text = element_text(size = 7))

levels(plot.LT$deep) <- c('deeper 2 km', 'shallower 2 km')
levels(plot.LT$far) <- c('further 250 km', 'closer 250 km')
Fig4B <- ggplot(data = plot.LT) +
          geom_histogram(aes(group, colour = real.trend, fill = trap.trend), stat = 'count', lwd = 0.5, alpha = 0.8) +
          scale_fill_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
          scale_colour_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
          theme_bw() +
          theme(axis.title.x = element_blank()) +
          facet_grid(cols = vars(deep), rows = vars(far)) +
          theme(text = element_text(size = 7),
            panel.grid = element_blank(),
            strip.text = element_text(size = 5, lineheight = 0.01),
            rect = element_rect(fill = 'transparent'))

Fig4 <- ggarrange(Fig4A, Fig4B,
          nrow = 2, heights = c(1, 2))

# extended data figures ####
# estimated age of core tops (ED Fig 1) ####
# estimate mean age of sediments in mixed layer based on sediment accumulation rate and bioturbation depth
# references Berger, W. H., & Heath, G. R. (1968).
# Vertical mixing in pelagic sediments.
# Journal of Marine Research, 26(2), 134–143.
# requires:
# 1) estimate of sedimentation rate
# 2) estimate of mixing depth

# 1) 
# get estimate of sedimentation rates for ForCenS samples
# based on relationship between water depth and sed. rate
# in Holocene sediments
# Burwitz et al., 2011
# GCA
# https://doi.org/10.1016/j.gca.2011.05.029

ForCens <- readRDS('forcens_trimmed.RDS')
forcens.water.depth <- ForCens$meta$Water_depth[!is.na(ForCens$meta$Water_depth)]     
get.sed.rate <- function(z){
  (0.117/(1+(z/200)^3) + 0.006/(1+(z/4000)^10))*1000
}
sed.rate.forcens <- get.sed.rate(forcens.water.depth)
quantile(sed.rate.forcens, probs = c(0.165, 0.5, 0.835))

# 2)
# estimated global bioturbation depth
# after Boudreau 1998
# Mean mixed depth of sediments: The wherefore and the why
# L&O
mean.bio.depth <- 9.8
sd.bio.depth <- 4.5

# plot
step.size <- 0.01
ymax <- 10
xmax <- 15
sed.rate <- seq(1, ymax, by = step.size)
bio.depth <- seq(1, xmax, by = step.size)

mean.age <- outer(bio.depth, sed.rate, FUN = '/')*1000

brks <-  c(250, 500, 750, 1000, 1500, 2000, 3000, 5000)

par(pty = "s") 
plot(NA,
     xlim = range(bio.depth),
     ylim = range(sed.rate),
     xlab = 'bioturbation depth [cm]',
     ylab = 'sediment accumulation rate [cm/kyr]')
rect(xleft = mean.bio.depth - sd.bio.depth,
     xright = mean.bio.depth + sd.bio.depth,
     ybottom = quantile(sed.rate.forcens, 0.165),
     ytop = quantile(sed.rate.forcens, 0.835),
     col = 'grey80',
     border = NA) 
contour(x = bio.depth, y = sed.rate, z = mean.age, add = TRUE, levels = brks)

# example of linear regression between dissimilarity and distance (ED Fig 2)####
FigED2 <- ggplot(mod.dat[[26]], aes(x, y)) +
            geom_abline(slope = lin.mod[[26]]$coefficients, lwd = 1, colour = 'firebrick3') +
            geom_point(size = 0.7, alpha = 0.7) +
            xlab('dissimilarity') +
            ylab('latitudinal distance [km]') +
            ylim(c(0, 2500)) +
            theme_bw() +
            theme(text = element_text(size = 7),
              panel.grid = element_blank(),
              strip.text = element_text(size = 5, lineheight = 0.01),
              rect = element_rect(fill = 'transparent')) +
            annotate('text', x = 1.5, y = 0, label = names(dat.sel)[26], size = 1.5)
     
# sensitivity to size fraction (ED Fig 3) ####
plot.LT$shell_size <- factor(plot.LT$shell_size, levels = c('small', 'large'),
                             ordered = TRUE, labels=c(expression('>125 ' *mu* 'm'), expression('>150 ' *mu* 'm')))


FigED3 <- ggplot(data = plot.LT) +
  geom_histogram(aes(group, colour = real.trend, fill = trap.trend), stat = 'count', lwd = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
  scale_colour_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  facet_grid(.~shell_size, labeller = label_parsed) +
  theme(text = element_text(size = 7),
        panel.grid = element_blank(),
        strip.text = element_text(size = 5, lineheight = 0.01),
        rect = element_rect(fill = 'transparent'))

# compare HadISST with ERSSTv5 (ED Fig 4) ####
had.plot <- readRDS('plot_LT_HadSST_1870-1899.RDS')
ER.plot <- readRDS('plot_LT_ERSST_1854-1883.RDS')

ER.plot$product <- 'ERSST'
had.plot$product <- 'HadISST'

plot_data <- rbind.data.frame(ER.plot, had.plot[, -which(names(had.plot) %in% 'shell_size')])
xlab1 <- expression('temperature change ['*degree*C*']')
xlab2 <- c('latitudinal\ndisplacement [103 km]')
ylab1 <- c('dissimilarity from\nnearest sample')
ylab2 <- c('dissimilarity from\nclosest analogue')


FigS4A <- ggplot(plot_data, aes(x = abs(real.dT), y = sqcd.nearest, colour = product)) +
  geom_point(aes(size = duration), alpha = 0.7) +
  stat_smooth(method = 'lm', mapping = aes(weight = duration)) +
  coord_cartesian(ylim = c(0, 2), xlim = c(0, 1.7)) +
  ylab(ylab1) +
  xlab(xlab1) +
  scale_size_continuous(range = c(0.5, 1.5), guide = 'none') +
  scale_colour_manual(values = c('dodgerblue4', 'firebrick3'), name = 'SST product') +
  theme_bw() +
  theme(text = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = 'bottom')

FigS4B <- ggplot(data = plot_data) +
  geom_histogram(aes(group, colour = real.trend, fill = trap.trend), stat = 'count', lwd = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
  scale_colour_manual(values = c('dodgerblue4', 'firebrick3'), guide = 'none') +
  ylim(c(0, 30)) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  coord_flip() +
  theme(text = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = 'right',
        rect = element_rect(fill = 'transparent'),
        legend.key.height=unit(0.3, 'cm')) +
  facet_grid(.~product)

FigED4 <- egg::ggarrange(plots = list(FigS4A, FigS4B), ncol=2)
