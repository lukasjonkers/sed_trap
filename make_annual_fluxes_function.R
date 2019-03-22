# impute sediment trap time series. Fill gaps with median flux of other years
# written by Lukas Jonkers
# last edit 04 June 2018
process.TS <- function(x, it){ # x = flux time series, it = imputatation threshold: fraction of year
  if(!require(lubridate)){install.packages('lubridate')}else{library(lubridate)}
  if(!require(plyr)){install.packages('plyr')}else{library(plyr)}
  dat <- x$dat
  dat <- dat[complete.cases(dat),]
  spp.flux <- dat[,-c(1,2)]
  # rename spp names
  rename.fun <- function(name){
    if(name == 'G_rubescens'){
      name <- 'G_rub1'
    }else if(name == 'G_conglobatus'){
      name <- 'G_con1'
    }else{
      name <- substr(name, start = 1, stop = 5)
    }
    name
  }
  names(spp.flux) <- sapply(names(spp.flux), rename.fun)
  open <- as.Date(dat$open)
  close <- as.Date(dat$close)
  TS.length <- as.numeric(difftime(tail(dat$close, 1), dat$open[1], units = 'days'))
  TS.years <- TS.length/365.25
  TS.whole.years <- floor(TS.years)
  TS.excess <- TS.years - TS.whole.years
  # gaps of one day are accepted because open and close are not always on same day
  gap.length <- as.numeric(difftime(dat$open[-1], dat$close[-length(dat$close)], units = 'days'))
  has.gaps <- any(gap.length >1)
  
  # make daily time series 
  collection.interval <- as.numeric(difftime(dat$close, dat$open, units = 'days'))
  # check if collection intervals OK
  if(!sum(collection.interval) + sum(gap.length) == TS.length){print('open and close dates one day apart, time missing')}
  
  TS.daily <- spp.flux[rep(seq_len(nrow(spp.flux)), collection.interval),]
  
  days.date <- do.call(c, mapply(function(i, j) i + 1:j-1, i = open, j = collection.interval, SIMPLIFY = FALSE))
  days.from.start <- as.numeric(days.date - days.date[1])
  yday.date <- yday(days.date)
  ydays.from.start <- yday(days.date - days(days.date[1]))
  
  # get annual fluxes
  if(TS.length <= 366){
    if(has.gaps){
      # linearly fill gaps
      desired.days <- seq(min(days.from.start), max(days.from.start))
      tem <- apply(TS.daily, 2, function(i) approx(days.from.start, i, xout = desired.days, method = 'linear')$y)
      annual.flux <- list(meta = x$meta,
                          flux = list(cbind.data.frame(days.from.start = desired.days, tem)),
                          year = year(days.date[1] + 183),
                          timeseries.length = TS.length)
    }
    else{
      # use as is
      annual.flux <- list(meta = x$meta,
                          flux = list(cbind.data.frame(days.from.start, TS.daily)),
                          year = year(days.date[1] + 183),
                          timeseries.length = TS.length)
    }
  }else{
    # make median year, interpolate gaps if necessary
    df <- cbind.data.frame(day = ydays.from.start, TS.daily)
    median.year <- ddply(df, .(day), colwise(median))
    # if length of median year is not 366, add first day to end
    if(length(median.year$day) != 366){
      median.year <- rbind(median.year, head(median.year, 1))
      median.year$day[length(median.year$day)] <- 367
    }
    median.year.smooth <- cbind.data.frame(day = 1:366,
        apply(median.year[,-1], 2, function(j) predict(loess(j ~ median.year$day, span = 0.1), newdata = 1:366)))
    median.year.smooth[median.year.smooth <0] <- 0
    
    
    # make daily time series without gaps for whole years + 1
    start.gap <- which(diff(days.from.start) >2)
    gap.len <- gap.length[gap.length >1]
    gaps <- mapply(function(i, j) as.vector(seq(from = i+1, length.out = j)), i = ydays.from.start[start.gap], j = gap.len)
    # and for missing end
    if(!is.list(gaps)){gaps <- list(gaps)}
    gaps[[length(gaps)+1]] <- (tail(ydays.from.start, 1)+1):365
    # in case gaps go across year boundary
    gaps <- lapply(gaps, function(i){
      if(any(i>365)){
        i[i>365] <- i[i>365]-365
        i
        }else{i}
    })
    # get fluxes for gaps
    gap.fluxes <- lapply(gaps, function(i) median.year.smooth[i,])
    # add impute flag
    gap.fluxes <- lapply(gap.fluxes, function(i) cbind.data.frame(i, imputed = rep(TRUE, times = nrow(i))))
    dfw <- cbind.data.frame(df, imputed = rep(FALSE, times = nrow(df)))
    # split and rbind imputed gaps
    TS.cut <- split(dfw, cumsum(1:nrow(dfw) %in% (start.gap+1)))
    TS.imputed <- do.call(rbind, mapply(rbind, TS.cut, gap.fluxes, SIMPLIFY = FALSE))
    # split into years
    #split.indx <- which(TS.imputed$day == 1)[-1]
    split.indx <- which(diff(TS.imputed$day) < 0) +1
    TS.imputed.yearly <- split(TS.imputed, cumsum(1:nrow(TS.imputed) %in% (split.indx)))
    #year.mid <- year(days.date[ydays.from.start == 1] +183)
    year.mid <- seq( from = year(days.date[1] +183), to = year(days.date[1] +183) + TS.whole.years)
    if(TS.excess <= it){
        annual.flux <- list(meta = x$meta,
                            flux = TS.imputed.yearly[1:TS.whole.years],
                            year = year.mid[1:TS.whole.years],
                            timeseries.length = TS.length)
    }
    else{
      annual.flux <- list(meta = x$meta,
                          flux = TS.imputed.yearly,
                          year = year.mid,
                          timeseries.length = TS.length)
    } 
    
  }
  annual.flux
}

