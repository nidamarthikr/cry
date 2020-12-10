readpd_rtv <- function(filename, messages=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  c_lcif <- lcif[!grepl('^#|^.*#', lcif)]
  l_list <- grep("loop_",c_lcif)
  l_list1 <- l_find(l_list, length(c_lcif))
  mat<-zoo::rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) c_lcif[(x[1]+1):(x[2]-1)])
  diffraction <- lapply(ch, ucheck, pattern="_pd_proc_point_id")
  reflection <- lapply(ch, ucheck, pattern="_refln_index_h")

  intro <- r_summ(c_lcif)
  diffractions <- if (is.na(nanona(diffraction)) == FALSE) clean(r_diffractions(nanona(diffraction))) else NULL
  reflections <- if (is.na(nanona(reflection)) == FALSE) clean(r_peakreflns(nanona(reflection))) else NULL
  CIF = list(HEADER=intro,DIFF=diffractions,REFL_PEAK=reflections)
  close(f)
  return(CIF)
}



### accessory functions ####

l_find <- function(a,n){
  if (length(a) > 1){
    a1 <- append(a,(n-1))
  } else
  { a1 <- c(2,a,(n-1))
  return(a1)
  }
}

stack<-function(x){
  j <- c(x[1],x[2])
  return(j)
}

ucheck <- function(x,pattern){
  r <- unlist(x)
  if (length(grep(pattern,r))>0){
    piece <- r
  } else
  { piece <- NA
  return(piece)
  }
}

ansnull <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

nanona <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NA
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

recheck <- function(r1){
  r2 <- gsub("[:):]","",gsub("[:(:]",",",r1))
  r3 <- as.numeric((strsplit(r2, ",")[[1]])[1])
  return(r3)
}

check <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  r2 <- if (length(grep("[:(:]",r1,value = TRUE)>0) == TRUE) recheck(r1) else r1
  return(r2)
}

check1 <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "'")[[1]])[2] else NA
  return(r1)
}

clean1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- as.data.frame(x)
  return(out)
  }
}


clean <- function(x){
  co1 <- data.frame(gsub ("[()]","",as.matrix(x),perl=T))
  ref <- data.frame(gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(x),perl=T))
  ref1 <- data.frame(gsub("[()]","",as.matrix(ref),perl=T))
  ref1[ref1==""]<-NA
  ref2 <- clean1(ref1)
  return(list(VAL=co1,STD=ref2))
}

reap1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- as.numeric(x)
  return(out)
  }
}

reap <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  v <- as.numeric(gsub ("[()]","",as.matrix(r1),perl=T))
  s <- gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(r1),perl=T)
  s1 <- gsub("[()]","",as.matrix(s),perl=T)
  s2 <- reap1(s1)
  return(list(VAL=v,STD=s2))
}


r_diffractions <- function (x){
  data <- unlist(x)
  data <- data[!grepl("_number_of_points", data)]
  nskip <- length((grep("_pd",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_pd_",data,value=TRUE))
  colnames(res) <- c(gsub("_pd_","",l_l))
  return(res)
}

r_peakreflns <- function (x){
  data <- unlist(x)
  data <- data[!is.na(data)]
  nskip <- length((grep("_refln",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_refln",data,value=TRUE))
  colnames(res) <- c(gsub("_refln_","",l_l))
  return(res)
}

r_summ <- function(x){
  data <- unlist(x)
  id <- check("_pd_block_id",data)
  pdm <- check("_pd_calc_method",data)
  mtmi <- reap("_meas_2theta_range_min",data)
  mtma <- reap("_meas_2theta_range_max",data)
  mti <- reap("_meas_2theta_range_inc",data)
  ptmi <- reap("_proc_2theta_range_min",data)
  ptma <- reap("_proc_2theta_range_max",data)
  pti <- reap("_proc_2theta_range_inc",data)
  rp <- check("_diffrn_radiation_probe ",data)
  wl <-reap("_diffrn_radiation_wavelength ",data)
  rf <- reap("_proc_ls_prof_R_factor",data)
  wr <- reap("_proc_ls_prof_wR_factor",data)
  e_wr <- reap("_proc_ls_prof_wR_expected",data)
  Theta_r <- list(Meas_min=mtmi,Meas_max=mtma,Meas_inc=mti,Proc_min=ptmi,Proc_max=ptma,Proc_inc=pti)
  R_fit <- list(Rp=rf,Rwp=wr,Rexp=e_wr)
  summ_c <- list(ENTRY=id,THETA_RANGE=Theta_r,PROBE_TYPE=rp,WAVELENGTH=wl,R_PRO_FIT=R_fit)
  return(summ_c)
}
