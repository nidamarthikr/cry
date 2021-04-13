#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output an CIF file
#'
#' @param filename A character string. The path to a valid CIF file.
#' @param message A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the cif file.
#' @return A named list. Each name correspond to a valid field in the SF cif.
#' @examples
#' #filename <- "/path_to_file/example.cif"
#' #lcif <- readSF_CIF(filename)
#' #print(names(lcif))
#' #print(lcif$INTRO$CELL)
#' #print(lcif$INTRO$HALL)
#' #print(lcif$INTRO$HM)
#' #print(lcif$REFL)
#' @export
readSF_CIF <- function(filename, message=FALSE){
  print("helo")
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- l_find(l_list, length(lcif))
  #l_list1 <- append(l_list,(length(lcif)-1))
  mat<-zoo::rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  crystal_summary <- lapply(ch, ucheck, pattern="_audit.revision_id ")
  symmetry <- lapply(ch, ucheck, pattern="_symmetry_equiv.id")
  reflection <- lapply(ch, ucheck, pattern="_refln.index_h")


  intro <- r_summ(lcif)
  symm <- if (is.na(nanona(symmetry)) == FALSE) clean(r_symm(nanona(symmetry))) else NULL
  reflections <- if (is.na(nanona(reflection)) == FALSE) clean(r_reflections(nanona(reflection))) else NULL
  CIF = list(HEADER=intro,SYMM=symm,REFL=reflections)
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
r_symm <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_symmetry",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_symmetry",data,value=TRUE))
  colnames(res) <- c(gsub("_symmetry_","",l_l))
  return(res)
}

r_reflections <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_refln",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_refln",data,value=TRUE))
  colnames(res) <- c(gsub("_refln.","",l_l))
  return(res)
}

r_summ <- function(x){
  data <- unlist(x)
  id <- check("_cell.entry_id",data)
  a <- reap("_cell.length_a",data)
  b <- reap("_cell.length_b",data)
  c <- reap("_cell.length_c",data)
  al <- reap("_cell.angle_alpha",data)
  be <- reap("_cell.angle_beta",data)
  ga <- reap("_cell.angle_gamma",data)
  w_i <- check("_diffrn_radiation_wavelength.id ",data)
  wl <-reap("_diffrn_radiation_wavelength.wavelength ",data)
  hr <- reap("_diffrn_reflns.pdbx_d_res_high",data)
  lr <- reap("_diffrn_reflns.pdbx_d_res_low",data)
  sg_n <- as.numeric(check("_symmetry.Int_Tables_number",data))
  sg_hall<- check1("_symmetry.space_group_name_Hall",data)
  sg_HM <- check1("_symmetry.space_group_name_H-M",data)
  cell <- list(A=a,B=b,C=c,ALPHA=al,BETA=be,GAMMA=ga)
  summ_c <- list(ENTRY=id,CELL=cell,WAVELENGTHID=w_i,WAVELENGTH=wl,HIGH_RES=hr,LOW_RES=lr,SGN=sg_n,HALL=sg_hall,HM=sg_HM)
  return(summ_c)
}
