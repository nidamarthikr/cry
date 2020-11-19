#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output an CIF file
#'
#' @param filename A character string. The path to a valid CIF file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the cif file.
#' @return A named list. Each name correspond to a valid field in the cif.
#' @examples
#' filename <- "/path_to_file/example.cif"
#' lcif <- readCIF(filename)
#' print(names(lcif))
#' print(lcif$INTRO$CELL)
#' print(lcif$INTRO$HALL)
#' print(lcif$INTRO$HM)
#' print(lcif$SYMM)
#' @export
readCIF <- function(filename, messages=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- append(l_list,length(lcif))
  mat<-rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  crystal_summary <- lapply(ch, ucheck, pattern="_publ_author_name")
  symmetry <- lapply(ch, ucheck, pattern="_space_group_symop_operation")
  reflection <- lapply(ch, ucheck, pattern="_refln.index_h")
  coordinates <- lapply(ch, ucheck, pattern="_atom_site_fract_x")
  anisotropy <- lapply(ch, ucheck, pattern="_atom_site_aniso_label")
  symbol <- lapply(ch, ucheck, pattern="_atom_type_symbol")
  geom_angle <- lapply(ch, ucheck, pattern="geom_angle_atom_site_label_1")
  geom_distance <- lapply(ch, ucheck, pattern="_geom_bond_atom_site_label_1")
  geom_hbond <- lapply(ch, ucheck, pattern="geom_hbond_atom_site_label_D")
  geom_torsion <- lapply(ch, ucheck, pattern="geom_torsion_atom_site_label_1")

  intro <- r_summ(lcif)
  symm <- ansnull(symmetry)
  #reflections <- if (is.na(nanona(reflection)) == FALSE) r_angle(reflection) else NULL
  reflections <- ansnull(reflection)
  coordinate <- if (is.na(nanona(coordinates)) == FALSE) r_positions(nanona(coordinates)) else NULL
  #coordinate <- r_positions(ansnull(coordinates))
  anisotropies <- if (is.na(nanona(anisotropy)) == FALSE) r_aniso(nanona(anisotropy)) else NULL
  #anisotropies <- r_aniso(ansnull(anisotropy))
  symbolics <- ansnull(symbol)
  angle <- if (is.na(nanona(geom_angle)) == FALSE) r_angle(nanona(geom_angle)) else NULL
  #angle <- r_angle(ansnull(geom_angle))
  distance <- if (is.na(nanona(geom_distance)) == FALSE) r_dist(nanona(geom_distance)) else NULL
  #distance <- r_dist(ansnull(geom_distance))
  hbond <- if (is.na(nanona(geom_hbond)) == FALSE) r_hbond(nanona(geom_hbond)) else NULL
  #hbond <- r_hbond(ansnull(geom_hbond))
  torsion <- if (is.na(nanona(geom_torsion)) == FALSE) r_torsion(nanona(geom_torsion)) else NULL
  #torsion <- r_torsion(ansnull(geom_torsion))
  #CIF=list(INTRO=intro,SYMM=symm,COOR=coordinate,ANISO=anisotropies)
  CIF = list(HEADER=intro,SYMM=symm,REFL=reflections,COOR=coordinate,ANISO=anisotropies,SYMB=symbolics,ANGLE=angle,DIST=distance,HBOND=hbond,TORSION=torsion)
  close(f)
  return(CIF)
}


### accessory functions ####

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

r_positions <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_atom",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_atom",data,value=TRUE))
  colnames(res) <- c(gsub("_atom_site_","",l_l))
  return(res)
}

r_aniso <- function(x){
  data <- unlist(x)
  nskip <- length((grep("_atom",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_atom",data,value=TRUE))
  colnames(res) <- c(gsub("_atom_site_aniso_","",l_l))
  return(res)
}


r_angle <- function(x){
  data <- unlist(x)
  nskip <- length((grep("_geom_angle",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_geom_angle",data,value=TRUE))
  colnames(res) <- c(gsub("_geom_angle_","",l_l))
  return(res)
}

r_dist <- function(x){
  data <- unlist(x)
  nskip <- length((grep("_geom_bond",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_geom_bond",data,value=TRUE))
  colnames(res) <- c(gsub("_geom_bond_","",l_l))
  return(res)
}

r_hbond <- function(x){
  data <- unlist(x)
  nskip <- length((grep("_geom_hbond",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_geom_hbond",data,value=TRUE))
  colnames(res) <- c(gsub("_geom_hbond_","",l_l))
  return(res)
}

r_torsion <- function(x){
  data <- unlist(x)
  nskip <- length((grep("_geom_torsion",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_geom_torsion",data,value=TRUE))
  colnames(res) <- c(gsub("_geom_torsion_","",l_l))
  return(res)
}



r_summ <- function(x){
  data <- unlist(x)
  formula <- check1("_chemical_formula_sum",data)
  a <- check("_cell_length_a",data)
  b <- check("_cell_length_b",data)
  c <- check("_cell_length_c",data)
  al <- check("_cell_angle_alpha",data)
  be <- check("_cell_angle_beta",data)
  ga <- check("_cell_angle_gamma",data)
  v <- as.numeric(check("_cell_volume",data))
  den <-as.numeric(check("_exptl_crystal_density_diffrn",data))
  c_sy <- check("_symmetry_cell_setting",data)
  sg_n <- as.numeric(check("_space_group_IT_number",data))
  sg_hall<- check1("_symmetry_space_group_name_Hall",data)
  sg_HM <- check1("_symmetry_space_group_name_H-M",data)
  cell <- as.numeric(list(a,b,c,al,be,ga))
  #sg <- list(c_sy,sg_n,sg_hall,sg_HM)
  #prop <- list(v,den)
  summ_c <- list(FORMULA=formula,CELL=cell,VOL=v,DEN=den,CrysSys=c_sy,SGN=sg_n,HALL=sg_hall,HM=sg_HM)
  return(summ_c)
}


