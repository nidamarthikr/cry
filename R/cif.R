#
# This file is part of the cry package
#
readcif <- function(filename, messages=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- append(l_list,length(lcif))
  mat<-rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  crystal_summary <- lapply(ch, ucheck, pattern="_publ_author_name")
  symmetry <- lapply(ch, ucheck, pattern="_space_group_symop_operation")
  reflection <- lapply(ch, ucheck, pattern="_refln.index_h")
  coordinates <- lapply(ch, ucheck, pattern="_atom_site_label")
  anisotropy <- lapply(ch, ucheck, pattern="_atom_site_aniso_label")
  symbol <- lapply(ch, ucheck, pattern="_atom_type_symbol")
  #geom_angle <- lapply(ch, ucheck, pattern="geom_angle_atom_site_label_1")
  #geom_distance <- lapply(ch, ucheck, pattern="_geom_bond_atom_site_label_1")
  #geom_hbond <- lapply(ch, ucheck, pattern="geom_hbond_atom_site_label_D")
  #geom_torsion <- lapply(ch, ucheck, pattern="geom_torsion_atom_site_label_1")

  intro <- r_summ(lcif)
  symm <- nanona(symmetry)
  reflections <- nanona(reflection)
  coordinate <- r_positions(nanona(coordinates))
  anisotropies <- r_aniso(nanona(anisotropy))
  symbolics <- nanona(symbol)
  #angle <- r_angle(geom_angle)
  #distance <- r_dist(geom_distance)
  #hbond <- r_hbond(geom_hbond)
  #torsion <- r_torsion(geom_torsion)
  CIF=list(INTRO=intro,SYMM=symm,COOR=coordinate,ANISO=anisotropies)
  #CIF = list(INTRO=intro,SYMM=symm,REFL=reflections,COOR=coordinate,ANISO=anisotropies,SYMB=symbolics,ANGLE=angle,DIST=distance,HBOND=hbond,TORSION=torsion)
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

nanona <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NA
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

check <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  return(r1)
}

check1 <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "'")[[1]])[2] else NA
  return(r1)
}

r_positions <- function (x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=4,col.names=c("label","x","y","z")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_aniso <- function(x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=7,col.names=c("label","U_11","U_22","U_33","U_12","U_13","U_23")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}


r_angle <- function(x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=4,col.names=c("label_1","label_2","label_3","angle")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_dist <- function(x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=3,col.names=c("label_1","label_2","distance")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_hbond <- function(x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=8,col.names=c("D","H","A","symmetry_A","DH","HA","DA","DHA")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_torsion <- function(x){
  data <- unlist(x)
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=5,col.names=c("label_1","label_2","label_3","label_4","torsion")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
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
  v <- check("_cell_volume",data)
  den <-check("_exptl_crystal_density_diffrn",data)
  c_sy <- check("_symmetry_cell_setting",data)
  sg_n <- check("_space_group_IT_number",data)
  sg_hall<- check1("_symmetry_space_group_name_Hall",data)
  sg_HM <- check1("_symmetry_space_group_name_H-M",data)
  cell <- as.numeric(list(a,b,c,al,be,ga))
  sg <- list(c_sy,sg_n,sg_hall,sg_HM)
  prop <- list(v,den)
  summ_c <- list(CELL=cell,SG=sg,PROP=prop)
  return(summ_c)
}


