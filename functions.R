library(tidyverse)
library(readxl)
library(tools)
library(reshape2)
library(pheatmap)
library(viridis)
library(here)
library(scales)
library(grid)
library(ggthemes)
library(plater)
library(networkD3)


import_all <- function(input_dir) {
  files = list.files(path = input_dir, pattern = "*.xlsx")
  
  ## load and bind rows of all files ##
  big_tbl <- file.path(input_dir, files) %>% map_df(~import_wb(.))
  
  return(big_tbl)
}

import_wb <- function(input_file) {
  # The consolidated table
  wb <- tibble()
  
  # Create a data frame of sheet names in the workbook with columns for day and wether it's dead or live
  sheets <- excel_sheets(input_file)
  mdata_d <- matrix(ncol = 2, unlist(strsplit(sheets, '_')), byrow = TRUE)
  mdata <-
    data.frame(
      sheet = sheets,
      Day = as.integer(gsub('^D', '', mdata_d[, 1])),
      LD = mdata_d[, 2],
      stringsAsFactors = FALSE
    )
  
  # Iterate across all worksheets in the workbook and add them to the consolidated table
  # For every sheet, two columns are added: one for the day and one for the Live/dead status
  for (i in seq(from = 1, to = nrow(mdata))) {
    tws <- read_excel(input_file,
                 sheet = as.character(mdata[i, "sheet"]),
                 col_types = 'text')
    names(tws) <- gsub(" ?(?:Live|Dead) ?", '', names(tws), perl = TRUE) #Some column names are different in the _Dead and _Alive sheets. This removes references to Live and Dead
    tws$Day <- mdata[i, "Day"]
    tws$LD <- mdata[i, "LD"]
    wb <- rbind(wb, tws)
  }
  names(wb) <- gsub("^- ", "", names(wb))
  wb$ExpID <- basename(tools::file_path_sans_ext(input_file))
  return(wb)
}

clean_input <- function(wb) {
  #Fixing known inconsitencies
  wb$`Cell Type` <- gsub(
    pattern = " CRC",
    replacement = "",
    x = wb$`Cell Type`)
  wb$Compound[is.na(wb$Compound)] <- "Media"
  wb$Compound[wb$Compound == 'control'] <- "Control"
  wb$Compound[wb$Compound == 'Control'] <- "Media"
  wb$Compound[wb$Compound == 'media'] <- "Media"
  wb[wb$Compound == 'STA 5' & is.null(wb$Concentration == 'STA 5'),"Concentration"] <- 5
  wb[wb$Compound == 'STA 5' & is.null(wb$Concentration == 'STA 5'),"Compound"] <- "STA"
  
  
  #Setting types
  wb <- readr::type_convert(wb)
  
  wb$Concentration <- as.factor(wb$Concentration)
  wb$`Cell Type` <- as.factor(wb$`Cell Type`)
  wb$Compound <- as.factor(wb$Compound)
  
  return(wb)
}

fix_vars <- function(dat) {
  dat <- dat %>%
    filter(
      !variable %in% c(
        "Object.ID",
        "Object.No.in.Organoids",
        "Object.No.in.Organoids.Selected",
        "Object.No.in.Organoids.Selected.Selected"
      )
    )
  names(dat) <- str_replace(string = names(dat),
                            pattern = "\\.",
                            replacement = " ")
  names(dat)[names(dat) == "Cell Type"] <- "Patient"
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = "\\.",
      replacement = " "
    ) %>%
    str_squish()
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = " 1 px",
      replacement = ""
    )
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = "1 2",
      replacement = "1.2"
    )
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = "2 2",
      replacement = "2.2"
    )
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = "Region.",
      replacement = ""
    )
  
  return(dat)
}

fix_vars2<-function(dat) {
  dat <- dat %>%
    filter(
      !variable %in% c(
        "Object ID",
        "Object No in Organoids",
        "Object No in Organoids Selected",
        "Object No in Organoids Selected Selected"
      )
    )
  
  names(dat)[names(dat) == "Cell Type"] <- "Patient"
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = " 1 px",
      replacement = ""
    ) %>%
    str_squish()
  
  dat$variable <-
    str_replace_all(
      string = dat$variable,
      pattern = "Region.",
      replacement = ""
    )
  
  return(dat)
}


make_map_pdcc_ld <- function(dat, LD_filter = NULL) {
  if (!is.null(LD_filter)) {
    dat <- dat %>% filter(LD == LD_filter)
    
    dat <-
      dat %>% group_by(Compound, Concentration, Patient, Day, LD, variable) %>%
      dplyr::summarise(mean_val = mean(value)) %>%
      spread(variable, mean_val) %>% ungroup() %>%
      unite(Patient,
            Day,
            Compound,
            Concentration,
            col = P_D_C_C,
            remove = FALSE)
    
    drug_colnames <- names(dat)
    mat <-
      dat  %>% column_to_rownames("P_D_C_C") %>% dplyr::select(-Patient,-Day,-Compound,-Concentration,-LD) %>%
      mapply(FUN = scale, center = FALSE) %>%  as.matrix() %>% t()
    colnames(mat) <- dat$P_D_C_C
    
    
    annot <-
      dat  %>%  column_to_rownames("P_D_C_C") %>% dplyr::select(Day, Compound) %>% as.data.frame()
    return(list(mat = mat, annot = annot))
  }
  dat <-
    dat %>% group_by(Compound, Concentration, Patient, Day, LD, variable) %>%
    dplyr::summarise(mean_val = mean(value)) %>%
    spread(variable, mean_val) %>% ungroup() %>%
    unite(Patient,
          Day,
          Compound,
          Concentration,
          LD,
          col = P_D_C_C_LD,
          remove = FALSE)
  
  mat <-
    dat  %>% column_to_rownames("P_D_C_C_LD") %>% dplyr::select(-Patient,-Day,-Compound,-Concentration,-LD) %>%
    as.matrix() %>% t()
  
  annot <-
    dat  %>%  column_to_rownames("P_D_C_C_LD") %>% dplyr::select(Day, Compound, LD) %>% as.data.frame()
  return(list(mat = mat, annot = annot))
  
}

make_map_pdr <- function(dat) {
  dat <- dat  %>% 
    group_by(ExpID, Compound, Patient, Day, variable) %>% 
    dplyr::summarise(mean_val = mean(value)) %>%
    spread(variable, mean_val) %>% ungroup() %>%
    dplyr::select(-Compound) %>% 
    unite(Patient, Day, ExpID, col = P_D_R, remove = FALSE)
  
  mat<-
    dat %>% column_to_rownames("P_D_R") %>% dplyr::select(-Patient,-Day,-ExpID) %>%
    as.matrix() %>% t()
  
  annot <-
    dat %>%  column_to_rownames("P_D_R") %>% dplyr::select(Patient, Day) %>% as.data.frame()
  return(list(mat = mat, annot = annot))
}

#subset_measurements
subset_measurements <-   function(measurements,
                                  feature = 'counts',
                                  missing_as_zero = TRUE) {
  sumfun <-
    paste0('mean(utils::type.convert(`', feature, '`,na.strings=""),na.rm=TRUE)') #The summarization function (mean of the feature)
  vname <-
    paste0('mean_', make.names(feature, allow_ = TRUE)) #The name of the summarized variable
  
  #This is a bit of a helper step: data are imported as characters
  #But concentration is a number. This converts it to the appropriate type
  #TODO: Concentration is a grouping variable, so it needs to be modified early on.
  #There might be a better place for this call.
  measurements<-measurements %>% mutate(Concentration=utils::type.convert(Concentration,na.string="")) 
  
  wd <- tibble()
  
  if (feature == 'counts') {
    #counts are handled differently, as it's not a data column.
    vname = 'counts'
    wd <-
      measurements %>% group_by(Day, LD, Compound, Concentration, `Cell Type`, Row, Column) %>%
      summarise(counts = n())
  } else{
    wd <-
      measurements %>% group_by(Day, LD, Compound, Concentration, `Cell Type`, Row, Column) %>%
      summarise_(.dots = setNames(sumfun, vname))
  }
  
  # This reshapes the table so that for every well the values for dead and live are in columns
  dwd <-
    dcast(wd,
          Day + Compound + Concentration + `Cell Type` + Row + Column ~ LD ,
          value.var = vname)
  
  #TODO: Should cases where Dead=0 set to 0 or leave NA?
  #In some cases, it makes sense to set missing data to 0 (e.g. counts, area). In some others, it would be inocrrect (e.g. roundness)
  if (missing_as_zero) {
    dwd[is.na(dwd$Dead), "Dead"] <- 0
    dwd[is.na(dwd$Live), "Live"] <- 0
  }
  return(dwd)
}

unify<-function(l_wb,r_wb,l_suffix,r_suffix){
  #Merging by UID
  merged_data<-left_join(l_wb,r_wb,by='UID',suffix=c(l_suffix,r_suffix))
  
  duplicated_cols<-names(merged_data)
  l_cols<-grep(paste0(l_suffix,"$"),duplicated_cols)
  r_cols<-grep(paste0(r_suffix,"$"),duplicated_cols)
  unique_cols<-duplicated_cols[-1*c(l_cols,r_cols)]
  duplicated_cols<-duplicated_cols[c(l_cols,r_cols)]
  duplicated_cols<-sub(paste0(l_suffix,'$'),'',duplicated_cols)
  duplicated_cols<-sub(paste0(r_suffix,'$'),'',duplicated_cols)
  duplicated_cols<-unique(duplicated_cols)
  merged_data<-ungroup(merged_data)
  wb_t<-merged_data[,unique_cols]
  tnames<-unique_cols
  
  for(fid in duplicated_cols)  {
    if(all(merged_data[,paste0(fid,r_suffix)] == merged_data[,paste0(fid,l_suffix)],na.rm=TRUE)){
      wb_t<-cbind(wb_t,merged_data[,paste0(fid,r_suffix)])
      tnames<-c(tnames,fid)
    }else{
      wb_t<-cbind(wb_t,merged_data[,paste0(fid,r_suffix)],merged_data[,paste0(fid,l_suffix)])
      tnames<-c(tnames,paste0(fid,r_suffix),paste0(fid,l_suffix))
    }
  }
  names(wb_t)<-tnames
  return(wb_t)  
}

restack<-function(merged_data,foi="LD",l_suffix,r_suffix){
  l_tbl<-merged_data[,-which(names(merged_data)==paste0(foi,r_suffix))]
  l_names<-names(l_tbl)
  names(l_tbl)[which(l_names==paste0(foi,l_suffix))]<-foi
  l_tbl[,"Treatment"]<-paste0(l_tbl[,"Compound"],l_suffix)
  l_tbl[,"Compound"]<-paste0(l_tbl[,"Compound"],l_suffix)
  l_tbl[,"Assay"]<-l_suffix
  
  r_tbl<-merged_data[,-which(names(merged_data)==paste0(foi,l_suffix))]
  r_names<-names(r_tbl)
  names(r_tbl)[which(r_names==paste0(foi,r_suffix))]<-foi
  r_tbl[,"Treatment"]<-paste0(r_tbl[,"Compound"],r_suffix)
  r_tbl[,"Compound"]<-paste0(r_tbl[,"Compound"],r_suffix)
  r_tbl[,"Assay"]<-r_suffix
  return(rbind(l_tbl,r_tbl))
}




plot_curves <-
  function(measurements,
           feature,
           patient = NULL,
           metric = NULL,
           normalized = TRUE,
           error_bars = TRUE,
           freescale = FALSE,
           i_aes=NULL) {
    y_label <- switch(metric,
                      proportion = "Proportion Live/Total",
                      ratio = "Ratio Live/Dead")
    
    measurements <- make_treatments(measurements)
    
    if (normalized) {
      measurements <- measurements %>% rename(val = nlive, dev = nlive_sd)
      y_label <- paste("Normalized", y_label)
    } else{
      measurements <- measurements %>% rename(val = live_mean, dev = live_sd)
    }
    measurements <- measurements %>% filter(abs(val) != Inf)
    pd <- position_identity()
    days <- sort(as.numeric(unique(measurements$Day)))
    
    if(is.null(i_aes)){
      cols <-
        make_aes(measurements)$cols
      shapes <-
        make_aes(measurements)$shapes
    } else{
      cols<-i_aes$cols
      shapes<-i_aes$shapes
    }
    
    pl <- ggplot(data = measurements) + aes(
      x = Day,
      y = val,
      col = Treatment,
      shape = Compound
    )
    
    if (error_bars) {
      pd <- position_dodge(0.1)
      pl <-
        pl + geom_errorbar(aes(ymin = val - dev, ymax = val + dev),
                           width = 1,
                           position = pd)
    }
    
    pl <- pl  +
      geom_line(lwd = 1) +
      geom_point(position = pd) +
      scale_x_continuous(breaks = days) +
      scale_color_manual(values = cols) +
      ylab(label = y_label) +
      labs(title = ifelse(is.null(patient), '', paste("Patient", patient)), subtitle =
             feature) +
      scale_shape_manual(values = unique(shapes), guide=FALSE) +
      guides(colour = guide_legend(override.aes = list(shape = shapes ))) +
      theme_linedraw()
    
    if (freescale) {
      pl <- pl +
        facet_grid(`Cell Type` ~ ., scales = "free_y")
    } else{
      if (error_bars) {
        maxy <- max(measurements$val + measurements$dev, na.rm = TRUE)
        miny <-
          ifelse(min(measurements$val - measurements$dev, na.rm = TRUE) > 0,
                 0,
                 min(measurements$val - measurements$dev, na.rm = TRUE))
      } else{
        maxy <- max(measurements$val, na.rm = TRUE)
        miny <-
          ifelse(min(measurements$val, na.rm = TRUE) > 0,
                 0,
                 min(measurements$val, na.rm = TRUE))
      }
      if(miny<0){
        miny<-0
      }else {
        miny <- floor(miny)
      }
      maxy <- ifelse(maxy < 2, round(maxy, 1), ceiling(maxy))
      scale_fac <- ifelse(maxy < 2, .2, round(maxy / 5, 1))
      pl <- pl +
        facet_grid(`Cell Type` ~ .) +
        scale_y_continuous(breaks = seq(miny, maxy, by = scale_fac)) +
        expand_limits(y = miny)
    }
    
    
  }

plot_c_curves <-
  function(measurements,
           feature,
           patient = NULL,
           metric = NULL,
           normalized = TRUE,
           error_bars = TRUE,
           freescale = FALSE,
           i_aes=NULL) {
    y_label <- switch(metric,
                      proportion = "Proportion Live/Total",
                      ratio = "Ratio Live/Dead")
    
    measurements <- make_treatments(measurements)
    
    if (normalized) {
      measurements <- measurements %>% rename(val = nlive, dev = nlive_sd)
      y_label <- paste("Normalized", y_label)
    } else{
      measurements <- measurements %>% rename(val = live_mean, dev = live_sd)
    }
    measurements <- measurements %>% filter(abs(val) != Inf)
    pd <- position_identity()
    days <- sort(as.numeric(unique(measurements$Day)))
    n_days <-length(days)
    
    if(is.null(i_aes)){
      cols <-
        make_aes(measurements)$cols
      shapes <-
        make_aes(measurements)$shapes
    } else{
      cols<-i_aes$cols
      shapes<-i_aes$shapes
    }
    
    cols=c('black','red')
    
    pl <- ggplot(data = measurements) + aes(
      x = Day,
      y = val,
      col = Compound,
      #      shape = Mode,
      linetype = Assay
    )
    
    if (error_bars) {
      pd <- position_dodge(0.1)
      pl <-
        pl + geom_errorbar(aes(ymin = val - dev, ymax = val + dev),
                           width = 1,
                           position = pd)
    }
    
    
    pl <- pl  +
      #      geom_line(size = 1, linetype=c(rep(1,n_days),rep(0,n_days),rep(2,n_days),rep(3,n_days) )) +
      geom_path(size = 1) +
      geom_point(position = pd) +
      scale_x_continuous(breaks = days) +
      scale_color_manual(values = cols) +
      ylab(label = y_label) +
      labs(title = ifelse(is.null(patient), '', paste("Patient", patient)), subtitle =
             feature) +
      #     scale_shape_manual(values = unique(shapes), guide=FALSE) +
      scale_linetype_manual(values = c(1,3,1,3)) +
      #      guides(colour = guide_legend(override.aes = list(shape = shapes ))) +
      
      theme_linedraw()
    
    if (freescale) {
      pl <- pl +
        facet_grid(`Cell Type` ~ ., scales = "free_y")
    } else{
      if (error_bars) {
        maxy <- max(measurements$val + measurements$dev, na.rm = TRUE)
        miny <-
          ifelse(min(measurements$val - measurements$dev, na.rm = TRUE) > 0,
                 0,
                 min(measurements$val - measurements$dev, na.rm = TRUE))
      } else{
        maxy <- max(measurements$val, na.rm = TRUE)
        miny <-
          ifelse(min(measurements$val, na.rm = TRUE) > 0,
                 0,
                 min(measurements$val, na.rm = TRUE))
      }
      miny <- floor(miny)
      maxy <- ifelse(maxy < 2, round(maxy, 1), ceiling(maxy))
      scale_fac <- ifelse(maxy < 2, .2, round(maxy / 5, 1))
      pl <- pl +
        facet_grid(`Cell Type` ~ .) +
        scale_y_continuous(breaks = seq(miny, maxy, by = scale_fac)) +
        expand_limits(y = miny)
    }
    
    
  }

normalize_data <- function(dwd, metric = 'proportion') {
  dwd$plive <- NA
  if (metric != 'proportion') {
    dwd$plive <- dwd$Live / dwd$Dead
  } else{
    dwd$plive <- dwd$Live / (dwd$Live + dwd$Dead)
    
  }
  
  #Normalizes the proportion/ratios by the corresponding value at day 0
  ndwd <- dwd %>%
    group_by(`Cell Type`, Compound, Concentration, Day) %>%
    summarise(
      N_wells = sum(!is.na(plive)),
      live_mean = mean(plive, na.rm = TRUE),
      live_sd = sd(plive, na.rm = TRUE)
    ) %>%
    mutate(
      nlive = live_mean / first(live_mean, order_by = 'Day'),
      nlive_sd = live_sd / first(live_mean, order_by = 'Day')
    )
  
  return(ndwd)
}


make_aes <- function(measurements) {
  treatments <-
    measurements %>% ungroup() %>% select(Compound, Concentration)  %>% distinct()
  nmedia <-
    unlist(treatments %>% filter(tolower(Compound) == "media" | tolower(Compound) == "control" ) %>% count())
  if (nmedia  == 1) {
    treatments <- treatments %>% filter(!(tolower(Compound) == "media" | tolower(Compound) == "control"))
  }
  
  ncats <-
    as.integer(treatments %>% select(Compound) %>% distinct() %>% count())
  dark_cols <- few_pal(palette = 'Dark')(ncats)
  light_cols <- few_pal(palette = 'Light')(ncats)
  
  cat_tbl <-
    data.frame(
      treatments %>% group_by(Compound) %>% summarise(n = n()),
      from = light_cols,
      to = dark_cols
    )
  cat_tbl$shape <-
    as.numeric(row.names(cat_tbl)) - 1 #if there are >20 categories, it gets messy
  
  cols <-
    unlist(apply(cat_tbl, 1, function(x)
      colorRampPalette(c(
        ifelse(as.numeric(x[2]) == 1, x[4], x[3]), x[4]
      ))(as.numeric(x[2]))))
  
  shapes <-
    unlist(apply(cat_tbl, 1, function(x)
      rep(as.numeric(x[5]), as.numeric(x[2]))))
  
  if (nmedia == 1) {
    cols <- c('#000000', cols)
    shapes <- c(32, shapes)
  }
  return(list(cols = cols, shapes = shapes))
}


make_treatments <- function(measurements) {
  
  #    if(all(is.na(measurements$Concentration))){
  #    measurements<-get_concentration(measurements)
  #  }
  
  treatments <-
    measurements %>% ungroup() %>% select(Compound, Concentration)  %>% distinct()
  
  
  #Separate controls and response curves. Controls alpha sorted come first,
  #then response cureves sorted alpha on compond, then increasing concentration
  k<-unlist(treatments %>% group_by(Compound) %>%summarise(n=n()) %>% filter(n==1) %>%select(Compound))
  concs <-
    treatments %>% filter(! Compound %in% k) %>% arrange(Compound, Concentration) #Compounds with multiple concs
  ctrls <-
    treatments %>% filter(Compound %in% k) %>% arrange(Compound)#Controls.
  
  #Move "Media" to the top of the list
  media <- ctrls %>% filter(tolower(Compound) == "media" | tolower(Compound) == "control")
  if (media %>% count() == 1) {
    ctrls <-
      ctrls %>% filter(!(tolower(Compound) == "media" | tolower(Compound) == "control")) %>% rbind(media, .)
  }
  
  #Merge the controls and curves
  treatments <- rbind(ctrls, concs)
  
  #Treatment is combination of compound and concentration
  treatments <-
    treatments %>% mutate(Treatment = paste0(Compound, ifelse(
      is.na(Concentration), '', paste0("_", Concentration)
    )))
  #TODO: Concentration should be rounded before paste
  
  measurements <-
    measurements %>% mutate(Treatment = paste0(Compound, ifelse(
      is.na(Concentration), '', paste0("_", Concentration)
    )))
  
  measurements$Treatment <-
    ordered(measurements$Treatment, levels = treatments$Treatment)
  measurements$Compound <-
    ordered(measurements$Compound, levels = unlist(treatments %>% select(Compound) %>% distinct()))
  
  return(measurements)
}

make_sankey<-function(dat){
  dat<- dat %>% mutate(group = paste(Timepoint,surv_status, sep = '_')) %>% 
     group_by(unique_ID)
  
  organoids_sankey_0_1_4_1 <- dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_L-VD" & nth(group,2) == "1_D-ML_D-VD" ~ "4_1")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "4_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_4_2 <- dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_L-VD" & nth(group,2) == "1_D-ML_L-VD" ~ "4_2")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "4_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_4_3 <- dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_L-VD" & nth(group,2) == "1_L-ML_D-VD" ~ "4_3")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "4_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_4_4 <- dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_L-VD" & nth(group,2) == "1_L-ML_L-VD" ~ "4_4")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "4_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_3_1 <- dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_D-VD" & nth(group,2) == "1_D-ML_D-VD" ~ "3_1")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "3_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_3_2 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_D-VD" & nth(group,2) == "1_D-ML_L-VD" ~ "3_2")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "3_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_3_3 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_D-VD" & nth(group,2) == "1_L-ML_D-VD" ~ "3_3")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "3_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_3_4 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_L-ML_D-VD" & nth(group,2) == "1_L-ML_L-VD" ~ "3_4")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "3_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_2_1 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_L-VD" & nth(group,2) == "1_D-ML_D-VD" ~ "2_1")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "2_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_2_2 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_L-VD" & nth(group,2) == "1_D-ML_L-VD" ~ "2_2")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "2_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_2_3 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_L-VD" & nth(group,2) == "1_L-ML_D-VD" ~ "2_3")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "2_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_0_1_2_4 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_L-VD" & nth(group,2) == "1_L-ML_L-VD" ~ "2_4")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "2_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_0_1_1_1 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_D-VD" & nth(group,2) == "1_D-ML_D-VD" ~ "1_1")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "1_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_1_2 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_D-VD" & nth(group,2) == "1_D-ML_L-VD" ~ "1_2")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "1_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_1_3 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_D-VD" & nth(group,2) == "1_L-ML_D-VD" ~ "1_3")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "1_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_0_1_1_4 <-  dat %>% 
    mutate(source_target = case_when(first(group) == "0_D-ML_D-VD" & nth(group,2) == "1_L-ML_L-VD" ~ "1_4")) %>% 
    arrange(unique_ID)  %>% do(head(., 2)) %>% filter(source_target == "1_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  #### 1-3 ####
  
  organoids_sankey_1_3_4_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_L-VD" & nth(group,3) == "3_D-ML_D-VD" ~ "4_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "4_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_4_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_L-VD" & nth(group,3) == "3_D-ML_L-VD" ~ "4_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "4_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_4_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_L-VD" & nth(group,3) == "3_L-ML_D-VD" ~ "4_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "4_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_4_4 <- dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_L-VD" & nth(group,3) == "3_L-ML_L-VD" ~ "4_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "4_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_3_1 <- dat %>% 
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_D-VD" & nth(group,3) == "3_D-ML_D-VD" ~ "3_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "3_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_3_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_D-VD" & nth(group,3) == "3_D-ML_L-VD" ~ "3_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "3_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_3_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_D-VD" & nth(group,3) == "3_L-ML_D-VD" ~ "3_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "3_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_3_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_L-ML_D-VD" & nth(group,3) == "3_L-ML_L-VD" ~ "3_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "3_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_2_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_L-VD" & nth(group,3) == "3_D-ML_D-VD" ~ "2_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "2_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_2_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_L-VD" & nth(group,3) == "3_D-ML_L-VD" ~ "2_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "2_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_2_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_L-VD" & nth(group,3) == "3_L-ML_D-VD" ~ "2_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "2_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_1_3_2_4 <- dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_L-VD" & nth(group,3) == "3_L-ML_L-VD" ~ "2_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "2_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_1_3_1_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_D-VD" & nth(group,3) == "3_D-ML_D-VD" ~ "1_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "1_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_1_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_D-VD" & nth(group,3) == "3_D-ML_L-VD" ~ "1_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "1_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_1_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_D-VD" & nth(group,3) == "3_L-ML_D-VD" ~ "1_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "1_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_1_3_1_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,2) == "1_D-ML_D-VD" & nth(group,3) == "3_L-ML_L-VD" ~ "1_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 2:3)) %>% filter(source_target == "1_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  #### 3-7 ####
  
  organoids_sankey_3_7_4_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_L-VD" & nth(group,4) == "7_D-ML_D-VD" ~ "4_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "4_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_4_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_L-VD" & nth(group,4) == "7_D-ML_L-VD" ~ "4_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "4_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_4_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_L-VD" & nth(group,4) == "7_L-ML_D-VD" ~ "4_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "4_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_4_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_L-VD" & nth(group,4) == "7_L-ML_L-VD" ~ "4_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "4_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_3_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_D-VD" & nth(group,4) == "7_D-ML_D-VD" ~ "3_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "3_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_3_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_D-VD" & nth(group,4) == "7_D-ML_L-VD" ~ "3_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "3_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_3_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_D-VD" & nth(group,4) == "7_L-ML_D-VD" ~ "3_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "3_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_3_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_L-ML_D-VD" & nth(group,4) == "7_L-ML_L-VD" ~ "3_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "3_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_2_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_L-VD" & nth(group,4) == "7_D-ML_D-VD" ~ "2_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "2_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_2_2 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_L-VD" & nth(group,4) == "7_D-ML_L-VD" ~ "2_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "2_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_2_3 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_L-VD" & nth(group,4) == "7_L-ML_D-VD" ~ "2_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "2_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_3_7_2_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_L-VD" & nth(group,4) == "7_L-ML_L-VD" ~ "2_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "2_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  organoids_sankey_3_7_1_1 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_D-VD" & nth(group,4) == "7_D-ML_D-VD" ~ "1_1")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "1_1") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_1_2 <- dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_D-VD" & nth(group,4) == "7_D-ML_L-VD" ~ "1_2")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "1_2") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_1_3 <- dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_D-VD" & nth(group,4) == "7_L-ML_D-VD" ~ "1_3")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "1_3") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  organoids_sankey_3_7_1_4 <-  dat %>%
    mutate(source_target = case_when(nth(group,3) == "3_D-ML_D-VD" & nth(group,4) == "7_L-ML_L-VD" ~ "1_4")) %>% 
    arrange(unique_ID)  %>% do(slice(., 3:4)) %>% filter(source_target == "1_4") %>% 
    group_by(group) %>% dplyr::summarise(count = n())
  
  
  
  
  
  bound_sankey <- bind_rows(organoids_sankey_0_1_4_1, organoids_sankey_0_1_4_2, organoids_sankey_0_1_4_3, organoids_sankey_0_1_4_4, 
                            organoids_sankey_0_1_3_1, organoids_sankey_0_1_3_2, organoids_sankey_0_1_3_3, organoids_sankey_0_1_3_4, 
                            organoids_sankey_0_1_2_1, organoids_sankey_0_1_2_2, organoids_sankey_0_1_2_3, organoids_sankey_0_1_2_4, 
                            organoids_sankey_0_1_1_1, organoids_sankey_0_1_1_2, organoids_sankey_0_1_1_3, organoids_sankey_0_1_1_4,
                            organoids_sankey_1_3_4_1, organoids_sankey_1_3_4_2, organoids_sankey_1_3_4_3, organoids_sankey_1_3_4_4, 
                            organoids_sankey_1_3_3_1, organoids_sankey_1_3_3_2, organoids_sankey_1_3_3_3, organoids_sankey_1_3_3_4, 
                            organoids_sankey_1_3_2_1, organoids_sankey_1_3_2_2, organoids_sankey_1_3_2_3, organoids_sankey_1_3_2_4, 
                            organoids_sankey_1_3_1_1, organoids_sankey_1_3_1_2, organoids_sankey_1_3_1_3, organoids_sankey_1_3_1_4,
                            organoids_sankey_3_7_4_1, organoids_sankey_3_7_4_2, organoids_sankey_3_7_4_3, organoids_sankey_3_7_4_4, 
                            organoids_sankey_3_7_3_1, organoids_sankey_3_7_3_2, organoids_sankey_3_7_3_3, organoids_sankey_3_7_3_4, 
                            organoids_sankey_3_7_2_1, organoids_sankey_3_7_2_2, organoids_sankey_3_7_2_3, organoids_sankey_3_7_2_4, 
                            organoids_sankey_3_7_1_1, organoids_sankey_3_7_1_2, organoids_sankey_3_7_1_3, organoids_sankey_3_7_1_4)
  
  source <- bound_sankey[seq(from = 1, to = length(bound_sankey$group), by = 2),] %>% rename(source = group)
  target <- bound_sankey[seq(from = 2, to = length(bound_sankey$group), by = 2),] %>% dplyr::select(target = group)
  links <- bind_cols(source,target)
  
  k<-as.data.frame(source$source %>% str_match("_(.)-ML_(.)"))
  k$group<-"type_a"
  k[k[,2]==k[,3],4]<-"type_b"
  
  links$group<-k$group
    
  
  # STA concoordance
  # links$group = as.factor(c("type_b","type_b","type_b","type_b","type_a","type_a","type_a","type_a","type_a","type_b",
  #                           "type_b","type_b","type_b","type_b","type_b","type_a","type_a","type_a","type_a","type_a",
  #                           "type_a","type_b","type_b","type_b","type_b","type_b","type_a","type_a","type_a","type_a",
  #                           "type_b"))
  
  nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
  
  links$IDsource=match(links$source, nodes$name)-1 
  links$IDtarget=match(links$target, nodes$name)-1
  
  
  ### group node coloring for concoordinate vs discoordinate STA
  nodes$group=as.factor(c("a","b","b","a","a","b","b","a","a","b","b","a","a","b","b","a"))
  return(list(links=links,nodes=nodes))
}

theme_Publication <-
  function(base_size = 14,
           base_family = "Helvetica") {
    (
      theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2),
          hjust = 0.5
        ),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 1),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.4, "cm"),
        # for spacing
        legend.spacing = unit(0, "cm"),
        legend.title = element_text(face = "italic"),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        strip.text = element_text(face = "bold")
      )
    )
  }

scale_fill_Publication <- function(...) {
  discrete_scale("fill", "Publication", manual_pal(
    values = c(
      "#386cb0",
      "#fdb462",
      "#7fc97f",
      "#ef3b2c",
      "#662506",
      "#a6cee3",
      "#fb9a99",
      "#984ea3",
      "#ffff33"
    )
  ), ...)
  
}

scale_color_Publication <- function(...) {
  discrete_scale("colour", "Publication", manual_pal(
    values = c(
      "#386cb0",
      "#fdb462",
      "#7fc97f",
      "#ef3b2c",
      "#662506",
      "#a6cee3",
      "#fb9a99",
      "#984ea3",
      "#ffff33"
    )
  ), ...)
  
}
## ggplot with Publication

ggp <-
  function(...)
    ggplot(...) + theme_Publication() + scale_fill_Publication() + scale_color_Publication()
ggp_bare <- function(...)
  ggplot(...) + theme_Publication()
