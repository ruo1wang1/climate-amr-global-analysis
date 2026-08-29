###############################################################################
# Analysis of climate-AMR associations under alternative reanalysis-product
# configurations
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse); library(mgcv); library(ggplot2); library(patchwork)
  library(grid); library(scales); library(zoo); library(gridExtra)
  library(cowplot); library(viridis); library(splines)
})

# Pragmatic thresholds for robustness classification
reversal_min_abs_dlogor <- 0.05
reversal_min_secondary_abs_dlogor <- 0.03
robust_curve_r_threshold <- 0.80
minor_curve_r_threshold <- 0.60
robust_shape_agreement_threshold <- 0.50
minor_shape_agreement_threshold <- 0.33
robust_mag_ratio_lower <- 0.50
robust_mag_ratio_upper <- 2.00
minor_mag_ratio_lower <- 0.25
minor_mag_ratio_upper <- 3.00

# ── Paths ──
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root,"outputs/historical_associations","model_ready_inputs")
lag_summary_path <- file.path(revision_root,
  "data",
  "source_data",
  "lag_selection",
  "historical_model_c",
  "01_csv",
  "historical_lag_summary_model_c.csv")
output_root <- file.path(revision_root,"outputs","reanalysis_product_configurations")

model_ready_input_name <- function(code) {
  paste0(code, "_model_ready_data.csv")
}

dirs <- list(
  root=output_root,
  exposure =file.path(output_root,"01_exposure_comparison"),
  model    =file.path(output_root,"02_model_performance"),
  curve    =file.path(output_root,"03_curve_comparison"),
  robust   =file.path(output_root,"04_robustness_assessment"),
  summaries=file.path(output_root,"05_model_summaries"),
  bridge   =file.path(output_root,"06_bridge_table"))
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive=TRUE)

# ── Bacteria ──
bact_specs <- list(
  list(code="3GCREC",title="3GCR-Ec",fn=model_ready_input_name("3GCREC")),
  list(code="3GCRKP",title="3GCR-Kp",fn=model_ready_input_name("3GCRKP")),
  list(code="CRAB",  title="CR-Ab",  fn=model_ready_input_name("CRAB")),
  list(code="CREC",  title="CR-Ec",  fn=model_ready_input_name("CREC")),
  list(code="CRKP",  title="CR-Kp",  fn=model_ready_input_name("CRKP")),
  list(code="CRPA",  title="CR-Pa",  fn=model_ready_input_name("CRPA")))

# ── Climate source column mappings ──
clim_cols <- list(
  Primary=list(TMP="TMP",PREC="PREC",WET="WET",HUM="HUM"),
  ERA5   =list(TMP="ERA5_annual_mean_t2m_c",PREC="ERA5_annual_precip_mm",
               WET="ERA5_wet_days",HUM="HUM"),
  MERRA2 =list(TMP="MERRA2_annual_mean_t2m_c",PREC="MERRA2_annual_precip_mm",
               WET="MERRA2_wet_days",HUM="MERRA2_annual_mean_rh_pct"))

# ── Colorblind-safe palette with strong separation ──
src_labels <- c(Primary="Primary (CRU TS + ERA5 RH)", ERA5="ERA5", MERRA2="MERRA-2")
src_colors <- c(Primary="#B24745", ERA5="#2F6DB2", MERRA2="#138A72")
src_ltys   <- c(Primary="solid",   ERA5="dashed",  MERRA2="dotdash")
src_shape_values <- c(Primary=21, ERA5=22, MERRA2=24)
src_ribbon_alpha <- 0.055
legend_ribbon_alpha <- 0.17
pair_ribbon_alpha <- 0.11
rmse_fill_palette <- c("#FFF7EC", "#FDD49E", "#F67E4B", "#C4411C")
bias_div_palette <- c("#3C6EBA", "#F7F6F3", "#B24745")
smooth_fill_palette <- c("#F7FBFF", "#DCEAF6", "#9ECAE1", "#4292C6", "#084594")

bacteria_order <- c("3GCR-Ec","3GCR-Kp","CR-Ab","CR-Ec","CR-Kp","CR-Pa")
curve_variable_order <- c("Temperature (\u00B0C)", "Relative Humidity (%)", "Precipitation (mm)", "Wet Days (d)")
pair_levels <- c("Primary vs ERA5", "Primary vs MERRA-2", "ERA5 vs MERRA-2")
pair_colors <- c("Primary vs ERA5"="#7E5AA6", "Primary vs MERRA-2"="#2A9D8F", "ERA5 vs MERRA-2"="#4C78A8")

clim_labels <- c(TMP="Temperature (\u00B0C)",PREC="Precipitation (mm)",
                 WET="Wet Days (d)",HUM="Relative Humidity (%)")
cv2label    <- c(Temperature="Temperature (\u00B0C)",Precipitation="Precipitation (mm)",
                 WetDays="Wet Days (d)",Humidity="Relative Humidity (%)")
vn2label    <- c(TMP_scaled_lag="Temperature (\u00B0C)",HUM_scaled_lag="Relative Humidity (%)",
                 PREC_scaled_lag="Precipitation (mm)",WET_scaled_lag="Wet Days (d)")

is_rh_same <- function(v,s1,s2) grepl("HUM|Humidity",v,ignore.case=TRUE) && setequal(c(s1,s2),c("Primary","ERA5"))

pair_label <- function(s1, s2) {
  if (setequal(c(s1, s2), c("Primary", "ERA5"))) return("Primary vs ERA5")
  if (setequal(c(s1, s2), c("Primary", "MERRA2"))) return("Primary vs MERRA-2")
  if (setequal(c(s1, s2), c("ERA5", "MERRA2"))) return("ERA5 vs MERRA-2")
  paste(sort(c(s1, s2)), collapse = " vs ")
}

theme_s <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill="#F0F0F0",colour="grey70"),
  strip.text=element_text(face="bold",size=10),
  plot.title=element_text(hjust=0.5,face="bold",size=12),
  plot.subtitle=element_text(hjust=0.5,size=9,color="grey40"),
  legend.position="bottom",legend.box="horizontal")

# ── Lag settings ──
lag_settings <- {
  df<-read.csv(lag_summary_path,stringsAsFactors=FALSE)
  out<-setNames(vector("list",nrow(df)),df$Display_Name)
  for(i in seq_len(nrow(df))) out[[df$Display_Name[i]]]<-list(
    temp_lag=as.integer(df$TMP_lag[i]),precip_lag=as.integer(df$PREC_lag[i]),
    humid_lag=as.integer(df$HUM_lag[i]),wetdays_lag=as.integer(df$WET_lag[i]))
  out}

clim_vars <- c("TMP_scaled_lag","PREC_scaled_lag","HUM_scaled_lag","WET_scaled_lag")

# Error log
error_log <- list()

# All helper functions

get_pls <- function(d) {
  cc<-paste0("PLS_Comp",1:4); p<-cc[cc%in%names(d)]
  p[sapply(p,function(x)!all(is.na(d[[x]])))]
}

# ── prepare_data: IDENTICAL to main Model C ──
prepare_data <- function(file_path, bname, source="Primary", sample_mode="full") {
  lc <- lag_settings[[bname]]
  if(is.null(lc)) lc<-list(temp_lag=3,precip_lag=3,humid_lag=1,wetdays_lag=1)
  data <- read.csv(file_path)

  # ═══ Column swap ═══
  if (source != "Primary") {
    sc<-clim_cols[[source]]
    data$TMP<-data[[sc$TMP]]; data$PREC<-data[[sc$PREC]]
    data$WET<-data[[sc$WET]]; data$HUM<-data[[sc$HUM]]
  }

  # ═══ Common-sample filter (data-driven) ═══
  if (sample_mode == "common") {
    check_cols <- c("TMP","PREC","WET","HUM",
      "ERA5_annual_mean_t2m_c","ERA5_annual_precip_mm","ERA5_wet_days",
      "MERRA2_annual_mean_t2m_c","MERRA2_annual_precip_mm",
      "MERRA2_annual_mean_rh_pct","MERRA2_wet_days")
    avail <- intersect(check_cols, names(data))
    valid_row <- apply(data[,avail,drop=FALSE],1,function(r){
      v<-suppressWarnings(as.numeric(r)); all(!is.na(v)&is.finite(v))})
    if ("CRUTS_country_match_type" %in% names(data))
      valid_row <- valid_row & (data$CRUTS_country_match_type != "unresolved")
    data <- data[valid_row, ]
  }

  # ═══ From here: IDENTICAL to main Model C ═══
  data <- data %>%
    mutate(year=as.numeric(as.character(year)), Region=factor(Region),
           climate_zone=factor(case_when(abs(lat)>66.5~"Polar Zone",
                                         abs(lat)>23.5~"Temperate Zone",TRUE~"Tropical Zone"))) %>%
    group_by(NAME) %>% mutate(location_id=cur_group_id()) %>% ungroup()
  data <- data %>% mutate(HUM=pmin(HUM,100))

  sp <- data %>% summarise(across(c(TMP,PREC,HUM,WET),list(mean=~mean(.,na.rm=T),sd=~sd(.,na.rm=T))))

  data <- data %>%
    mutate(TMP_orig=TMP,PREC_orig=PREC,HUM_orig=HUM,WET_orig=WET) %>%
    group_by(climate_zone) %>%
    mutate(across(c(TMP,PREC,HUM,WET),\(x)as.vector(scale(x)),.names="{.col}_scaled")) %>%
    group_by(location_id) %>% arrange(year) %>%
    mutate(TMP_scaled_lag=lag(TMP_scaled,lc$temp_lag),PREC_scaled_lag=lag(PREC_scaled,lc$precip_lag),
           HUM_scaled_lag=lag(HUM_scaled,lc$humid_lag),WET_scaled_lag=lag(WET_scaled,lc$wetdays_lag)) %>%
    filter(!is.na(TMP_scaled_lag)&!is.na(PREC_scaled_lag)&!is.na(HUM_scaled_lag)&!is.na(WET_scaled_lag)) %>%
    ungroup()

  list(data=data, scale_params=sp, n_obs=nrow(data), n_countries=n_distinct(data$NAME))
}

# ── build_gamm: IDENTICAL to main Model C ──
build_gamm <- function(data, bname) {
  ctrl<-gam.control(nthreads=4,maxit=1000,mgcv.tol=1e-7,mgcv.half=15)
  pls<-get_pls(data); if(length(pls)==0)stop("No PLS")
  pt<-paste0("s(",pls,",k=10,bs='cr')",collapse=" + ")
  fs<-paste0("logit_R~s(TMP_scaled_lag,k=5,bs='cr')+s(PREC_scaled_lag,k=10,bs='cr')+",
    "s(HUM_scaled_lag,k=10,bs='cr')+s(WET_scaled_lag,k=10,bs='cr')+",pt,
    "+s(lat,lon,bs='sos',k=20)+s(year,bs='cr',k=8)+s(Region,bs='re')+climate_zone")
  tryCatch(bam(as.formula(fs),data=data,family=gaussian(),method="REML",select=TRUE,control=ctrl),
    error=function(e){warning("bam->gam: ",e$message)
      gam(as.formula(fs),data=data,family=gaussian(),method="REML",select=TRUE)})
}

# ── save_model_summary: writes detailed txt (FIX for empty 05_model_summaries) ──
save_model_summary <- function(model, bt, src, mode, prep, out_dir) {
  fn <- file.path(out_dir, paste0("Sensitivity_",gsub("[^A-Za-z0-9]+","_",bt),
                                   "_",src,"_",mode,"_summary.txt"))
  ls <- lag_settings[[bt]]
  txt <- c(
    strrep("=",60),
    paste("Reanalysis Sensitivity Model Summary:", bt),
    strrep("=",60),
    paste("Climate source:", src),
    paste("Sample mode:", mode),
    paste("Sample size:", prep$n_obs),
    paste("Countries:", prep$n_countries),
    paste0("Lags: TMP=",ls$temp_lag,", PREC=",ls$precip_lag,
           ", HUM=",ls$humid_lag,", WET=",ls$wetdays_lag),
    "", "Model summary:",
    capture.output(print(summary(model))),
    "", "Model diagnostics:",
    capture.output(tryCatch(gam.check(model), error=function(e) cat("gam.check error\n")))
  )
  writeLines(txt, fn)
}

# ── extract_curve: IDENTICAL to main Model C ──
extract_curve <- function(model, data, sp, vn, bname) {
  ov<-str_extract(vn,"^[A-Z]+")
  vm<-sp[[paste0(ov,"_mean")]]; vs<-sp[[paste0(ov,"_sd")]]
  is_p<-grepl("PREC",vn);is_h<-grepl("HUM",vn);is_t<-grepl("TMP",vn);is_w<-grepl("WET",vn)
  vmin<-if(is_p)0 else if(is_h)30 else if(is_w)0 else -10
  vmax<-if(is_h)100 else if(is_p)3200 else if(is_w)300 else 40
  px<-seq((vmin-vm)/vs,(vmax-vm)/vs,length.out=1200); ox<-px*vs+vm
  if(is_p||is_w){ox[ox<=0]<-0.01;if(is_p)ox<-pmin(ox,3200);if(is_w)ox<-pmin(ox,300)}
  if(is_h) ox<-pmin(ox,100)
  pd<-data.frame(TMP_scaled_lag=mean(data$TMP_scaled_lag,na.rm=T),
    PREC_scaled_lag=mean(data$PREC_scaled_lag,na.rm=T),
    HUM_scaled_lag=mean(data$HUM_scaled_lag,na.rm=T),
    WET_scaled_lag=mean(data$WET_scaled_lag,na.rm=T),
    lat=mean(data$lat,na.rm=T),lon=mean(data$lon,na.rm=T),year=mean(data$year,na.rm=T),
    climate_zone=factor(names(which.max(table(data$climate_zone))),levels=levels(data$climate_zone)),
    Region=factor(names(which.max(table(data$Region))),levels=levels(data$Region)))
  for(c in get_pls(data)) pd[[c]]<-mean(data[[c]],na.rm=T)
  pd<-pd[rep(1,1200),]; pd[[vn]]<-px
  suppressWarnings({pred<-predict(model,pd,type="terms",se.fit=TRUE)})
  vc<-grep(vn,colnames(pred$fit)); fv<-pred$fit[,vc]; sv<-pred$se.fit[,vc]
  orv<-exp(fv);ol95<-exp(fv-1.96*sv);ou95<-exp(fv+1.96*sv)
  lf<-function(y,sp){x<-seq_along(y);predict(loess(y~x,span=sp))}
  d2<-diff(diff(orv));dvar<-var(d2,na.rm=T)
  span<-min(0.5,max(0.1,0.3-0.2*sqrt(dvar)))
  if(is_w){span<-min(span,0.22);if(bname=="CR-Kp")span<-min(span,0.18)
  }else if(bname=="CR-Kp"&&is_p){span<-min(span,0.2)
  }else if(bname=="3GCR-Kp"&&is_h){span<-min(span,0.2)
  }else if(bname=="CR-Ab"&&is_t){span<-min(span,0.18)
  }else if((bname=="3GCR-Kp"||bname=="CR-Ab")&&is_p){span<-min(span,0.18)
  }else if(bname=="CR-Pa"&&is_p){span<-min(span,0.19)
  }else if(bname=="CR-Pa"&&is_t){span<-min(span,0.15)}
  tibble(x_orig=ox,x_scaled=px,or=lf(orv,span),or_lower95=lf(ol95,span),
         or_upper95=lf(ou95,span),or_raw=orv,log_or=fv)
}

get_smooth <- function(model) {
  sm<-summary(model);st<-as.data.frame(sm$s.table);st$term<-rownames(st);rownames(st)<-NULL
  st%>%filter(grepl("TMP_scaled_lag|PREC_scaled_lag|HUM_scaled_lag|WET_scaled_lag",term))%>%
    mutate(CV=case_when(grepl("TMP",term)~"Temperature",grepl("PREC",term)~"Precipitation",
      grepl("HUM",term)~"Humidity",grepl("WET",term)~"WetDays"),
      EDF=edf,F_stat=F,Pval=`p-value`,Sig=Pval<0.05,
      SigL=case_when(Pval<0.001~"***",Pval<0.01~"**",Pval<0.05~"*",Pval<0.1~".",TRUE~"ns"),
      PenOut=EDF<0.01)%>%select(CV,EDF,F_stat,Pval,Sig,SigL,PenOut)
}

pct_contrast <- function(cv, data, vn) {
  ov<-str_extract(vn,"^[A-Z]+"); oc<-paste0(ov,"_orig")
  if(!oc%in%names(data))return(NULL)
  ob<-data[[oc]]; ob<-ob[!is.na(ob)&is.finite(ob)]; if(length(ob)<20)return(NULL)
  p10<-quantile(ob,.1);p50<-quantile(ob,.5);p90<-quantile(ob,.9)
  g<-function(xv){i<-which.min(abs(cv$x_orig-xv));cv$or[i]}
  o10<-g(p10);o50<-g(p50);o90<-g(p90)
  tibble(P10=as.numeric(p10),P50=as.numeric(p50),P90=as.numeric(p90),
    OR_P10=o10,OR_P50=o50,OR_P90=o90,ORratio=o90/o10,DlogOR=log(o90)-log(o10),Dir=sign(log(o90)-log(o10)))
}

classify_shape <- function(x,y) {
  if(length(x)<10||all(is.na(y))) return("?")
  lo<-tryCatch(loess(y~x,span=0.3),error=function(e)NULL); if(is.null(lo))return("?")
  ys<-predict(lo); if(any(is.na(ys)))ys[is.na(ys)]<-approx(x,y,xout=x[is.na(ys)])$y
  if(all(is.na(ys)))return("?"); dy<-diff(ys);sc<-sum(diff(sign(dy))!=0,na.rm=T)
  ot<-cor(x,ys,use="complete.obs"); if(is.na(ot))ot<-0
  if(sc==0){if(ot>0.3)"Mono+" else if(ot< -0.3)"Mono-" else "Flat"}
  else if(sc==1){fc<-which(diff(sign(dy))!=0)[1];if(dy[max(1,fc-1)]>0)"Inv-U" else "U"
  }else if(sc==2)"J/S" else "Complex"
}

# ── Custom legend panel showing both fitted line and 95% CI ribbon ──
make_legend_plot <- function() {
  df <- tibble(
    Source = factor(names(src_labels), levels = names(src_labels)),
    y = c(3.1, 2.05, 1.0),
    xmin = 0.10,
    xmax = 0.38,
    label_x = 0.46,
    label = unname(src_labels)
  ) %>%
    mutate(ymin = y - 0.21, ymax = y + 0.21)

  ggplot(df) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Source),
      alpha = legend_ribbon_alpha, color = NA, show.legend = FALSE
    ) +
    geom_segment(
      aes(x = xmin, xend = xmax, y = y, yend = y, color = Source, linetype = Source),
      linewidth = 1.45, lineend = "round", show.legend = FALSE
    ) +
    geom_text(
      aes(x = label_x, y = y, label = label),
      hjust = 0, size = 3.7, color = "black"
    ) +
    annotate(
      "text", x = 0.10, y = 3.84,
      label = "Line = fitted OR; shaded band = 95% CI",
      hjust = 0, size = 3.8, color = "#333333", fontface = "plain"
    ) +
    annotate(
      "segment", x = 0.71, xend = 0.85, y = 3.84, yend = 3.84,
      colour = "grey45", linetype = "dashed", linewidth = 0.62
    ) +
    annotate(
      "text", x = 0.87, y = 3.84,
      label = "OR = 1", hjust = 0, size = 3.55, color = "#555555"
    ) +
    scale_color_manual(values = src_colors) +
    scale_fill_manual(values = src_colors) +
    scale_linetype_manual(values = src_ltys) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.40, 4.06), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(4, 10, 4, 10),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

make_legend_grob <- function() {
  ggplotGrob(make_legend_plot())
}

get_eval_alt_sources <- function(cv) {
  if (cv == "Humidity") {
    # Primary and ERA5 both use ERA5 humidity, so only MERRA-2 is an independent
    # alternative source for assessing RH robustness.
    return("MERRA2")
  }
  c("ERA5", "MERRA2")
}

get_primary_alt_curve_metrics <- function(ctbl, bt, cv_label, mode, alt_source) {
  x <- ctbl %>%
    filter(
      Bacteria == bt,
      Variable == cv_label,
      Mode == mode,
      ((S1 == "Primary" & S2 == alt_source) | (S1 == alt_source & S2 == "Primary"))
    )
  if (nrow(x) == 0) {
    return(tibble(
      Alt_Source = alt_source,
      Curve_r = NA_real_,
      DirAgree = NA,
      ShapeAgree = NA
    ))
  }
  x %>%
    summarise(
      Alt_Source = alt_source,
      Curve_r = mean(Curve_r, na.rm = TRUE),
      DirAgree = ifelse(all(is.na(DirAgree)), NA, all(DirAgree, na.rm = TRUE)),
      ShapeAgree = ifelse(all(is.na(ShapeAgree)), NA, all(ShapeAgree, na.rm = TRUE))
    )
}

get_primary_alt_pct_metrics <- function(ptbl, bt, cv, mode, alt_source) {
  pri <- ptbl %>% filter(Bacteria == bt, CV == cv, Mode == mode, Source == "Primary")
  alt <- ptbl %>% filter(Bacteria == bt, CV == cv, Mode == mode, Source == alt_source)

  if (nrow(pri) == 0 || nrow(alt) == 0) {
    return(tibble(
      Alt_Source = alt_source,
      Primary_DlogOR = NA_real_,
      Alt_DlogOR = NA_real_,
      Direction_Agree = NA,
      Magnitude_Ratio = NA_real_
    ))
  }

  primary_d <- pri$DlogOR[1]
  alt_d <- alt$DlogOR[1]
  magnitude_ratio <- NA_real_
  if (!is.na(primary_d) && abs(primary_d) > 0.05 && !is.na(alt_d)) {
    magnitude_ratio <- abs(alt_d) / abs(primary_d)
  }

  tibble(
    Alt_Source = alt_source,
    Primary_DlogOR = primary_d,
    Alt_DlogOR = alt_d,
    Direction_Agree = ifelse(is.na(pri$Dir[1]) || is.na(alt$Dir[1]), NA, pri$Dir[1] == alt$Dir[1]),
    Magnitude_Ratio = magnitude_ratio,
    Reversal_Substantive = ifelse(
      is.na(primary_d) || is.na(alt_d) || is.na(pri$Dir[1]) || is.na(alt$Dir[1]),
      NA,
      (pri$Dir[1] != alt$Dir[1]) &&
        min(abs(primary_d), abs(alt_d)) >= reversal_min_secondary_abs_dlogor &&
        max(abs(primary_d), abs(alt_d)) >= reversal_min_abs_dlogor
    )
  )
}

classify_rob <- function(stbl, ctbl, ptbl, bt, cv, mode) {
  sm <- stbl %>% filter(Bacteria == bt, CV == cv, Mode == mode)
  if (nrow(sm) == 0) {
    return(tibble(
      Robustness = "Insufficient data",
      Robustness_Reason = "No smooth-term records available for this bacteria-variable-mode combination.",
      Eval_Alt_Sources = NA_character_,
      Primary_Supported = NA,
      Sig_Retained_Rate = NA_real_,
      Curve_r_Mean = NA_real_,
      Shape_Agreement_Rate = NA_real_,
      Direction_Agreement_Rate = NA_real_,
      Magnitude_Ratio_Mean = NA_real_,
      Reversal_Substantive_Rate = NA_real_
    ))
  }

  pr <- sm %>% filter(Source == "Primary")
  if (nrow(pr) == 0) {
    return(tibble(
      Robustness = "No primary",
      Robustness_Reason = "Primary source model was not available.",
      Eval_Alt_Sources = NA_character_,
      Primary_Supported = NA,
      Sig_Retained_Rate = NA_real_,
      Curve_r_Mean = NA_real_,
      Shape_Agreement_Rate = NA_real_,
      Direction_Agreement_Rate = NA_real_,
      Magnitude_Ratio_Mean = NA_real_,
      Reversal_Substantive_Rate = NA_real_
    ))
  }

  eval_alts <- get_eval_alt_sources(cv)
  alt_rows <- sm %>% filter(Source %in% eval_alts)
  cv_label <- cv2label[cv]

  curve_eval <- bind_rows(lapply(eval_alts, function(src) {
    get_primary_alt_curve_metrics(ctbl, bt, cv_label, mode, src)
  }))

  pct_eval <- bind_rows(lapply(eval_alts, function(src) {
    get_primary_alt_pct_metrics(ptbl, bt, cv, mode, src)
  }))

  primary_supported <- isTRUE(pr$Sig[1]) && !isTRUE(pr$PenOut[1])
  alt_supported <- alt_rows %>%
    mutate(Supported = Sig & !PenOut)

  sig_retained_rate <- if (nrow(alt_supported) > 0) {
    mean(alt_supported$Supported, na.rm = TRUE)
  } else {
    NA_real_
  }

  curve_r_mean <- if (nrow(curve_eval) > 0 && any(!is.na(curve_eval$Curve_r))) {
    mean(curve_eval$Curve_r, na.rm = TRUE)
  } else {
    NA_real_
  }

  shape_agreement_rate <- if (nrow(curve_eval) > 0 && any(!is.na(curve_eval$ShapeAgree))) {
    mean(curve_eval$ShapeAgree, na.rm = TRUE)
  } else {
    NA_real_
  }

  direction_agreement_rate <- if (nrow(pct_eval) > 0 && any(!is.na(pct_eval$Direction_Agree))) {
    mean(pct_eval$Direction_Agree, na.rm = TRUE)
  } else {
    NA_real_
  }

  magnitude_ratio_mean <- if (nrow(pct_eval) > 0 && any(!is.na(pct_eval$Magnitude_Ratio))) {
    mean(pct_eval$Magnitude_Ratio, na.rm = TRUE)
  } else {
    NA_real_
  }

  reversal_substantive_rate <- if (nrow(pct_eval) > 0 && any(!is.na(pct_eval$Reversal_Substantive))) {
    mean(pct_eval$Reversal_Substantive, na.rm = TRUE)
  } else {
    NA_real_
  }

  robustness <- "Attenuated"
  reason <- "Primary association is retained only partially across alternative source specifications."

  if (!primary_supported) {
    if (all(!(alt_supported$Supported), na.rm = TRUE)) {
      robustness <- "Robust (consistently NS)"
      reason <- "Primary model was non-significant or penalized out, and all independent alternative sources remained non-significant."
    } else {
      robustness <- "Product-sensitive (emergent)"
      reason <- "Primary model was non-significant, but at least one independent alternative source produced a retained significant association."
    }
  } else if (all(!(alt_supported$Supported), na.rm = TRUE)) {
    robustness <- "Product-sensitive (lost)"
    reason <- "Primary model was supported, but all independent alternative sources lost the association."
  } else if (
    !is.na(direction_agreement_rate) && direction_agreement_rate < 1 &&
    !is.na(reversal_substantive_rate) && reversal_substantive_rate > 0
  ) {
    robustness <- "Product-sensitive (reversal)"
    reason <- "The high-versus-low exposure contrast reversed direction in at least one independent alternative source with a non-trivial percentile-based effect magnitude."
  } else if (
    !is.na(sig_retained_rate) && sig_retained_rate == 1 &&
    (is.na(curve_r_mean) || curve_r_mean >= robust_curve_r_threshold) &&
    (is.na(shape_agreement_rate) || shape_agreement_rate >= robust_shape_agreement_threshold) &&
    (is.na(magnitude_ratio_mean) || (magnitude_ratio_mean >= robust_mag_ratio_lower && magnitude_ratio_mean <= robust_mag_ratio_upper))
  ) {
    robustness <- "Robust"
    reason <- "All independent alternative sources retained the association, with high overall curve agreement and comparable percentile-based effect magnitude."
  } else if (
    !is.na(sig_retained_rate) && sig_retained_rate >= 0.5 &&
    (is.na(curve_r_mean) || curve_r_mean >= minor_curve_r_threshold) &&
    (is.na(shape_agreement_rate) || shape_agreement_rate >= minor_shape_agreement_threshold) &&
    (is.na(magnitude_ratio_mean) || (magnitude_ratio_mean >= minor_mag_ratio_lower && magnitude_ratio_mean <= minor_mag_ratio_upper))
  ) {
    robustness <- "Robust (minor variation)"
    reason <- "The association was retained in at least one independent alternative source and the overall pattern remained broadly consistent, with moderate variation in magnitude or shape."
  }

  tibble(
    Robustness = robustness,
    Robustness_Reason = reason,
    Eval_Alt_Sources = paste(eval_alts, collapse = ", "),
    Primary_Supported = primary_supported,
    Sig_Retained_Rate = sig_retained_rate,
    Curve_r_Mean = curve_r_mean,
    Shape_Agreement_Rate = shape_agreement_rate,
    Direction_Agreement_Rate = direction_agreement_rate,
    Magnitude_Ratio_Mean = magnitude_ratio_mean,
    Reversal_Substantive_Rate = reversal_substantive_rate
  )
}

format_bridge_tag_components <- function(tag) {
  parts <- strsplit(tag, "_", fixed = TRUE)[[1]]
  source <- parts[1]
  sample <- parts[2]

  source_label <- switch(
    source,
    "Primary" = "Primary",
    "ERA5" = "ERA5",
    "MERRA2" = "MERRA-2",
    source
  )

  sample_label <- switch(
    sample,
    "full" = "Full",
    "common" = "Common",
    sample
  )

  c(source = source_label, sample = sample_label)
}

draw_bridge_figure_table <- function(df, pdf_path, png_path = NULL) {
  parsed <- t(vapply(df$Tag, format_bridge_tag_components, character(2)))
  df$Source <- parsed[, "source"]
  df$Sample <- parsed[, "sample"]

  bacteria_order <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
  source_order <- c("Primary", "ERA5", "MERRA-2")
  sample_order <- c("Full", "Common")

  tab <- df %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_order),
      Source = factor(Source, levels = source_order),
      Sample = factor(Sample, levels = sample_order)
    ) %>%
    arrange(Bacteria, Source, Sample) %>%
    mutate(
      Bacteria = as.character(Bacteria),
      Source = as.character(Source),
      Sample = as.character(Sample),
      Bacteria_Display = Bacteria,
      N_txt = ifelse(is.na(N), "N/A", as.character(N)),
      Countries_txt = ifelse(is.na(Countries), "N/A", as.character(Countries)),
      AIC_txt = ifelse(is.na(AIC), "N/A", format(round(AIC, 2), nsmall = 2)),
      Dev_txt = ifelse(is.na(DevExpl), "N/A", paste0(format(round(DevExpl * 100, 2), nsmall = 2), "%")),
      Rsq_txt = ifelse(is.na(Rsq), "N/A", format(round(Rsq, 2), nsmall = 2))
    )

  for (i in seq_len(nrow(tab))) {
    if (i > 1 && identical(tab$Bacteria[i], tab$Bacteria[i - 1])) {
      tab$Bacteria_Display[i] <- ""
    }
  }

  col_names <- c("AMR Phenotype", "Climate Source", "Sample", "N", "Countries", "AIC", "Explained Deviance (%)", "Adjusted R²")
  col_widths <- c(0.18, 0.14, 0.10, 0.07, 0.10, 0.12, 0.17, 0.12)
  pdf_width <- 11.69
  pdf_height <- 8.27

  make_plot <- function(device = c("pdf", "png"), file_path) {
    device <- match.arg(device)
    if (device == "pdf") {
      pdf(file_path, width = pdf_width, height = pdf_height, useDingbats = FALSE, family = "sans")
    } else {
      png(file_path, width = pdf_width, height = pdf_height, units = "in", res = 300, type = "cairo")
    }

    grid.newpage()
    pushViewport(viewport(x = 0.5, y = 0.5, width = 0.985, height = 0.965))

    title_y <- 0.975
    grid.text(
      "Reanalysis Sensitivity Analysis: Model Performance Bridge Table",
      x = 0.5, y = title_y, just = "center",
      gp = gpar(fontface = "bold", fontsize = 14, col = "#333333")
    )
    grid.text(
      "Primary exposure set vs ERA5 vs MERRA-2 across full and common samples",
      x = 0.5, y = title_y - 0.035, just = "center",
      gp = gpar(fontsize = 11, col = "#333333")
    )

    table_top <- 0.88
    table_bottom <- 0.12
    header_height <- 0.06
    table_width <- 0.965
    row_height <- (table_top - table_bottom - header_height) / nrow(tab)

    col_positions <- cumsum(c((1 - table_width) / 2, col_widths * table_width))
    left_x <- col_positions[1]
    right_x <- col_positions[length(col_positions)]
    bottom_y <- table_top - header_height - nrow(tab) * row_height

    grid.rect(x = 0.5, y = (table_top + bottom_y) / 2, width = table_width, height = table_top - bottom_y,
              gp = gpar(fill = NA, col = "#8EBAE5", lwd = 2))
    grid.rect(x = 0.5, y = table_top - header_height / 2, width = table_width, height = header_height,
              gp = gpar(fill = "#E6F3FF", col = "#8EBAE5", lwd = 1))

    for (i in seq_along(col_names)) {
      mid <- (col_positions[i] + col_positions[i + 1]) / 2
      grid.text(col_names[i], x = mid, y = table_top - header_height / 2,
                gp = gpar(fontface = "bold", fontsize = 9.8, col = "#333333"))
    }

    grid.lines(x = c(left_x, right_x), y = c(table_top - header_height, table_top - header_height),
               gp = gpar(col = "#8EBAE5", lwd = 1.5))

    for (i in 2:(length(col_positions) - 1)) {
      grid.lines(x = c(col_positions[i], col_positions[i]), y = c(table_top, bottom_y),
                 gp = gpar(col = "#8EBAE5", lwd = 1))
    }

    for (i in seq_len(nrow(tab))) {
      y_pos <- table_top - header_height - (i - 0.5) * row_height
      if (i %% 2 == 1) {
        grid.rect(x = 0.5, y = y_pos, width = table_width, height = row_height,
                  gp = gpar(fill = "#F5F9FC", col = NA))
      }

      row_values <- c(
        tab$Bacteria_Display[i],
        tab$Source[i],
        tab$Sample[i],
        tab$N_txt[i],
        tab$Countries_txt[i],
        tab$AIC_txt[i],
        tab$Dev_txt[i],
        tab$Rsq_txt[i]
      )

      for (j in seq_along(row_values)) {
        mid <- (col_positions[j] + col_positions[j + 1]) / 2
        grid.text(
          row_values[j],
          x = mid, y = y_pos, just = "center",
          gp = gpar(
            fontsize = 9.1,
            fontface = if (j == 1 && nzchar(row_values[j])) "bold" else "plain",
            col = "#333333"
          )
        )
      }

      if (i < nrow(tab)) {
        grid.lines(x = c(left_x, right_x),
                   y = c(table_top - header_height - i * row_height, table_top - header_height - i * row_height),
                   gp = gpar(col = "#8EBAE5", lwd = 1))
      }
    }

    grid.text(
      "Note: Primary = CRU TS temperature, precipitation, and wet days + ERA5 relative humidity.",
      x = 0.5, y = bottom_y - 0.038, just = "center",
      gp = gpar(fontface = "italic", fontsize = 8.8, col = "#333333")
    )
    grid.text(
      "Common sample = observations with complete climate exposure coverage across Primary, ERA5, and MERRA-2.",
      x = 0.5, y = bottom_y - 0.064, just = "center",
      gp = gpar(fontface = "italic", fontsize = 8.8, col = "#333333")
    )

    popViewport()
    dev.off()
  }

  make_plot("pdf", pdf_path)
  if (!is.null(png_path)) {
    make_plot("png", png_path)
  }
}

draw_robustness_figure_table <- function(df, pdf_path, png_path = NULL, mode_label = "common") {
  tab <- df %>%
    filter(Mode == mode_label) %>%
    mutate(
      Primary_F_txt = ifelse(is.na(PrimaryF), "N/A", format(round(PrimaryF, 2), nsmall = 2)),
      ERA5_F_txt = ifelse(is.na(ERA5F), "N/A", format(round(ERA5F, 2), nsmall = 2)),
      MERRA2_F_txt = ifelse(is.na(MERRA2F), "N/A", format(round(MERRA2F, 2), nsmall = 2)),
      Curve_r_txt = ifelse(is.na(Curve_r_Mean), "N/A", format(round(Curve_r_Mean, 2), nsmall = 2)),
      Dir_txt = ifelse(is.na(DirectionAgreementRate), "N/A", format(round(DirectionAgreementRate, 2), nsmall = 2)),
      Bacteria_Display = Bacteria
    )

  if (nrow(tab) == 0) {
    return(invisible(NULL))
  }

  for (i in seq_len(nrow(tab))) {
    if (i > 1 && identical(tab$Bacteria[i], tab$Bacteria[i - 1])) {
      tab$Bacteria_Display[i] <- ""
    }
  }

  col_names <- c("AMR Phenotype", "Climate Variable", "Primary F", "ERA5 F", "MERRA-2 F", "Curve r", "Dir. Agree.", "Robustness")
  col_widths <- c(0.16, 0.14, 0.09, 0.09, 0.10, 0.09, 0.10, 0.23)
  pdf_width <- 12
  pdf_height <- 8.5

  make_plot <- function(device = c("pdf", "png"), file_path) {
    device <- match.arg(device)
    if (device == "pdf") {
      pdf(file_path, width = pdf_width, height = pdf_height, useDingbats = FALSE, family = "sans")
    } else {
      png(file_path, width = pdf_width, height = pdf_height, units = "in", res = 300, type = "cairo")
    }

    grid.newpage()
    pushViewport(viewport(x = 0.5, y = 0.5, width = 0.985, height = 0.965))

    title_y <- 0.975
    grid.text(
      paste0("Reanalysis Sensitivity Robustness Summary (", tools::toTitleCase(mode_label), " Sample)"),
      x = 0.5, y = title_y, just = "center",
      gp = gpar(fontface = "bold", fontsize = 14, col = "#333333")
    )
    grid.text(
      "Rule-based robustness classification using primary-vs-alternative smooth support, curve agreement, and percentile contrast stability",
      x = 0.5, y = title_y - 0.032, just = "center",
      gp = gpar(fontsize = 10.5, col = "#333333")
    )

    table_top <- 0.89
    table_bottom <- 0.12
    header_height <- 0.06
    table_width <- 0.97
    row_height <- (table_top - table_bottom - header_height) / nrow(tab)

    col_positions <- cumsum(c((1 - table_width) / 2, col_widths * table_width))
    left_x <- col_positions[1]
    right_x <- col_positions[length(col_positions)]
    bottom_y <- table_top - header_height - nrow(tab) * row_height

    grid.rect(x = 0.5, y = (table_top + bottom_y) / 2, width = table_width, height = table_top - bottom_y,
              gp = gpar(fill = NA, col = "#8EBAE5", lwd = 2))
    grid.rect(x = 0.5, y = table_top - header_height / 2, width = table_width, height = header_height,
              gp = gpar(fill = "#E6F3FF", col = "#8EBAE5", lwd = 1))

    for (i in seq_along(col_names)) {
      mid <- (col_positions[i] + col_positions[i + 1]) / 2
      grid.text(col_names[i], x = mid, y = table_top - header_height / 2,
                gp = gpar(fontface = "bold", fontsize = 9.4, col = "#333333"))
    }

    grid.lines(x = c(left_x, right_x), y = c(table_top - header_height, table_top - header_height),
               gp = gpar(col = "#8EBAE5", lwd = 1.5))

    for (i in 2:(length(col_positions) - 1)) {
      grid.lines(x = c(col_positions[i], col_positions[i]), y = c(table_top, bottom_y),
                 gp = gpar(col = "#8EBAE5", lwd = 1))
    }

    for (i in seq_len(nrow(tab))) {
      y_pos <- table_top - header_height - (i - 0.5) * row_height

      if (i %% 2 == 1) {
        grid.rect(x = 0.5, y = y_pos, width = table_width, height = row_height,
                  gp = gpar(fill = "#F5F9FC", col = NA))
      }

      row_values <- c(
        tab$Bacteria_Display[i],
        tab$CV[i],
        tab$Primary_F_txt[i],
        tab$ERA5_F_txt[i],
        tab$MERRA2_F_txt[i],
        tab$Curve_r_txt[i],
        tab$Dir_txt[i],
        tab$Robustness[i]
      )

      for (j in seq_along(row_values)) {
        mid <- (col_positions[j] + col_positions[j + 1]) / 2
        grid.text(
          row_values[j],
          x = mid, y = y_pos, just = "center",
          gp = gpar(
            fontsize = 8.6,
            fontface = if (j == 1 && nzchar(row_values[j])) "bold" else "plain",
            col = "#333333"
          )
        )
      }

      if (i < nrow(tab)) {
        grid.lines(x = c(left_x, right_x),
                   y = c(table_top - header_height - i * row_height, table_top - header_height - i * row_height),
                   gp = gpar(col = "#8EBAE5", lwd = 1))
      }
    }

    grid.text(
      "Note: `Curve r` = mean Pearson correlation between Primary and independent alternative curves; `Dir. Agree.` = proportion of alternative sources preserving the percentile-contrast direction.",
      x = 0.5, y = bottom_y - 0.04, just = "center",
      gp = gpar(fontface = "italic", fontsize = 8.4, col = "#333333")
    )

    popViewport()
    dev.off()
  }

  make_plot("pdf", pdf_path)
  if (!is.null(png_path)) {
    make_plot("png", png_path)
  }
}

cat("Part 0 complete.\n")


# Part 1: Run all models
cat("\n",strrep("=",70),"\n PART 1: Model Fitting\n",strrep("=",70),"\n\n")

R <- list()

for (spec in bact_specs) {
  bt<-spec$title; fp<-file.path(input_data_dir,spec$fn)
  if(!file.exists(fp)){cat("  SKIP",bt,"\n");next}
  cat("\n  ===",bt,"===\n"); R[[bt]]<-list()

  for (src in names(clim_cols)) {
    R[[bt]][[src]] <- list()
    for (mode in c("full","common")) {
      tag <- paste0(src,"_",mode)
      cat("    ",tag,"... ")
      tryCatch({
        prep<-prepare_data(fp,bt,src,mode); md<-prep$data; sp<-prep$scale_params
        model<-build_gamm(md,bt)

        # ── Save model summary ──
        save_model_summary(model, bt, src, mode, prep, dirs$summaries)

        ss<-get_smooth(model); ss$Bacteria<-bt; ss$Source<-src; ss$Mode<-mode
        cvs<-list(); pcs<-list()
        for(vn in clim_vars){
          c<-extract_curve(model,md,sp,vn,bt); c$Source<-src; c$Bacteria<-bt; c$Variable<-vn; c$Mode<-mode
          cvs[[vn]]<-c
          pc<-pct_contrast(c,md,vn)
          if(!is.null(pc)){pc$Source<-src;pc$Bacteria<-bt;pc$Mode<-mode
            pc$CV<-case_when(grepl("TMP",vn)~"Temperature",grepl("PREC",vn)~"Precipitation",
              grepl("HUM",vn)~"Humidity",grepl("WET",vn)~"WetDays"); pcs[[vn]]<-pc}}
        sm<-summary(model)
        pf<-tibble(Bacteria=bt,Source=src,Mode=mode,AIC=AIC(model),DevExpl=sm$dev.expl,
          Rsq=sm$r.sq,N=prep$n_obs,Countries=prep$n_countries)
        R[[bt]][[src]][[mode]]<-list(model=model,data=md,sp=sp,smooth=ss,
          curves=cvs,pct=bind_rows(pcs),perf=pf)
        cat("N=",prep$n_obs,"(",prep$n_countries,"ctry) Dev=",round(sm$dev.expl,4),"\n")
      },error=function(e){
        msg<-conditionMessage(e); cat("ERR:",msg,"\n")
        error_log[[length(error_log)+1]] <<- tibble(Bacteria=bt,Source=src,Mode=mode,Error=msg)
        R[[bt]][[src]][[mode]]<<-NULL})
    }
  }
}
cat("\n  All models fitted.\n")

# Collect
coll <- function(what) {
  bind_rows(lapply(names(R),function(bt)
    bind_rows(lapply(names(R[[bt]]),function(src)
      bind_rows(lapply(names(R[[bt]][[src]]),function(mode){
        r<-R[[bt]][[src]][[mode]]; if(!is.null(r)) r[[what]] else NULL}))))))}

perf_all <- coll("perf"); smooth_all <- coll("smooth"); pct_all <- coll("pct")
write.csv(perf_all,file.path(dirs$model,"performance.csv"),row.names=FALSE)
write.csv(smooth_all,file.path(dirs$model,"smooth_terms.csv"),row.names=FALSE)
write.csv(pct_all,file.path(dirs$curve,"percentile_contrast.csv"),row.names=FALSE)

# Save error log
if(length(error_log)>0) write.csv(bind_rows(error_log),file.path(output_root,"error_log.csv"),row.names=FALSE)


# Part 2: Exposure comparison
cat("\n",strrep("=",70),"\n PART 2: Exposure Comparison\n",strrep("=",70),"\n\n")

vm <- list(TMP=c(Primary="TMP",ERA5="ERA5_annual_mean_t2m_c",MERRA2="MERRA2_annual_mean_t2m_c"),
           PREC=c(Primary="PREC",ERA5="ERA5_annual_precip_mm",MERRA2="MERRA2_annual_precip_mm"),
           WET=c(Primary="WET",ERA5="ERA5_wet_days",MERRA2="MERRA2_wet_days"),
           HUM=c(Primary="HUM",ERA5="HUM",MERRA2="MERRA2_annual_mean_rh_pct"))
pairs <- list(c("Primary","ERA5"),c("Primary","MERRA2"),c("ERA5","MERRA2"))

# ── Layer 1: Raw pairwise overlap ──
exp_list <- list()
exp_point_list <- list()
for(spec in bact_specs){
  bt<-spec$title; fp<-file.path(input_data_dir,spec$fn); if(!file.exists(fp))next
  raw<-read.csv(fp)
  for(vn in names(vm)) for(pr in pairs){
    c1<-vm[[vn]][pr[1]]; c2<-vm[[vn]][pr[2]]
    if(!all(c(c1,c2)%in%names(raw)))next
    x<-as.numeric(raw[[c1]]); y<-as.numeric(raw[[c2]])
    ok<-!is.na(x)&!is.na(y)&is.finite(x)&is.finite(y)
    if(sum(ok)<10) next
    raw_ok <- raw[ok,,drop=FALSE]
    x<-x[ok];y<-y[ok]
    exp_list[[length(exp_list)+1]]<-tibble(
      Bacteria=bt,Variable=clim_labels[vn],S1=pr[1],S2=pr[2],Layer="Raw",
      N_overlap=length(x),r=cor(x,y),rho=cor(x,y,method="spearman"),
      Bias=mean(y-x),RMSE=sqrt(mean((y-x)^2)),
      Same_Source=is_rh_same(clim_labels[vn],pr[1],pr[2]))
    exp_point_list[[length(exp_point_list)+1]] <- tibble(
      Bacteria=bt,
      Variable=clim_labels[vn],
      S1=pr[1],
      S2=pr[2],
      Pair=pair_label(pr[1], pr[2]),
      Source1Value=x,
      Source2Value=y,
      MeanValue=(x+y)/2,
      Diff=(y-x),
      Year=if("year"%in%names(raw_ok)) as.numeric(as.character(raw_ok$year)) else NA_real_,
      Latitude=if("lat"%in%names(raw_ok)) as.numeric(raw_ok$lat) else NA_real_,
      AbsLatitude=if("lat"%in%names(raw_ok)) abs(as.numeric(raw_ok$lat)) else NA_real_,
      Country=if("NAME"%in%names(raw_ok)) as.character(raw_ok$NAME) else NA_character_,
      Same_Source=is_rh_same(clim_labels[vn],pr[1],pr[2]))
  }
}
exp_df <- bind_rows(exp_list)
exp_point_df <- bind_rows(exp_point_list)
write.csv(exp_df,file.path(dirs$exposure,"exposure_agreement_raw.csv"),row.names=FALSE)

# ── Layer 2: Analysis-sample exposure agreement (post prepare_data) ──
exp_analysis <- list()
for(bt in names(R)) for(src in names(R[[bt]])){
  r <- R[[bt]][[src]][["full"]]
  if(is.null(r)) next
  md <- r$data
  for(vn_short in c("TMP","PREC","HUM","WET")){
    oc <- paste0(vn_short,"_orig")
    if(!oc %in% names(md)) next
    exp_analysis[[length(exp_analysis)+1]] <- tibble(
      Bacteria=bt, Source=src, Variable=clim_labels[vn_short],
      N=nrow(md), Mean=mean(md[[oc]],na.rm=T), SD=sd(md[[oc]],na.rm=T),
      P10=quantile(md[[oc]],0.1,na.rm=T), P90=quantile(md[[oc]],0.9,na.rm=T))
  }
}
exp_anal_df <- bind_rows(exp_analysis)
write.csv(exp_anal_df,file.path(dirs$exposure,"exposure_analysis_sample_summary.csv"),row.names=FALSE)

# ── Heatmap 1: Pearson r (with N_overlap) ──
if(nrow(exp_df)>0){
  hd<-exp_df%>%mutate(
    Pair=factor(mapply(pair_label, S1, S2), levels=pair_levels),
    Bacteria=factor(Bacteria, levels=bacteria_order),
    Variable=factor(Variable, levels=curve_variable_order),
    disp=ifelse(Same_Source,"same\nsource",paste0(sprintf("%.3f",r),"\n(N=",N_overlap,")")),
    fr=ifelse(Same_Source,NA_real_,r))
  p1<-ggplot(hd,aes(x=Pair,y=Bacteria,fill=fr))+
    geom_tile(color="white",linewidth=0.8)+
    geom_text(aes(label=disp),size=2.5,fontface=ifelse(hd$Same_Source,"italic","plain"))+
    facet_wrap(~Variable,nrow=1)+
    scale_fill_gradient2(low="#BC4B51",mid="#F7F4EF",high="#2C7C6D",midpoint=0.95,
      limits=c(0.8,1),name="Pearson r",oob=squish,na.value="grey85")+
    labs(title="Level 1a: Exposure Agreement \u2014 Correlation",
         subtitle="Grey = same source by design | N = pairwise overlap",x=NULL,y=NULL)+
    theme_s+theme(
      axis.text.x=element_text(size=8, angle=20, hjust=1),
      panel.grid=element_blank(),
      strip.background=element_rect(fill="#F3F5F7", colour=NA),
      legend.position="right")
  ggsave(file.path(dirs$exposure,"heatmap_correlation.pdf"),p1,width=14,height=5,dpi=300)
  ggsave(file.path(dirs$exposure,"heatmap_correlation.png"),p1,width=14,height=5,dpi=300)

  # ── Heatmap 2: Bias & RMSE ──
  hd2 <- exp_df %>%
    filter(!Same_Source) %>%
    mutate(
      Pair=factor(mapply(pair_label, S1, S2), levels=pair_levels),
      Bacteria=factor(Bacteria, levels=bacteria_order),
      Variable=factor(Variable, levels=curve_variable_order))
  p2<-ggplot(hd2,aes(x=Pair,y=Bacteria))+
    geom_tile(aes(fill=RMSE),color="white",linewidth=1.0)+
    geom_point(aes(size=abs(Bias),color=Bias),shape=16,alpha=0.94)+
    facet_wrap(~Variable,nrow=1,scales="free")+
    scale_fill_gradientn(colours=rmse_fill_palette,trans="sqrt",name="RMSE")+
    scale_color_gradient2(low=bias_div_palette[1],mid=bias_div_palette[2],high=bias_div_palette[3],
      midpoint=0,name="Mean bias")+
    scale_size_continuous(range=c(1.8,7.2),name="|Mean bias|")+
    labs(title="Level 1b: Exposure Agreement \u2014 Bias & RMSE",
         subtitle="Tile fill = RMSE; point color and size = signed and absolute mean bias on the overlap sample",
         x=NULL,y=NULL)+
    guides(fill=guide_colorbar(order=1),color=guide_colorbar(order=2),size=guide_legend(order=3))+
    theme_s+theme(
      axis.text.x=element_text(size=8, angle=20, hjust=1),
      panel.grid=element_blank(),
      strip.background=element_rect(fill="#F3F5F7", colour=NA),
      legend.position="right",
      legend.box="vertical")
  ggsave(file.path(dirs$exposure,"heatmap_bias_rmse.pdf"),p2,width=14,height=5.6,dpi=300,bg="white")
  ggsave(file.path(dirs$exposure,"heatmap_bias_rmse.png"),p2,width=14,height=5.6,dpi=300,bg="white")

  # ── Alternative view: correlation lollipop plot ──
  corr_alt <- exp_df %>%
    filter(!Same_Source) %>%
    mutate(
      Pair = factor(mapply(pair_label, S1, S2), levels = pair_levels),
      Bacteria = factor(Bacteria, levels = bacteria_order)
    ) %>%
    pivot_longer(c(r, rho), names_to = "Metric", values_to = "Correlation") %>%
    mutate(
      Metric = factor(Metric, levels = c("r", "rho"), labels = c("Pearson r", "Spearman \u03c1")),
      Variable = factor(Variable, levels = curve_variable_order)
    )

  p_corr_alt <- ggplot(corr_alt, aes(x = Bacteria, y = Correlation, color = Pair, shape = Pair)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95, ymax = 1.00,
             fill = "#ECECEC", alpha = 0.55) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "#A9A9A9", linewidth = 0.35) +
    geom_linerange(
      aes(ymin = 0.80, ymax = Correlation),
      position = position_dodge(width = 0.68),
      linewidth = 0.62, alpha = 0.80
    ) +
    geom_point(
      position = position_dodge(width = 0.68),
      size = 2.6, stroke = 0.8, fill = "white"
    ) +
    coord_flip() +
    facet_grid(Metric ~ Variable) +
    scale_color_manual(values = pair_colors, drop = FALSE) +
    scale_shape_manual(values = c(21, 22, 24), drop = FALSE) +
    scale_y_continuous(limits = c(0.80, 1.00), breaks = c(0.80, 0.90, 0.95, 1.00)) +
    labs(
      title = "Alternative View: Cross-product Exposure Agreement",
      subtitle = "Pearson captures linear agreement, whereas Spearman captures rank-order consistency; the shaded band marks very high agreement (r \u2265 0.95)",
      x = NULL, y = "Correlation coefficient", color = "Product pair", shape = "Product pair"
    ) +
    theme_s +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "#E9E9E9"),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
  ggsave(file.path(dirs$exposure,"correlation_lollipop.pdf"), p_corr_alt, width = 14, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(dirs$exposure,"correlation_lollipop.png"), p_corr_alt, width = 14, height = 7, dpi = 300, bg = "white")

  bias_rmse_profile <- hd2 %>%
    select(Bacteria, Variable, Pair, Bias, RMSE) %>%
    pivot_longer(c(RMSE, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = factor(Metric, levels = c("RMSE", "Bias"), labels = c("RMSE", "Mean bias")),
      Bacteria = factor(Bacteria, levels = bacteria_order),
      Variable = factor(Variable, levels = curve_variable_order),
      xmin = pmin(0, Value),
      xmax = pmax(0, Value)
    )

  p_bias_profile <- ggplot(bias_rmse_profile, aes(y = Bacteria, x = Value, color = Pair, shape = Pair)) +
    geom_vline(xintercept = 0, color = "#A9A9A9", linewidth = 0.38) +
    geom_linerange(
      aes(xmin = xmin, xmax = xmax),
      position = position_dodge(width = 0.68),
      linewidth = 0.72, alpha = 0.82
    ) +
    geom_point(
      position = position_dodge(width = 0.68),
      size = 2.7, stroke = 0.85, fill = "white"
    ) +
    facet_grid(Metric ~ Variable, scales = "free_x") +
    scale_color_manual(values = pair_colors, drop = FALSE) +
    scale_shape_manual(values = c(21, 22, 24), drop = FALSE) +
    labs(
      title = "Alternative View: Magnitude and Direction of Product Differences",
      subtitle = "Top row shows RMSE and bottom row shows signed mean bias on the overlap sample; free x-scales allow small and large discrepancies to be seen together",
      x = "Difference metric value", y = NULL, color = "Product pair", shape = "Product pair"
    ) +
    theme_s +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "#E9E9E9"),
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "#F3F5F7", colour = NA),
      legend.position = "bottom"
    )
  ggsave(file.path(dirs$exposure, "bias_rmse_profile.pdf"), p_bias_profile, width = 14, height = 7.2, dpi = 300, bg = "white")
  ggsave(file.path(dirs$exposure, "bias_rmse_profile.png"), p_bias_profile, width = 14, height = 7.2, dpi = 300, bg = "white")

  if (nrow(exp_point_df) > 0) {
    exp_point_plot <- exp_point_df %>%
      filter(!Same_Source) %>%
      mutate(
        Pair = factor(Pair, levels = pair_levels),
        Variable = factor(Variable, levels = curve_variable_order)
      )

    exp_scatter <- exp_point_plot %>%
      group_by(Variable, Pair) %>%
      group_modify(~ if (nrow(.x) > 2200) dplyr::slice_sample(.x, n = 2200) else .x) %>%
      ungroup()

    p_scatter <- ggplot(exp_scatter, aes(x = Source1Value, y = Source2Value)) +
      geom_bin2d(bins = 28) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white", linewidth = 0.55) +
      facet_grid(Variable ~ Pair, scales = "free") +
      scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Point count") +
      labs(
        title = "Alternative View: Pairwise Exposure Scatter-Density",
        subtitle = "Sampled overlap rows; dashed line denotes perfect 1:1 agreement between the paired products",
        x = "Source 1 value", y = "Source 2 value"
      ) +
      theme_s +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#F3F5F7", colour = NA),
        legend.position = "right"
      )
    ggsave(file.path(dirs$exposure, "scatter_density_products.pdf"), p_scatter, width = 13.5, height = 9.5, dpi = 300, bg = "white")
    ggsave(file.path(dirs$exposure, "scatter_density_products.png"), p_scatter, width = 13.5, height = 9.5, dpi = 300, bg = "white")

    bias_time <- exp_point_plot %>%
      filter(!is.na(Year)) %>%
      group_by(Variable, Pair, Year) %>%
      summarise(
        MedianDiff = median(Diff, na.rm = TRUE),
        Q25 = quantile(Diff, 0.25, na.rm = TRUE),
        Q75 = quantile(Diff, 0.75, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    write.csv(bias_time, file.path(dirs$exposure, "bias_over_time_summary.csv"), row.names = FALSE)

    p_bias_time <- ggplot(bias_time, aes(x = Year, y = MedianDiff, color = Pair, fill = Pair)) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45) +
      geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = pair_ribbon_alpha, color = NA) +
      geom_line(linewidth = 0.95) +
      geom_point(aes(size = N), alpha = 0.76) +
      facet_wrap(~Variable, ncol = 2, scales = "free_y") +
      scale_color_manual(values = pair_colors, drop = FALSE) +
      scale_fill_manual(values = pair_colors, drop = FALSE) +
      scale_size_continuous(range = c(0.8, 2.4), guide = "none") +
      labs(
        title = "Alternative View: Bias Stability Over Time",
        subtitle = "Lines show annual median differences (source 2 \u2212 source 1); ribbons denote the interquartile range across overlap rows",
        x = "Year", y = "Difference in exposure value", color = "Product pair", fill = "Product pair"
      ) +
      theme_s +
      theme(
        panel.grid.major.y = element_line(linewidth = 0.25, color = "#E8E8E8"),
        panel.grid.major.x = element_line(linewidth = 0.25, color = "#EFEFEF"),
        strip.background = element_rect(fill = "#F3F5F7", colour = NA)
      )
    ggsave(file.path(dirs$exposure, "bias_over_time.pdf"), p_bias_time, width = 13.2, height = 7.6, dpi = 300, bg = "white")
    ggsave(file.path(dirs$exposure, "bias_over_time.png"), p_bias_time, width = 13.2, height = 7.6, dpi = 300, bg = "white")

    bias_lat <- exp_point_plot %>%
      filter(!is.na(AbsLatitude)) %>%
      mutate(LatBandMid = pmin(floor(AbsLatitude / 10) * 10 + 5, 85)) %>%
      group_by(Variable, Pair, LatBandMid) %>%
      summarise(
        MedianDiff = median(Diff, na.rm = TRUE),
        Q25 = quantile(Diff, 0.25, na.rm = TRUE),
        Q75 = quantile(Diff, 0.75, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    write.csv(bias_lat, file.path(dirs$exposure, "bias_by_latitude_summary.csv"), row.names = FALSE)

    p_bias_lat <- ggplot(bias_lat, aes(x = LatBandMid, y = MedianDiff, color = Pair, fill = Pair)) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45) +
      geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = pair_ribbon_alpha, color = NA) +
      geom_line(linewidth = 0.95) +
      geom_point(aes(size = N), alpha = 0.76) +
      facet_wrap(~Variable, ncol = 2, scales = "free_y") +
      scale_color_manual(values = pair_colors, drop = FALSE) +
      scale_fill_manual(values = pair_colors, drop = FALSE) +
      scale_size_continuous(range = c(0.8, 2.4), guide = "none") +
      scale_x_continuous(
        breaks = seq(5, 85, by = 10),
        labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90")
      ) +
      labs(
        title = "Alternative View: Bias Profile Across Latitude Bands",
        subtitle = "Lines show median differences within 10\u00b0 absolute-latitude bands; ribbons denote the interquartile range",
        x = "Absolute latitude band (\u00b0)", y = "Difference in exposure value", color = "Product pair", fill = "Product pair"
      ) +
      theme_s +
      theme(
        panel.grid.major.y = element_line(linewidth = 0.25, color = "#E8E8E8"),
        panel.grid.major.x = element_line(linewidth = 0.25, color = "#EFEFEF"),
        strip.background = element_rect(fill = "#F3F5F7", colour = NA)
      )
    ggsave(file.path(dirs$exposure, "bias_by_latitude.pdf"), p_bias_lat, width = 13.2, height = 7.6, dpi = 300, bg = "white")
    ggsave(file.path(dirs$exposure, "bias_by_latitude.png"), p_bias_lat, width = 13.2, height = 7.6, dpi = 300, bg = "white")
  }
}


# Part 3: Model performance and bridge table
cat("\n",strrep("=",70),"\n PART 3: Model Performance\n",strrep("=",70),"\n\n")

pf <- perf_all %>% filter(Mode=="full")
if(nrow(pf)>0){
  p_perf<-pf%>%
    pivot_longer(c(AIC,DevExpl,Rsq),names_to="Metric",values_to="Value")%>%
    mutate(Metric=recode(Metric,AIC="AIC",DevExpl="Deviance Explained",Rsq="Adjusted R\u00B2"),
           Source=factor(Source,levels=c("Primary","ERA5","MERRA2")))%>%
    ggplot(aes(x=Bacteria,y=Value,fill=Source))+
    geom_col(position=position_dodge(0.7),width=0.6)+
    geom_text(data=pf%>%mutate(Source=factor(Source,levels=c("Primary","ERA5","MERRA2")),
                                Metric="AIC",Value=AIC),
              aes(label=paste0("N=",N)),position=position_dodge(0.7),vjust=-0.3,size=2.5)+
    facet_wrap(~Metric,scales="free_y",nrow=1)+
    scale_fill_manual(values=src_colors,labels=src_labels,name="Climate Data Source")+
    labs(title="Level 2: Model Performance (Full Sample)",x=NULL,y=NULL)+
    theme_s+theme(axis.text.x=element_text(angle=45,hjust=1))
  ggsave(file.path(dirs$model,"performance_barplot.pdf"),p_perf,width=14,height=5,dpi=300)
  ggsave(file.path(dirs$model,"performance_barplot.png"),p_perf,width=14,height=5,dpi=300)
}

# EDF heatmaps
for(m in c("full","common")){
  st<-smooth_all%>%filter(Mode==m); if(nrow(st)==0)next
  p_edf<-st%>%
    mutate(
      Source=factor(Source,levels=c("Primary","ERA5","MERRA2")),
      Bacteria=factor(Bacteria, levels=bacteria_order),
      CV=factor(CV, levels=c("Temperature","Humidity","Precipitation","WetDays"),
        labels=c("Temperature","Humidity","Precipitation","Wet Days")),
      lab=ifelse(PenOut,"Out",paste0("edf ",round(EDF,1),"\nF ",ifelse(F_stat>0.1,round(F_stat,1),"<0.1"),SigL)),
      fv=ifelse(PenOut,NA_real_,EDF),
      lab_face=ifelse(PenOut,"italic","plain"),
      lab_col=ifelse(PenOut,"#767676","#1F1F1F"))%>%
    ggplot(aes(x=Source,y=CV,fill=fv))+
    geom_tile(color="white",linewidth=1.1)+
    geom_text(aes(label=lab,fontface=lab_face,color=lab_col),size=2.35,lineheight=0.88,show.legend=FALSE)+
    facet_wrap(~Bacteria,nrow=1)+
    scale_fill_gradientn(
      colours=smooth_fill_palette,
      values=rescale(c(0,0.5,1.5,3,6)),
      na.value="#F1F1F1",limits=c(0,6),name="EDF",oob=squish)+
    scale_color_identity()+
    labs(title=paste0("Level 2: Smooth Terms (",tools::toTitleCase(m)," Sample)"),
         subtitle="Grey tiles indicate penalised-out smooths (EDF < 0.01); text shows edf and F statistic with significance markers",
         x=NULL,y=NULL)+
    guides(fill=guide_colorbar(barheight=unit(3.2,"cm")))+
    theme_s+theme(
      axis.text.x=element_text(size=8,angle=20,hjust=1),
      panel.grid=element_blank(),
      panel.spacing=unit(1.5,"mm"),
      strip.background=element_rect(fill="#EEF3F7",colour=NA),
      legend.position="right")
  ggsave(file.path(dirs$model,paste0("smooth_heatmap_",m,".pdf")),p_edf,width=16,height=5.4,dpi=300,bg="white")
  ggsave(file.path(dirs$model,paste0("smooth_heatmap_",m,".png")),p_edf,width=16,height=5.4,dpi=300,bg="white")
}

# Bridge table
bridge <- perf_all %>% mutate(Tag=paste0(Source,"_",Mode)) %>%
  select(Bacteria,Tag,N,Countries,AIC,DevExpl,Rsq)
write.csv(bridge,file.path(dirs$bridge,"bridge_table.csv"),row.names=FALSE)
draw_bridge_figure_table(
  bridge,
  pdf_path = file.path(dirs$bridge, "bridge_table_figurestyle.pdf"),
  png_path = file.path(dirs$bridge, "bridge_table_figurestyle.png")
)


# Part 4: Curve and percentile comparison
cat("\n",strrep("=",70),"\n PART 4: Curve Comparison\n",strrep("=",70),"\n\n")

curve_stats_list <- list()

for(m in c("full","common")){
  cat("  Mode:",m,"\n")
  for(bt in names(R)){
    pl<-list()
    for(vn in clim_vars){
      cdl<-list()
      for(src in names(R[[bt]])){
        r<-R[[bt]][[src]][[m]]
        if(!is.null(r)&&!is.null(r$curves[[vn]])) cdl[[src]]<-r$curves[[vn]]
      }
      if(length(cdl)<2) next
      cc<-bind_rows(cdl)%>%mutate(Source=factor(Source,levels=c("Primary","ERA5","MERRA2")))

      # Pairwise
      sn<-names(cdl)
      for(i in seq_along(sn))for(j in seq_along(sn)){
        if(j<=i)next; c1<-cdl[[sn[i]]];c2<-cdl[[sn[j]]]
        cx<-seq(max(min(c1$x_orig,na.rm=T),min(c2$x_orig,na.rm=T)),
                min(max(c1$x_orig,na.rm=T),max(c2$x_orig,na.rm=T)),length.out=200)
        o1<-approx(c1$x_orig,c1$or,xout=cx)$y;o2<-approx(c2$x_orig,c2$or,xout=cx)$y
        ok<-!is.na(o1)&!is.na(o2);if(sum(ok)<20)next
        o1v<-o1[ok];o2v<-o2[ok];cxv<-cx[ok]
        s1<-classify_shape(cxv,o1v);s2<-classify_shape(cxv,o2v)
        curve_stats_list[[length(curve_stats_list)+1]]<-tibble(
          Bacteria=bt,Variable=vn2label[vn],S1=sn[i],S2=sn[j],Mode=m,
          Curve_r=cor(o1v,o2v),MaxDiff=max(abs(o1v-o2v)),MeanDiff=mean(abs(o1v-o2v)),
          DirAgree=sign(cor(cxv,o1v))==sign(cor(cxv,o2v)),
          Shape1=s1,Shape2=s2,ShapeAgree=(s1==s2),
          SameSrc=is_rh_same(vn2label[vn],sn[i],sn[j]))
      }

      # ── Per-variable curve plot (NO legend; will add shared one) ──
      p<-ggplot(cc,aes(x=x_orig,y=or,color=Source,linetype=Source,fill=Source))+
        geom_ribbon(aes(ymin=or_lower95,ymax=or_upper95),alpha=src_ribbon_alpha,color=NA)+
        geom_line(linewidth=1.05,lineend="round")+
        geom_hline(yintercept=1,linetype="dashed",color="grey50",linewidth=0.4)+
        scale_color_manual(values=src_colors,labels=src_labels)+
        scale_fill_manual(values=src_colors,labels=src_labels)+
        scale_linetype_manual(values=src_ltys,labels=src_labels)+
        labs(title=vn2label[vn],x=vn2label[vn],y="OR (95% CI)")+
        theme_s+theme(
          legend.position="none",
          plot.title=element_text(size=10.5,margin=margin(b=4)),
          panel.grid.major=element_line(linewidth=0.20,color="#EEEEEE"),
          axis.title.x=element_text(size=8.6),
          axis.title.y=element_text(size=8.6)
        )
      pl[[vn]]<-p
    }

    # ── Per-bacteria combined figure WITH LEGEND ──
    avv<-intersect(c("TMP_scaled_lag","HUM_scaled_lag","PREC_scaled_lag","WET_scaled_lag"),names(pl))
    if(length(avv)>=4){
      pg<-plot_grid(pl[["TMP_scaled_lag"]],pl[["HUM_scaled_lag"]],
                    pl[["PREC_scaled_lag"]],pl[["WET_scaled_lag"]],ncol=2,align="hv")
      legend_grob <- make_legend_grob()
      fp<-plot_grid(pg, legend_grob, ncol=1, rel_heights=c(1, 0.16))
      title <- ggdraw() + draw_label(
        paste0(bt," \u2014 Exposure-Response Curves (",tools::toTitleCase(m)," Sample)"),
        fontface="bold", size=13, x=0.5, hjust=0.5, color="black")
      fp_titled <- plot_grid(title, fp, ncol=1, rel_heights=c(0.04, 1))
      ggsave(file.path(dirs$curve,paste0("curves_",gsub("-","",bt),"_",m,".pdf")),
             fp_titled,width=10,height=10,dpi=300,bg="white")
      ggsave(file.path(dirs$curve,paste0("curves_",gsub("-","",bt),"_",m,".png")),
             fp_titled,width=10,height=10,dpi=300,bg="white")
    }
  }
}

curve_stats_df <- bind_rows(curve_stats_list)
write.csv(curve_stats_df,file.path(dirs$curve,"curve_agreement.csv"),row.names=FALSE)

# Alternative view: curve-correlation lollipop plot
if (nrow(curve_stats_df) > 0) {
  for (m in c("full", "common")) {
    curve_alt <- curve_stats_df %>%
      filter(Mode == m, !SameSrc) %>%
      mutate(
        Pair = factor(mapply(pair_label, S1, S2), levels = pair_levels),
        ShapeFlag = factor(ifelse(ShapeAgree, "Shape-concordant", "Shape-different"),
                           levels = c("Shape-concordant", "Shape-different")),
        Bacteria = factor(Bacteria, levels = bacteria_order),
        Variable = factor(Variable, levels = curve_variable_order)
      )
    if (nrow(curve_alt) == 0) next

    p_curve_alt <- ggplot(curve_alt, aes(x = Bacteria, y = Curve_r, color = Pair, shape = ShapeFlag)) +
      geom_hline(
        data = tibble(yint = c(0.60, 0.80), lty = c("dotted", "dashed"),
                      col = c("#C9C9C9", "#9E9E9E"), lwd = c(0.35, 0.45)),
        aes(yintercept = yint, linetype = lty, color = col, linewidth = lwd),
        inherit.aes = FALSE, show.legend = FALSE
      ) +
      geom_linerange(
        aes(ymin = 0.55, ymax = Curve_r),
        position = position_dodge(width = 0.68),
        linewidth = 0.62, alpha = 0.82
      ) +
      geom_point(
        position = position_dodge(width = 0.68),
        size = 2.8, stroke = 0.8, fill = "white"
      ) +
      coord_flip() +
      facet_wrap(~Variable, nrow = 1) +
      scale_color_manual(values = pair_colors, drop = FALSE) +
      scale_shape_manual(values = c("Shape-concordant" = 21, "Shape-different" = 4), drop = FALSE) +
      scale_y_continuous(limits = c(0.55, 1.00), breaks = c(0.60, 0.80, 1.00)) +
      labs(
        title = paste0("Alternative View: Cross-product Curve Agreement (", tools::toTitleCase(m), " Sample)"),
        subtitle = "Point position shows mean curve correlation against an alternative product; open circles denote shape concordance and crosses denote shape differences",
        x = NULL, y = "Mean curve correlation", color = "Product pair", shape = NULL
      ) +
      theme_s +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, color = "#E9E9E9"),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom"
      )
    ggsave(file.path(dirs$curve, paste0("curve_correlation_lollipop_", m, ".pdf")),
           p_curve_alt, width = 14, height = 5.8, dpi = 300, bg = "white")
    ggsave(file.path(dirs$curve, paste0("curve_correlation_lollipop_", m, ".png")),
           p_curve_alt, width = 14, height = 5.8, dpi = 300, bg = "white")
  }
}

# Percentile heatmaps (full + common)
for (m in c("full", "common")) {
  pct_m <- pct_all %>% filter(Mode == m)
  if (nrow(pct_m) == 0) next

  p_pct <- pct_m %>%
    mutate(
      Source = factor(Source, levels = c("Primary", "ERA5", "MERRA2")),
      lab = sprintf("%.2f", ORratio)
    ) %>%
    ggplot(aes(x = Source, y = Bacteria, fill = DlogOR)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = lab), size = 2.8) +
    facet_wrap(~CV, nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "\u0394log(OR)\nP90-P10", limits = c(-1, 1), oob = squish
    ) +
    labs(
      title = paste0("Level 3: Percentile Effect Contrast (", tools::toTitleCase(m), " Sample)"),
      subtitle = "Cell = OR(P90)/OR(P10); fill = \u0394log(OR)",
      x = NULL, y = NULL
    ) +
    theme_s +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1), panel.grid = element_blank())

  ggsave(file.path(dirs$curve, paste0("percentile_heatmap_", m, ".pdf")), p_pct, width = 14, height = 5, dpi = 300)
  ggsave(file.path(dirs$curve, paste0("percentile_heatmap_", m, ".png")), p_pct, width = 14, height = 5, dpi = 300)

  if (m == "full") {
    ggsave(file.path(dirs$curve, "percentile_heatmap.pdf"), p_pct, width = 14, height = 5, dpi = 300)
    ggsave(file.path(dirs$curve, "percentile_heatmap.png"), p_pct, width = 14, height = 5, dpi = 300)
  }

  # Alternative view: directional percentile-contrast dot plot
  p_pct_alt <- pct_m %>%
    mutate(
      Source = factor(Source, levels = c("Primary", "ERA5", "MERRA2")),
      Bacteria = factor(Bacteria, levels = bacteria_order),
      CV = factor(CV, levels = c("Temperature", "Humidity", "Precipitation", "WetDays"),
                  labels = c("Temperature", "Humidity", "Precipitation", "Wet Days"))
    ) %>%
    ggplot(aes(x = Bacteria, y = DlogOR, color = Source, shape = Source)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45) +
    geom_linerange(
      aes(ymin = 0, ymax = DlogOR),
      position = position_dodge(width = 0.68),
      linewidth = 0.68, alpha = 0.85
    ) +
    geom_point(
      position = position_dodge(width = 0.68),
      size = 2.6, stroke = 0.8, fill = "white"
    ) +
    coord_flip() +
    facet_wrap(~CV, nrow = 1) +
    scale_color_manual(values = src_colors, labels = src_labels) +
    scale_shape_manual(values = src_shape_values, labels = src_labels) +
    labs(
      title = paste0("Alternative View: Percentile-based Effect Contrast (", tools::toTitleCase(m), " Sample)"),
      subtitle = "Segments extend from 0 to \u0394log(OR) between the 10th and 90th percentiles; positive values indicate higher fitted OR at P90 than at P10",
      x = NULL, y = "\u0394log(OR), P90 \u2212 P10", color = "Climate data source", shape = "Climate data source"
    ) +
    theme_s +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "#E9E9E9"),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
  ggsave(file.path(dirs$curve, paste0("percentile_contrast_dotplot_", m, ".pdf")),
         p_pct_alt, width = 14, height = 5.8, dpi = 300, bg = "white")
  ggsave(file.path(dirs$curve, paste0("percentile_contrast_dotplot_", m, ".png")),
         p_pct_alt, width = 14, height = 5.8, dpi = 300, bg = "white")
}


# Part 5: Robustness assessment
cat("\n",strrep("=",70),"\n PART 5: Robustness\n",strrep("=",70),"\n\n")

rob_list<-list(); cvns<-c("Temperature","Precipitation","Humidity","WetDays")
for(bt in unique(smooth_all$Bacteria)) for(cv in cvns) for(m in c("full","common")){
  cls<-classify_rob(smooth_all,curve_stats_df,pct_all,bt,cv,m)
  pr<-smooth_all%>%filter(Bacteria==bt,CV==cv,Source=="Primary",Mode==m)
  er<-smooth_all%>%filter(Bacteria==bt,CV==cv,Source=="ERA5",Mode==m)
  mr<-smooth_all%>%filter(Bacteria==bt,CV==cv,Source=="MERRA2",Mode==m)
  pp<-pct_all%>%filter(Bacteria==bt,CV==cv,Source=="Primary",Mode==m)
  ep<-pct_all%>%filter(Bacteria==bt,CV==cv,Source=="ERA5",Mode==m)
  mp<-pct_all%>%filter(Bacteria==bt,CV==cv,Source=="MERRA2",Mode==m)
  rob_list[[length(rob_list)+1]]<-tibble(
    Bacteria=bt,CV=cv,Mode=m,
    PrimaryF=if(nrow(pr)>0)round(pr$F_stat[1],2)else NA,
    PrimarySig=if(nrow(pr)>0)pr$SigL[1]else NA,
    PrimaryPenOut=if(nrow(pr)>0)pr$PenOut[1]else NA,
    ERA5F=if(nrow(er)>0)round(er$F_stat[1],2)else NA,
    ERA5Sig=if(nrow(er)>0)er$SigL[1]else NA,
    ERA5PenOut=if(nrow(er)>0)er$PenOut[1]else NA,
    MERRA2F=if(nrow(mr)>0)round(mr$F_stat[1],2)else NA,
    MERRA2Sig=if(nrow(mr)>0)mr$SigL[1]else NA,
    MERRA2PenOut=if(nrow(mr)>0)mr$PenOut[1]else NA,
    PrimaryOR=if(nrow(pp)>0)round(pp$ORratio[1],3)else NA,
    ERA5OR=if(nrow(ep)>0)round(ep$ORratio[1],3)else NA,
    MERRA2OR=if(nrow(mp)>0)round(mp$ORratio[1],3)else NA,
    EvalAltSources=cls$Eval_Alt_Sources,
    PrimarySupported=cls$Primary_Supported,
    SigRetainedRate=round(cls$Sig_Retained_Rate, 3),
    Curve_r_Mean=round(cls$Curve_r_Mean, 3),
    ShapeAgreementRate=round(cls$Shape_Agreement_Rate, 3),
    DirectionAgreementRate=round(cls$Direction_Agreement_Rate, 3),
    MagnitudeRatioMean=round(cls$Magnitude_Ratio_Mean, 3),
    ReversalSubstantiveRate=round(cls$Reversal_Substantive_Rate, 3),
    Robustness=cls$Robustness,
    RobustnessReason=cls$Robustness_Reason
  )
}
rob_df<-bind_rows(rob_list)
write.csv(rob_df,file.path(dirs$robust,"robustness_table.csv"),row.names=FALSE)

rcm<-c("Robust"="#1A9850","Robust (minor variation)"="#66BD63",
  "Robust (consistently NS)"="#A6D96A","Attenuated"="#FEE08B",
  "Product-sensitive (lost)"="#F46D43","Product-sensitive (reversal)"="#D73027",
  "Product-sensitive (emergent)"="#F17C67",
  "Insufficient data"="grey80","No primary"="grey90")

for(m in c("full","common")){
  rd<-rob_df%>%filter(Mode==m); if(nrow(rd)==0)next
  p_rob<-rd%>%ggplot(aes(x=CV,y=Bacteria,fill=Robustness))+
    geom_tile(color="white",linewidth=1.2)+
    geom_text(aes(label=paste0(
      ifelse(is.na(PrimarySig),"-",PrimarySig),"/",
      ifelse(is.na(ERA5Sig),"-",ERA5Sig),"/",
      ifelse(is.na(MERRA2Sig),"-",MERRA2Sig))),size=2.8)+
    scale_fill_manual(values=rcm,name="Robustness",drop=FALSE)+
    labs(title=paste0("Level 4: Robustness (",tools::toTitleCase(m)," Sample)"),
      subtitle="Cell: significance (Primary / ERA5 / MERRA-2)",x=NULL,y=NULL)+
    theme_s+theme(panel.grid=element_blank(),legend.position="right",
      legend.key.size=unit(0.5,"cm"),legend.text=element_text(size=8))
  ggsave(file.path(dirs$robust,paste0("robustness_",m,".pdf")),p_rob,width=13,height=6,dpi=300)
  ggsave(file.path(dirs$robust,paste0("robustness_",m,".png")),p_rob,width=13,height=6,dpi=300)

  draw_robustness_figure_table(
    rob_df,
    pdf_path = file.path(dirs$robust, paste0("robustness_table_", m, "_figurestyle.pdf")),
    png_path = file.path(dirs$robust, paste0("robustness_table_", m, "_figurestyle.png")),
    mode_label = m
  )
}

cat("\n  FULL:\n"); rf<-rob_df%>%filter(Mode=="full");rc<-table(rf$Robustness)
for(c in names(rc)) cat(sprintf("    %s: %d\n",c,rc[c]))
cat("  COMMON:\n"); rf2<-rob_df%>%filter(Mode=="common");rc2<-table(rf2$Robustness)
for(c in names(rc2)) cat(sprintf("    %s: %d\n",c,rc2[c]))


# Part 6: Combined curve figures
cat("\n",strrep("=",70),"\n PART 6: Combined Figures\n",strrep("=",70),"\n\n")

bo<-c("3GCR-Ec","3GCR-Kp","CR-Ab","CR-Ec","CR-Kp","CR-Pa")
vo<-c("TMP_scaled_lag","HUM_scaled_lag","PREC_scaled_lag","WET_scaled_lag")
vt<-c("Temperature","Humidity","Precipitation","Wet Days")

for(m in c("full","common")){
  acp<-list()
  for(bt in names(R)) for(vn in clim_vars){
    cdl<-list()
    for(src in names(R[[bt]])){
      r<-R[[bt]][[src]][[m]]
      if(!is.null(r)&&!is.null(r$curves[[vn]])) cdl[[src]]<-r$curves[[vn]]
    }
    if(length(cdl)==0)next
    cc<-bind_rows(cdl)%>%mutate(Source=factor(Source,levels=c("Primary","ERA5","MERRA2")))
    p<-ggplot(cc,aes(x=x_orig,y=or,color=Source,linetype=Source))+
      geom_ribbon(aes(ymin=or_lower95,ymax=or_upper95,fill=Source),alpha=src_ribbon_alpha,color=NA)+
      geom_line(linewidth=0.78,lineend="round")+
      geom_hline(yintercept=1,linetype="dashed",color="grey40",linewidth=0.4)+
      scale_color_manual(values=src_colors)+scale_fill_manual(values=src_colors)+
      scale_linetype_manual(values=src_ltys)+
      labs(x=NULL,y=NULL)+theme_bw(base_size=8.8)+
      theme(
        legend.position="none",
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(linewidth=0.18,color="#EEEEEE"),
        axis.text=element_text(size=8.8, colour="#343434"),
        axis.title.y=element_text(size=10.6, face="bold", margin=margin(r=5)),
        plot.margin=margin(3,5,3,5))
    acp[[paste(bt,vn,sep="_")]]<-p
  }
  if(length(acp)==0)next

  ordered_keys <- unlist(lapply(bo, function(bt) paste(bt, vo, sep="_")), use.names = FALSE)
  ordered_plots <- lapply(seq_along(ordered_keys), function(i) {
    k <- ordered_keys[[i]]
    p <- if (k %in% names(acp)) acp[[k]] else ggplot() + theme_void()
    if (((i - 1) %% length(vo)) == 0) p <- p + labs(y = "OR")
    p
  })
  aligned_plots <- align_plots(plotlist = ordered_plots, align = "hv", axis = "tblr")

  # Title row aligned to the plotting region
  tr<-plot_grid(
    ggplot()+theme_void(),
    plot_grid(
      plotlist=lapply(vt,function(t)
        ggplot()+theme_void()+ggtitle(t)+
          theme(plot.title=element_text(hjust=0.5,size=11.8,face="bold"))),
      ncol=4, rel_widths=rep(1,4)),
    ncol=2, rel_widths=c(0.08,0.92))

  # Bacteria rows
  rps<-list(tr)
  for(i in seq_along(bo)){
    bt <- bo[[i]]
    row <- aligned_plots[((i - 1) * length(vo) + 1):(i * length(vo))]
    lr<-plot_grid(
      ggplot()+theme_void()+
        annotate("text",x=.5,y=.5,label=bt,fontface="bold",size=3.8)+
        coord_cartesian(clip="off"),
      plot_grid(plotlist=row,ncol=4,align="hv",axis="tblr",rel_widths=rep(1,4)),
      ncol=2, rel_widths=c(0.08,0.92))
    rps[[length(rps)+1]]<-lr
  }

  # ── LEGEND: generous height allocation ──
  legend_grob <- make_legend_grob()
  rps[[length(rps)+1]] <- legend_grob

  # Assemble with custom ribbon-aware legend space
  mf <- plot_grid(plotlist=rps, ncol=1,
                  rel_heights=c(0.08, rep(1, length(bo)), 0.18))

  # Add overall title
  title_text <- paste0("Exposure\u2013response curves across climate products (",
                        tools::toTitleCase(m)," Sample)")
  mf_titled <- plot_grid(
    ggdraw()+draw_label(title_text,fontface="bold",size=14,x=0.5,hjust=0.5,color="black"),
    mf, ncol=1, rel_heights=c(0.03, 1))

  ggsave(file.path(dirs$robust,paste0("climate_product_curves_",m,".pdf")),
         mf_titled, width=16, height=20, dpi=300, bg="white")
  ggsave(file.path(dirs$robust,paste0("climate_product_curves_",m,".png")),
         mf_titled, width=16, height=20, dpi=300, bg="white")
  cat("  Curve figure saved:",m,"\n")
}


# Final analysis log
log_file <- file.path(output_root, "analysis_log.txt")
log_con <- file(log_file, open="wt")
sink(log_con)
on.exit({ sink(); close(log_con) }, add=TRUE)

cat("ANALYSIS OF CLIMATE-AMR ASSOCIATIONS UNDER ALTERNATIVE REANALYSIS-PRODUCT CONFIGURATIONS\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("DESIGN:\n")
cat("  Track A (Full): same sample as main Model C\n")
cat("  Track B (Common): data-driven row filtering (all 11 climate cols valid)\n")
cat("  Pipeline: identical to main Model C\n\n")
cat("SOURCES:\n  Primary = CRU TS TMP/PREC/WET + ERA5 HUM\n")
cat("  ERA5 = ERA5 TMP/PREC/WET + ERA5 HUM (note: HUM same as Primary)\n")
cat("  MERRA-2 = all four from MERRA-2\n\n")
cat("LAGS:\n")
for(bt in names(lag_settings)){ls<-lag_settings[[bt]]
  cat(sprintf("  %s: TMP=%d PREC=%d HUM=%d WET=%d\n",bt,ls$temp_lag,ls$precip_lag,ls$humid_lag,ls$wetdays_lag))}
cat("\nBRIDGE TABLE:\n"); if(exists("bridge")) print(as.data.frame(bridge))
cat("\nROBUSTNESS (Full):\n")
if(exists("rob_df")) print(as.data.frame(rob_df%>%filter(Mode=="full")%>%select(Bacteria,CV,Robustness)))
cat("\nROBUSTNESS (Common):\n")
if(exists("rob_df")) print(as.data.frame(rob_df%>%filter(Mode=="common")%>%select(Bacteria,CV,Robustness)))
if(length(error_log)>0){cat("\nERRORS:\n");print(as.data.frame(bind_rows(error_log)))}
cat("\n"); print(sessionInfo())

sink(); close(log_con)

cat("\n",strrep("=",60),"\n  ANALYSIS COMPLETE\n  ",output_root,"\n",strrep("=",60),"\n")
