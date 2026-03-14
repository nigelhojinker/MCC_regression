###########################################
#### CellChat and CellPhoneDB analysis ####
###########################################
run_cellchatv2 <- function(obj, mean, database, threshold){
  cellchat    <- createCellChat(object = obj, group.by = "labels")

  # Use ORIGINAL databases
  if (database == "CCDB") {
    cellchat@DB <- CellChatDB.human
  } else if(database == "CPDB") {
    cellchat@DB <- CellPhoneDB_ori.human
  }

  cellchat <- subsetData(cellchat) %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = mean, trim = threshold) %>%
    filterCommunication(min.cells = 10)

  return(cellchat)
}

run_cellphonedbv2 <- function(obj, mean, database, threshold) {
  name <- deparse(substitute(obj))
  ad <- import("anndata")
  cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

  # Load ORIGINAL database
  if (database == "CPDB") {
    cpdbPath <- "data/cellphonedb.zip"
  } else if(database == "CCDB") {
    cpdbPath <- "data/cellchat_ori.zip"
  }

  prefix <- file.path(tempdir(), name)

  dir.create(prefix %>% normalizePath(winslash = "/"), recursive = TRUE)

  meta <- obj@meta.data %>%
    rownames_to_column(var = "barcode") %>%
    select(barcode, labels) %>%
    write_delim(file = file.path(prefix, "input_meta.tsv"), delim = "\t")

  normcounts <- obj[["RNA"]]$data
  genes <- data.frame(gene = rownames(normcounts)) %>%
    column_to_rownames("gene")
  barcode <- data.frame(barcode = colnames(normcounts)) %>%
    column_to_rownames("barcode")

  adata <- ad$AnnData(X = normcounts, obs = genes, var = barcode)
  adata <- adata$transpose()
  adata$write_h5ad(file.path(prefix, "input.h5ad"))

  cpdb <- cpdb.analysis$call(cpdb_file_path = cpdbPath,
                             meta_file_path = file.path(prefix, "input_meta.tsv"),
                             counts_file_path = file.path(prefix, "input.h5ad"),
                             counts_data = "hgnc_symbol",
                             output_path = prefix,
                             avg_method = mean,
                             threshold = threshold,
                             threads = 1L,
                             debug_seed = 42L)

  class(cpdb) <- "CellPhoneDB"

  return(cpdb)
}

#################################################
#### Compare CellChat and CellPhoneDB output ####
#################################################
crosscheckv2 <- function(obj1 = NULL, obj2 = NULL,
                         comp.type = c("CC_CP", "CP_CC", "CC_CC", "CP_CP"),
                         threshold = 0.05, return.all = FALSE) {

  if(is.null(obj1)) stop("obj1 must be supplied.")
  if(is.null(obj2)) stop("obj2 must be supplied.")

  comp.type <- match.arg(comp.type)

  obj1.name <- deparse(substitute(obj1))
  obj2.name <- deparse(substitute(obj2))

  if(comp.type == "CC_CP") {
    if(! inherits(obj1, "CellChat"))    stop(paste(obj1.name, "is not a CellChat object."))
    if(! inherits(obj2, "CellPhoneDB")) stop(paste(obj2.name, "is not a CellPhoneDB object."))

    return(crosscheck_CC_CPDBv2(obj1, obj2, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CP_CC") {
    if(! inherits(obj1, "CellPhoneDB")) stop(paste(obj1.name, "is not a CellPhoneDB object."))
    if(! inherits(obj2, "CellChat"))    stop(paste(obj2.name, "is not a CellChat object."))

    return(crosscheck_CC_CPDBv2(obj2, obj1, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CC_CC") {
    if(! inherits(obj1, "CellChat"))    stop(paste(obj1.name, "is not a CellChat object."))
    if(! inherits(obj2, "CellChat"))    stop(paste(obj2.name, "is not a CellChat object."))

    return(crosscheck_CC_CCv2(obj1, obj2, name1 = obj1.name, name2 = obj2.name, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CP_CP") {
    if(! inherits(obj1, "CellPhoneDB")) stop(paste(obj1.name, "is not a CellPhoneDB object."))
    if(! inherits(obj2, "CellPhoneDB")) stop(paste(obj2.name, "is not a CellPhoneDB object."))

    return(crosscheck_CPDB_CPDBv2(obj1, obj2, name1 = obj1.name, name2 = obj2.name, threshold = threshold, return.all = return.all))
  }
}

crosscheck_CC_CPDBv2 <- function(cellchat, cellphonedb, threshold = 0.05, return.all = FALSE) {

  cellchat_res <- pull_netslot(cellchat) %>%
    left_join(step1_map %>% select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  cellphonedb_res <- cellphonedb$pvalues %>%
    pivot_longer(cols = contains("|"),
                 names_to = "cell_type_pair",
                 values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    mutate(partner_a = gsub(".*:", "", partner_a),
           partner_b = gsub(".*:", "", partner_b)) %>%
    left_join(step1_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  levels <- c("Sig", "Not_Sig", "Not_Found")

  combine.all <- cellchat_res %>% full_join(cellphonedb_res %>%
                                              select(id_cp_interaction, interacting_pair, source, target,  unique_int, pval),
                                            by = "unique_int",
                                            suffix = c("_cc", "_cpdb")) %>%
    distinct() %>%
    mutate(category = case_when(pval_cc < threshold & pval_cpdb < threshold ~ "Sig_Both",
                                pval_cc < threshold & is.na(pval_cpdb) ~ "Sig_CellChat_Not_Found_CellPhoneDB",
                                is.na(pval_cc) & pval_cpdb < threshold ~ "Sig_CellPhoneDB_Not_Found_CellChat",
                                pval_cc < threshold & pval_cpdb >= threshold ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                pval_cc >= threshold & pval_cpdb < threshold ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
           CellPhoneDB_status = case_when(pval_cpdb < threshold ~ "Sig",
                                          pval_cpdb >= threshold ~ "Not_Sig", is.na(pval_cpdb) ~ "Not_Found"),
           CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
           CellChat_status = case_when(pval_cc < threshold ~ "Sig",
                                       pval_cc >= threshold ~ "Not_Sig",
                                       is.na(pval_cc) ~ "Not_Found"),
           CellChat_status = factor(CellChat_status, levels = levels))

  summary <- table(CellChat = combine.all$CellChat_status,
                   CellPhoneDB = combine.all$CellPhoneDB_status) %>% addmargins()
  if (!return.all) {
    res <- combine.all %>% filter(!is.na(category))
  }
  cat("Comparative analysis done successfully for p-value",
      threshold, "\n", "scpeakeR summary: \n")
  print(summary)
  return(list(result = res, summary = summary))

}

crosscheck_CC_CCv2 <- function(cellchat1, cellchat2, name1, name2, threshold = 0.05, return.all = FALSE) {

  prob1 <- paste0("prob_", name1)
  pval1 <- paste0("pval_", name1)
  prob2 <- paste0("prob_", name2)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  levels <- c("Sig", "Not_Sig", "Not_Found")

  ## Process first CellChat object
  cellchat_res1 <- pull_netslot(cellchat1) %>%
    left_join(step1_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob1 := prob, !!pval1 := pval)

  ## Process second CellChat object
  cellchat_res2 <- pull_netslot(cellchat2) %>%
    left_join(step1_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob2 := prob, !!pval2 := pval)

  ## Combining results
  res <- full_join(cellchat_res1, cellchat_res2, by = "unique_int", suffix = c("_cc1", "_cc2")) %>%
    distinct() %>%
    mutate(interaction_name = coalesce(interaction_name_cc1, interaction_name_cc2),
           pathway_name     = coalesce(pathway_name_cc1, pathway_name_cc2),
           source     = coalesce(source_cc1,   source_cc2),
           target     = coalesce(target_cc1,   target_cc2),
           ligand     = coalesce(ligand_cc1,   ligand_cc2),
           receptor   = coalesce(receptor_cc1, receptor_cc2),
           annotation = coalesce(annotation_cc1, annotation_cc2),
           evidence   = coalesce(evidence_cc1, evidence_cc2)) %>%
    select(interaction_name, pathway_name, source, target, ligand, receptor,
           unique_int, starts_with("prob"), starts_with("pval"), annotation, evidence) %>%
    mutate(category = case_when(!!sym(pval1) < threshold  & !!sym(pval2) < threshold  ~ "Sig_Both",
                                !!sym(pval1) < threshold  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                is.na(!!sym(pval1))  & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                !!sym(pval1) < threshold  & !!sym(pval2) >= threshold ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                !!sym(pval1) >= threshold & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
           !!status1 := factor(
             case_when(!!sym(pval1) < threshold ~ "Sig",
                       !!sym(pval1) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval1)) ~ "Not_Found"),
             levels = levels),
           !!status2 := factor(
             case_when(!!sym(pval2) < threshold ~ "Sig",
                       !!sym(pval2) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval2)) ~ "Not_Found"),
             levels = levels))

  summary <- table(res[[status1]], res[[status2]]) %>%
    addmargins()

  dimnames(summary) <- dimnames(summary) %>%
    setNames(c(name1, name2))

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "scpeakeR summary: \n")

  print(summary)

  return(list(result = res, summary = summary))
}

crosscheck_CPDB_CPDBv2 <- function(cellphonedb1, cellphonedb2, name1, name2, threshold = 0.05, return.all = FALSE) {

  pval1 <- paste0("pval_", name1)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  levels <- c("Sig", "Not_Sig", "Not_Found")

  ## CellPhoneDB object 1 processing
  cellphonedb_res1 <- cellphonedb1$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    mutate(partner_a = gsub(".*:", "", partner_a),
           partner_b = gsub(".*:", "", partner_b)) %>%
    left_join(step1_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval1 := pval)

  ## CellPhoneDB object 2 processing
  cellphonedb_res2 <- cellphonedb2$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    mutate(partner_a = gsub(".*:", "", partner_a),
           partner_b = gsub(".*:", "", partner_b)) %>%
    left_join(step1_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval2 := pval)

  ## Combining results
  res <- full_join(cellphonedb_res1, cellphonedb_res2, by = "unique_int", suffix = c("_cp1", "_cp2")) %>%
    distinct() %>%
    mutate(id_cp_interaction = coalesce(id_cp_interaction_cp1, id_cp_interaction_cp2),
           interacting_pair  = coalesce(interacting_pair_cp1, interacting_pair_cp2),
           source     = coalesce(source_cp1,   source_cp2),
           target     = coalesce(target_cp1,   target_cp2),
           partner_a  = coalesce(partner_a_cp1, partner_a_cp2),
           partner_b  = coalesce(partner_b_cp1, partner_b_cp2),
           gene_a = coalesce(gene_a_cp1, gene_b_cp2),
           gene_b = coalesce(gene_b_cp1, gene_b_cp2),
           secreted = coalesce(secreted_cp1, secreted_cp2),
           receptor_a = coalesce(receptor_a_cp1, receptor_b_cp2),
           receptor_b = coalesce(receptor_b_cp1, receptor_b_cp2),
           annotation_strategy = coalesce(annotation_strategy_cp1, annotation_strategy_cp2),
           is_integrin = coalesce(is_integrin_cp1, is_integrin_cp2),
           directionality = coalesce(directionality_cp1, directionality_cp2),
           classification = coalesce(classification_cp1, classification_cp2)) %>%
    select(id_cp_interaction, interacting_pair, source, target, partner_a, partner_b,
           gene_a, gene_b, starts_with("pval"), secreted, receptor_a, receptor_b,
           annotation_strategy, is_integrin, directionality, classification, unique_int) %>%
    mutate(category = case_when(!!sym(pval1) < threshold  & !!sym(pval2) < threshold  ~ "Sig_Both",
                                !!sym(pval1) < threshold  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                is.na(!!sym(pval1))  & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                !!sym(pval1) < threshold  & !!sym(pval2) >= threshold ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                !!sym(pval1) >= threshold & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
           !!status1 := factor(
             case_when(!!sym(pval1) < threshold ~ "Sig",
                       !!sym(pval1) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval1)) ~ "Not_Found"),
             levels = levels),
           !!status2 := factor(
             case_when(!!sym(pval2) < threshold ~ "Sig",
                       !!sym(pval2) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval2)) ~ "Not_Found"),
             levels = levels))

  summary <- table(res[[status1]], res[[status2]]) %>%
    addmargins()

  dimnames(summary) <- dimnames(summary) %>%
    setNames(c(name1, name2))

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "scpeakeR summary: \n")

  print(summary)

  return(list(result = res, summary = summary))
}

#############################
#### Compute correlation ####
#############################
parse_comptype <- function(obj1, obj2) {
  if (class(obj1) == "CellChat" & class(obj2) == "CellPhoneDB") {
    return("CC_CP")
  } else if (class(obj1) == "CellPhoneDB" & class(obj2) == "CellChat") {
    return("CP_CC")
  } else if (class(obj1) == "CellChat" & class(obj2) == "CellChat") {
    return("CC_CC")
  } else if (class(obj1) == "CellPhoneDB" & class(obj2) == "CellPhoneDB") {
    return("CP_CP")
  }
}

compute_corr <- function(x) {
  step1 <- expand.grid(obj1 = names(x), obj2 = names(x))

  # Fix to allow left_join of columns with incompatible types
  objs <- x[grepl("cellphone.*_CCDB$", names(x))]
  objs <- purrr::map(objs, function(x){
    x$pvalues$directionality <- as.character(x$pvalues$directionality)
    x$pvalues$classification <- as.character(x$pvalues$classification)
    x
  })

  x[names(objs)] <- objs

  cor_list <- apply(step1, 1, function(row) {
    name1 <- row["obj1"]
    name2 <- row["obj2"]

    obj1 <- x[[name1]]
    obj2 <- x[[name2]]

    comp.type <- parse_comptype(obj1, obj2)

    out <- crosscheckv2(obj1, obj2, comp.type)

    res <- out$summary %>%
      as.data.frame()

    tp <- res$Freq[1]
    fp <- res$Freq[2] + res$Freq[3]
    tn <- res$Freq[6] + res$Freq[7] + res$Freq[10]
    fn <- res$Freq[5] + res$Freq[9]

    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

    list(mcc = mcc)

  })

  df <- bind_rows(cor_list)

  step1_cor <- bind_cols(step1, df)
}

#####################
#### MCC Heatmap ####
#####################
create_heatmap <- function(res_df, obj_list, colname) {
  # Create matrix with dimnames
  jac_mat <- matrix(res_df[[colname]],
                    nrow = length(obj_list), ncol = length(obj_list),
                    byrow = FALSE,
                    dimnames = list(names(obj_list), names(obj_list)))

  palette <- colorblind_pal()(8)

  database_order <- factor(
    c(rep("CPDB", 6),
      rep("CCDB", 6)),
    levels = c("CPDB", "CCDB"))

  mean_order <- factor(
    rep(c(
      rep("thresholdedMean_10", 2),
      rep("thresholdedMean_25", 2),
      rep("triMean", 2)),
      2),
    levels = c("thresholdedMean_10", "thresholdedMean_25", "triMean"))

  method_order <- factor(
    rep(c("cellphonedb", "cellchat"), 6),
    levels = c("cellphonedb", "cellchat"))

  db_annotation_row <- HeatmapAnnotation(Method    = method_order,
                                         Mean_Type = mean_order,
                                         Database  = database_order,
                                         col = list(Method = c("cellphonedb" = palette[2], "cellchat" = palette[3]),
                                                    Mean_Type = c("triMean" = palette[5], "thresholdedMean_10" = palette[6], "thresholdedMean_25" = palette[7]),
                                                    Database = c("CPDB" = palette[2], "CCDB" = palette[3])),
                                         which = "row", show_annotation_name = FALSE)

  db_annotation_col <- HeatmapAnnotation(Database  = database_order,
                                         Mean_Type = mean_order,
                                         Method    = method_order,
                                         col = list(Database = c("CPDB" = palette[2], "CCDB" = palette[3]),
                                                    Mean_Type = c("triMean" = palette[5], "thresholdedMean_10" = palette[6], "thresholdedMean_25" = palette[7]),
                                                    Method = c("cellphonedb" = palette[2], "cellchat" = palette[3])),
                                         which = "column", show_legend = FALSE, show_annotation_name = FALSE)

  if (grepl("jaccard", colname)) {
    plot_title = "Jaccard score"
  } else {
    plot_title = "Matthews Correlation\nCoefficient (MCC)"
  }

  ComplexHeatmap::Heatmap(jac_mat,
                          column_names_gp = grid::gpar(fontsize = 12),
                          column_names_rot = 45,
                          row_names_gp = grid::gpar(fontsize = 12),
                          row_names_side = "left",
                          row_names_max_width = unit(20, "cm"),
                          row_title = NULL, column_title = NULL,
                          heatmap_legend_param = list(title = plot_title, fontsize = 10),
                          col = colorRamp2(c(0, 1), c("white", "red")),
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          bottom_annotation = db_annotation_col,
                          left_annotation = db_annotation_row,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", jac_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                          })
}

########################################
#### Linear regression data reshape ####
########################################
format_data <- function(df) {
  obj <- df %>%
    rowwise() %>%
    mutate(objs = paste(sort(c(obj1, obj2)), collapse = "|")) %>%
    ungroup() %>%
    filter(obj1 != obj2) %>%
    distinct(objs, .keep_all = TRUE) %>%
    select(-objs) %>%
    separate(obj1, into = c("method1", "mean1", "database1"), sep = "_") %>%
    separate(obj2, into = c("method2", "mean2", "database2"), sep = "_") %>%
    mutate(Method    = ifelse(method1 == method2,     "Same", "Diff"),
           Mean_Type = ifelse(mean1 == mean2,         "Same", "Diff"),
           Database  = ifelse(database1 == database2, "Same", "Diff"),
           across(c(Method, Mean_Type, Database), ~ factor(.x, levels = c("Same", "Diff")))
    )
}

