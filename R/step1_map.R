pacman::p_load(scpeakeR, ggthemes, ComplexHeatmap, circlize)
setwd(this.path::here())
rm(list = ls())

load("../data/CellPhoneDB_ori_human.rda") # CellChat-compatible CellPhoneDB database in rda format
data(CellChatDB.human)

# Perform Step 1 of Database Harmonization --------------------------------
map <- CellChatDB.human$geneInfo %>%
  select(Symbol, EntryID.uniprot) %>%
  drop_na() %>%
  deframe()

complex <- CellChatDB.human$complex

# Perfrom gene-to-uniprot conversion
complex$subunit_1 <- map[complex$subunit_1]
complex$subunit_2 <- map[complex$subunit_2]
complex$subunit_3 <- map[complex$subunit_3]
complex$subunit_4 <- map[complex$subunit_4]
complex$subunit_5 <- map[complex$subunit_5]

complex <- apply(complex, 1, function(x) paste(na.omit(x), collapse = "_"))

map <- c(map, complex)

# Load Original databases -------------------------------------------------

# CellChatDB
cellchatdb <- CellChatDB.human$interaction

cellchatdb$ligand   <- map[cellchatdb$ligand]
cellchatdb$receptor <- map[cellchatdb$receptor]

cellchatdb %>% filter(is.na(ligand))   %>% pull(interaction_name)
cellchatdb %>% filter(is.na(receptor)) %>% pull(interaction_name)

cellchatdb$LR <- paste(cellchatdb$ligand, cellchatdb$receptor, sep = "|")

unzip(zipfile = "../data/cellchat_ori.zip", files = c("interaction_table.csv", "multidata_table.csv"))

multidata <- read_csv("multidata_table.csv")
cellchatdb_cp <- read_csv("interaction_table.csv") %>%
  left_join(multidata %>%
              select(multidata_1_id = id_multidata, name)) %>%
  rename(ligand = name) %>%
  left_join(multidata %>%
              select(multidata_2_id = id_multidata, name)) %>%
  rename(receptor = name)

# match the NA interaction partners in CellChatDB
cellchatdb_cp <- cellchatdb_cp %>%
  mutate(ligand   = ifelse(ligand   == "P16619", NA, ligand),
         receptor = ifelse(receptor == "P46091", NA, receptor))

cellchatdb_cp$complex_ligand   <- cellchatdb_cp$ligand   %in% names(complex)
cellchatdb_cp$complex_receptor <- cellchatdb_cp$receptor %in% names(complex)

cellchatdb_cp <- cellchatdb_cp %>%
  mutate(ligand   = ifelse(complex_ligand,   complex[cellchatdb_cp$ligand],   ligand),
         receptor = ifelse(complex_receptor, complex[cellchatdb_cp$receptor], receptor),
         LR = paste(ligand, receptor, sep = "|"))

identical(sort(cellchatdb$LR), sort(cellchatdb_cp$LR))

# CellChatDB mapping table
ccdb_map <- full_join(cellchatdb    %>% select(interaction_name,  LR),
                      cellchatdb_cp %>% select(id_cp_interaction, LR)) %>%
  distinct(interaction_name, .keep_all = TRUE)

# CellPhoneDB
cellphonedb_cc <- CellPhoneDB_ori.human$interaction

map <- CellPhoneDB_ori.human$geneInfo %>%
  select(Symbol, EntryID.uniprot) %>%
  drop_na() %>%
  deframe()

cellphonedb_cc$complex_ligand   <- cellphonedb_cc$ligand   %in% rownames(CellPhoneDB_ori.human$complex)
cellphonedb_cc$complex_receptor <- cellphonedb_cc$receptor %in% rownames(CellPhoneDB_ori.human$complex)

cellphonedb_cc <- cellphonedb_cc %>%
  mutate(ligand   = ifelse(! complex_ligand,   map[cellphonedb_cc$ligand],   ligand),
         receptor = ifelse(! complex_receptor, map[cellphonedb_cc$receptor], receptor),
         LR = paste(ligand, receptor, sep = "|"))

multidata   <- read_csv("../data/cellphonedb/multidata_table.csv")
cellphonedb <- read_csv("../data/cellphonedb/interaction_table.csv") %>%
  left_join(multidata %>%
              select(multidata_1_id = id_multidata, name)) %>%
  rename(ligand = name) %>%
  left_join(multidata %>%
              select(multidata_2_id = id_multidata, name)) %>%
  rename(receptor = name) %>%
  mutate(LR = paste(ligand, receptor, sep = "|"))

identical(sort(cellphonedb$LR), sort(cellphonedb_cc$LR))
# [1] FALSE

cellphonedb <- cellphonedb %>%
  mutate(receptor = ifelse(grepl("^GALR", receptor), paste0(receptor, "_GPR151"), receptor),
         LR = paste(ligand, receptor, sep = "|"))

cellphonedb_cc <- cellphonedb_cc %>%
  mutate(ligand = case_when(
    ligand == "P04439" ~ "HLAA",
    ligand == "P01889" ~ "HLAB",
    ligand == "P10321" ~ "HLAC",
    ligand == "P01562" ~ "IFNA1",
    ligand == "P58400" ~ "Q9ULB1",
    .default = ligand),
    receptor = ifelse(receptor == "P58400", "Q9ULB1", receptor),
    LR = paste(ligand, receptor, sep = "|"))

identical(sort(cellphonedb$LR), sort(cellphonedb_cc$LR))

cpdb_map <- full_join(cellphonedb_cc %>% select(interaction_name, LR),
                      cellphonedb %>% select(id_cp_interaction, LR))

ggvenn::ggvenn(list(CellChat = ccdb_map$LR, CellPhoneDB = cpdb_map$LR))

step1_map <- bind_rows(ccdb_map, cpdb_map) %>%
  distinct()

rm(ccdb_map, cpdb_map, cellchatdb, cellchatdb_cp, cellphonedb, cellphonedb_cc, multidata, complex, map)

save(step1_map, file = "../data/step1_map.rda")

