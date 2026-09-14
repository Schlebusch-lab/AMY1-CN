# Genetic relatedness

```{ruby}
# PLINK file set comprised of .bed, .bim and .fam files
DB=ddPCR_data_set

# Run KING
/king -b ${DB}.bed --kinship --prefix KING/${DB} > KING/${DB}_king_info.txt

# Sort individuals based on their kinship coefficient
sort -k9 -r -h KING/${DB}.kin > KING/${DB}_kinship_list.txt

REL_N=20 # Ex. amount of related pairs (20)

# Identify (unique) related individuals
head -n ${REL_N} KING/${DB}_kinship_list.txt | cut -f3 | sort | uniq > related_list.txt

# Get their full names from the .fam file
cat related_list.txt | while read line; \
    do grep ${line}' ' ${DB}.fam >> related_list_full_names.txt; \
    done

# Remove them from dataset
plink --bfile ${DB} --remove related_list_full_names.txt --make-bed --out ${DB}_unrelated
```

# ILR transformation

```{ruby }
# Select compositional variables
comp_vars <- ddPCR_data_set[,c("agriculture", "pastoralism", "fishing", 
                               "gathering", "hunting")]

# Replace zeros if present
comp_vars_no_zeros <- if (any(comp_vars == 0)) {
  cmultRepl(comp_vars, method = "CZM", output = "p-counts")} else {comp_vars}

# Convert to compositional class
comp <- acomp(comp_vars_no_zeros)

# Number of variables
D <- ncol(comp)

# FUNCTION: Create ILR contrast vectors and compute ILR coordinates (new variables)
create_ilr_vector <- function(i, D) {
  v <- rep(-sqrt(1 / (D * (D - 1))), D)
  v[i] <- sqrt((D - 1) / D)
  return(v)}

log_comp <- log(comp)                                                        # Log of variables
log_comp_mat <- matrix(unclass(log_comp), nrow = nrow(log_comp), ncol = D)   # Matrix format
colnames(log_comp_mat) <- colnames(comp)                                     # Column names

ilr_coords <- matrix(NA, nrow = nrow(ddPCR_data_set), ncol = D)              # Empty matrix
colnames(ilr_coords) <- paste0(colnames(comp_vars), "_ILR")                  # Column names

# Apply function to calculate IRL transformations
for (i in 1:D) {
  v <- create_ilr_vector(i, D)
  ilr_coords[, i] <- log_comp_mat %*% v}

# Append ILR coordinates to original dataset
data <- cbind(ddPCR_data_set, ilr_coords)
```


# Models formulation

Below are two examples of the `brms` and `glmmTMB` syntax used in modeling *AMY1* CN in the **ddPCR data set**.  

```{ruby }
# Example of glmmTMB model synthax
model_null <- glmmTMB(AMY1_CN ~ PC1 + PC2 + PC3 + PC4 + (1|population),
                      data=ddPCR_data_set,
                      family=gaussian)

# Example of brms model synthax
model_null <- brm(AMY1_CN ~ (1|gr(sample,cov=VCV)) + (1|population),
                  data = brms_data,
                  data2=list(VCV=VCV),
                  prior=set_prior("normal(0,0.5)", class="sd", group="sample"),
                  family=student(),
                  save_pars=save_pars(all=TRUE),
                  sample_prior=TRUE,
                  chains=4, iter=4000, warmup=3000,
                  control=list(adapt_delta=0.999, max_treedepth=12))
```

# Model assumptions and stability

## `glmmTMB` specifics

```{ruby }
# Models assumptions and stability (glmmTMB)
for (name in model_names) {
  # 'models' is a list of all models run with the same data set
  # (e.g. model_null, model_agriculture, model_pastoralism)
  model <- models[[name]]
  
  # Coefficients summary
  summary(model)$coefficients
  
  # Multicollinearity check
  check_collinearity(model)

  # Overdispersion test
  sim <- simulateResiduals(model)
  testDispersion(sim)

  # DHARMa diagnostic plots
  plot(sim)
}
```

## `brms` specifics

```{ruby }
# Models assumptions and stability (brms)
for (name in model_names) {
  model <- models[[name]]
  
  # Summary, LOO and Pareto K
  summary(model)
  loo(model)

  # Sampler diagnostics
  nuts_params <- rstan::get_sampler_params(model$fit, inc_warmup = FALSE)
  nuts_df <- do.call(rbind, nuts_params)

  # Divergences
  if ("divergent__" %in% colnames(nuts_df)) {
    cat("Divergences:", sum(nuts_df[, "divergent__"]), "\n")}

  # Treedepth
  if ("treedepth__" %in% colnames(nuts_df)) {
    max_td <- attr(model$fit, "stan_args")[[1]]$control$max_treedepth
    td_exceeded <- sum(nuts_df[, "treedepth__"] >= max_td)
    cat("Max treedepth exceeded:", td_exceeded, "(", max_td, ")\n")}

  # E-BFMI
  ebfmi <- sapply(nuts_params, function(chain) {
    sum(diff(chain[, "energy__"])^2) / (var(chain[, "energy__"]) * (nrow(chain) - 1))})
  print(ebfmi)

  # MCMC diagnostics
  mcmc_plot(model, type = "areas", prob = 0.95)
  pairs(model)

  # Multicollinearity
  check_collinearity(model)
}
```

# Model comparisons

## `glmmTMB` specifics

The **Likelihood ratio test** (LRT) was performed with `anova` function of the `stats` package between nested models run with the same data.

```{ruby }
# LRT
anova(null_model, Agr.vs.NonAgr_model)
```

## `brms` specifics

```{ruby }
# ROPE
rope_result <- rope(model)
plot(rope_result)

equivalence_test(model, verbose = FALSE)
```

## Make IBD matrix

```{ruby }
# PLINK file set comprised of .bed, .bim and .fam files
DB=ddPCR_data_set_unrelated

# Filter variants not in chr 1, with any level of missingness or MAF<0.2
plink --bfile ${DB} --allow-no-sex --chr 1 --geno 0 --maf 0.2 --make-bed --out ${DB}_filtered

# Make IBS-based distance matrix in square format
plink --bfile ${DB}_filtered --distance square 1-ibs --out ${DB}_filtered
```

## Make NJ tree

```{ruby }
# Read 1-IBS matrix and format for later use
IBS_mdist <- read_table("ddPCR_data_set_unrelated_filtered.mdist", col_names=FALSE) 
IBS_mdist.id <- read_table("ddPCR_data_set_unrelated_filtered.mdist.id", col_names=FALSE) 
rownames(IBS_mdist) <- IBS_mdist.id$X2
colnames(IBS_mdist) <- IBS_mdist.id$X2

# Reorder IBS matrix to match order in the dataframe containing ddPCR data
ddPCR_data_set <- ddPCR_data_set[match(rownames(IBS_mdist), ddPCR_data_set$sample),]   

# Make tree from 1-IBS matrix with NJ algorithm
tree <- nj(as.dist(IBS_mdist))

# Give tree tip labels and node numbers to identify best node for rooting
tree$tip.label <- ddPCR_modeling_data$pop
tree$node.label <- as.character(1:tree$Nnode)

# Plot tree to visualize node for rooting
plot(tree, "fan", show.node.label=TRUE, use.edge.length=FALSE, align.tip.label=TRUE)

# Select node for rooting
node.tree <- as_tibble(tree)[as_tibble(tree)$label=="130", "node"][[1]] # Ex. node: 130

# Root tree
tree_rooted <- root(tree, node=node.tree, resolve.root=TRUE)

# Chronometricize tree
tree_chrono <- chronos(tree_rooted, model="correlated")
```

## Phylogenetic variance-covariance (VCV) matrix

```{ruby }
# Obtain VCV matrix from rooted phylogenetic tree
VCV <- vcv.phylo(tree_chrono)
```

## Principal component analysis

```{ruby }
# PLINK file set comprised of .bed, .bim and .fam files
DB=ddPCR_data_set_unrelated

# Keep variants in chr 1
plink --bfile ${DB} --chr 1 --out ${DB}_chr1

# Identify variants in LD: 50 kb window size, 10 kb step size, 0.8 r2 threshold
plink --bfile ${DB}_chr1 --indep-pairwise 50 10 0.8 --out ${DB}_LD_results_50_10_0.8

# Remove variables in LD (keep variables in prune.in file)
plink --bfile ${DB}_chr1 --extract ${DB}_LD_results_50_10_0.8.prune.in --make-bed \
    --out ${DB}_chr1_LD_filtered_50_10_0.8

# Run PCR
plink --bfile ${DB}_chr1_LD_filtered_50_10_0.8 --pca \
    --out ${DB}_chr1_LD_filtered_50_10_0.8
```

# **Ancient Eurasians**

## Supervised ADMIXTURE

```{ruby }
# Supervised ADMIXTURE with 3 sources
admixture --cv=10 -j3 --supervised ${FILE}.bed 3
```


# **Out of Africa dispersal**

## Geographical distances between non-Sub-Saharan populations and East Africa

```{ruby }
# Migration waypoints
migration_origin <- c(39.5, 9.0)         
arabian_peninsula <- c(38.5, 35.0)     
turkey <- c(27.0, 38.0)              
thailand <- c(95, 25)               
singapore <- c(103, 1)             
pakistan <- c(68, 28)               
northern_china <- c(125, 50)         
central_alaska <- c( -147.7, 64.8 )   
bering_russia <- c( -171.0, 65.7 )     
washington_center <- c(-120.5, 47.5)

# FUNCTION: calculate geographical distance between specified waypoints
compute_path_dist <- function(waypoints) {
  total_distance <- 0
  for (i in 1:(length(waypoints)-1)) {
    # Add distances between all waypoints
    total_distance <- total_distance + distGeo(waypoints[[i]], waypoints[[i + 1]])}
  # Transform to kilometers
  return(total_distance/1000)}

n_min=5   # Min number of inds per pop

# Select data set to analyse
OOA_dataset <- RD_metadata %>%
  filter(region != "AFR") %>% # Keep non-Africans
  group_by(pop) %>%
  filter(n() >= n_min) %>% # Keep pops with ≥5 inds
  ungroup()

# Get population coordinates
OOA_pop_coords <- OOA_dataset %>%
  group_by(pop) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            region = first(region))

# Calculate distance between East Africa and each population
OOA_geo_distances <- OOA_pop_coords %>%
  rowwise() %>% mutate(
    OOA_pop_coords = list(c(lon, lat)), # Pop coordinates
    migration_distance_km = case_when(  # Calculate OOA distances, varying routes:
      
      # AMERICA
      region=="AMR" ~ 
        compute_path_dist(list(migration_origin, arabian_peninsula, pakistan, 
                               northern_china, bering_russia, central_alaska, 
                               washington_center, OOA_pop_coords)),
      
      # OCEANIA
      region=="OCN" ~ compute_path_dist(list(migration_origin, arabian_peninsula, 
                                             pakistan, thailand, singapore, 
                                             OOA_pop_coords)),
      
      # SOUTH ASIA
      region=="SA" ~ compute_path_dist(list(migration_origin, arabian_peninsula, 
                                            pakistan, OOA_pop_coords)),
      
      # EAST ASIA
      region=="EA" ~ 
        compute_path_dist(list(migration_origin, arabian_peninsula, pakistan, 
                               thailand, OOA_pop_coords)),
      
      # CENTRAL ASIA and SIBERIA
      region=="CAS" ~ 
        compute_path_dist(list(migration_origin, arabian_peninsula, OOA_pop_coords)),
      
      # MIDDLE EAST
      pop %in% c("Druze","Palestinian","Mozabite","Bedouin") ~ 
        compute_path_dist(list(migration_origin, arabian_peninsula, OOA_pop_coords)),
      
      # WESTERN EURASIA
      region=="WEA" ~ 
        compute_path_dist(list(migration_origin, arabian_peninsula, turkey, 
                               OOA_pop_coords))
      
    )) %>% ungroup()
```

## Measures of diversity at the *AMY1* locus

```{ruby }
# Calculate diversity metrics (no rarefication)
OOA_diversity_stats <- OOA_dataset %>%
  group_by(pop) %>%  
  summarise(
    n = n(),
    avg_AMY1 = mean(AMY1_rd),
    SD_AMY1 = sd(AMY1_rd),
    Var_AMY1 = var(AMY1_rd),
    CV_AMY1 = sd(AMY1_rd)/mean(AMY1_rd),
    Range_AMY1 = diff(range(AMY1_rd)),
    UniqueValues_AMY1 = length(unique(AMY1_rd)),
    ShannonEntropy_AMY1 = calc_shannon(AMY1_rd),
    .groups = "drop")

# Calculate diversity metrics (with rarefication)
OOA_diversity_stats_RAREFIED <- OOA_dataset %>%
  group_by(pop) %>%
  do({
    results <- replicate(500, {   # Num of replicates
      samp <- sample_n(., n_min)  # Each replicate based on n_min individuals
      data.frame(
        n = n_min,
        avg_AMY1 = mean(samp$AMY1_rd),
        SD_AMY1 = sd(samp$AMY1_rd),
        Var_AMY1 = var(samp$AMY1_rd),
        CV_AMY1 = sd(samp$AMY1_rd)/mean(samp$AMY1_rd),
        Range_AMY1 = diff(range(samp$AMY1_rd)),
        UniqueValues_AMY1 = length(unique(samp$AMY1_rd)),
        ShannonEntropy_AMY1 = calc_shannon(samp$AMY1_rd))},
      simplify = FALSE)
    bind_rows(results) %>%
      summarise(across(everything(), mean))}) %>%
  ungroup()
```

# **Measures of phylogenetic signal**

```{ruby }
# PREPARE OBJECT: 
# AMY1 CN + POPULATION METADATA to be combined with PHYLOGENETIC TREE
tree <- comparative.data(phy=tree_chrono, 
                         data=RD_modeling_data, 
                         names.col="sample", vcv=TRUE, na.omit=FALSE)
# Create column with sample names
tree$data$sample <- phy$phy$tip.label
```

## Phylogenetic signal

```{ruby }
# Phylogenetic signal
phyloSignal(phylo4d(tree, data["AMY1_CN"]), reps=1000)
```

## Phylogenetic correlogram

```{ruby }
# Phylogenetic correlograms
phyloCorrelogram(phylo4d(tree, data["AMY1_CN"]), trait="AMY1_CN")
```

## Local indicators of phylogenetic association

```{ruby }
# Local indicators of phylogenetic association
lipaMoran(phylo4d(tree, data["AMY1_CN"]), reps=1000, alternative="greater")
lipaMoran(phylo4d(tree, data["AMY1_CN"]), reps=1000, alternative="less")
```


