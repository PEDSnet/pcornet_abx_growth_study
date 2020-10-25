library(tidyverse) # loads up most functions in one cal
library(tableone) # for creating output table
library(scales)
library(DBI)
library(RPostgres)


postgre_table <- function(src, schema, table) {
  # function to get database link, name of schema, name of
  # table, return all

  print(paste("SELECT * FROM", paste(schema, table, sep = ".")))
  paste("SELECT * FROM", paste(schema, table, sep = ".")) %>%
    sql() %>%
    tbl(src = src)
}

# connection details deleted

con <- try(dbConnect(
  RPostgres::Postgres(),
  dbname = my_dbname,
  host = my_host, port = 5432, user = my_user, password = my_password,
  options = "-c search_path=c0_pool"
))

# bring in tables, modify where necessary

cohort <- postgre_table(con, "c0_pool", "c0_cv_demog") %>%
  mutate(race = case_when(race == "02" ~ "Asian", race == "03" ~ "Black", race == "05" ~ "White", race %in% c("01", "04", "06", "08", "OT") ~ "Other or multiracial", race %in%
    c("07", "NI", "UN") ~ "Missing"), hispanic = case_when(hispanic == "Y" ~ 1, TRUE ~ 0)) %>%
  dplyr::select(
    schema, patid,
    race, sex, hispanic, age, age_order, days_btw_rec, ht, wt, cdc_ht_z, cdc_wt_z, cdc_bmi_z
  ) %>%
  compute()


matt_final <- postgre_table(con, "c0_pool", "c0_matt_final_cohort") %>%
  dplyr::select(schema, patid, uniqueID, age, ht, wt, cdc_ht_z, cdc_wt_z) %>%
  mutate(matt_final = 1) %>%
  distinct() %>%
  compute()

cohort <- cohort %>%
  left_join(matt_final) %>%
  distinct() %>%
  compute()

cv_dates <- postgre_table(con, "c0_pool", "c0_table_dates") %>%
  filter(age < 72) %>%
  compute() # created table with all dates

cohort <- cohort %>%
  left_join(cv_dates) %>%
  distinct() %>%
  compute()
# old_bmi<-postgre_table(con, "c0_pool", "c0_cv_demog") %>% filter(age>=24,between(cdc_ht_z, -5, 4),between(cdc_wt_z, -5, 8),between(cdc_ht_z, -5, 4)) %>% dplyr::select(schema,patid,age,ht,wt,cdc_ht_z,cdc_wt_z,bmi,cdc_bmi_z) %>% compute()

# cohort <- inner_join(cohort,old_bmi) %>% distinct() %>% compute()

cohortpresc24 <- postgre_table(con, "c0_pool", "c0_cohortpresc") %>%
  filter(age < 24) %>%
  compute()

cohortencounter24 <- postgre_table(con, "c0_pool", "c0_cohortencounter") %>%
  filter(age < 24, enc_type %in%
    c("AV", "ED", "EI", "IP")) %>%
  compute()

abx_ht_after <- abx_measures %>%
  filter(between(who_ht_z, -5, 4), between(datediff, 14, 168)) %>%
  group_by(schema, patid) %>%
  summarise(n = n()) %>%
  compute()
abx_ht_after %>%
  filter(n >= 2) %>%
  ungroup() %>%
  summarise(n = n())
abx_wt_after <- abx_measures %>%
  filter(between(who_wt_z, -5, 8), between(datediff, 14, 168)) %>%
  group_by(schema, patid) %>%
  summarise(n = n()) %>%
  compute()
abx_wt_after %>%
  filter(n >= 2) %>%
  ungroup() %>%
  summarise(n = n())
abx_both_after <- abx_measures %>%
  filter(between(who_wt_z, -5, 8), between(who_ht_z, -5, 4), between(datediff, 14, 168)) %>%
  group_by(schema, patid) %>%
  summarise(n = n()) %>%
  compute()
abx_both_after %>%
  filter(n >= 2) %>%
  ungroup() %>%
  summarise(n = n())



cohortdx <- postgre_table(con, "c0_pool", "c0_cohortdx") %>%
  dplyr::select(schema, patid, age, age_order, days_btw_rec, dx_group) %>%
  distinct() %>%
  compute()

tiers <- cohortdx %>%
  mutate(Tier = case_when(dx_group == "Tier1_DX" ~ 1, dx_group == "Tier2_DX" ~ 2, dx_group == "Tier3_DX" ~ 3, TRUE ~ 9)) %>%
  compute()
tiers <- tiers %>%
  group_by(schema, patid, age, age_order, days_btw_rec) %>%
  summarise(Tier = min(Tier)) %>%
  compute()

count_tiers <- tiers %>%
  filter(age < 24) %>%
  group_by(Tier, schema, patid) %>%
  summarise(n = n()) %>%
  compute()
count1 <- count_tiers %>%
  ungroup() %>%
  filter(Tier == 1) %>%
  rename(count_tier_1 = n) %>%
  dplyr::select(-Tier) %>%
  compute()
count2 <- count_tiers %>%
  ungroup() %>%
  filter(Tier == 2) %>%
  rename(count_tier_2 = n) %>%
  dplyr::select(-Tier) %>%
  compute()
count3 <- count_tiers %>%
  ungroup() %>%
  filter(Tier == 3) %>%
  rename(count_tier_3 = n) %>%
  dplyr::select(-Tier) %>%
  compute()

cohort <- left_join(cohort, count1) %>% mutate(count_tier_1 = case_when(is.na(count_tier_1) ~ 0, TRUE ~ round(count_tier_1)))
cohort <- left_join(cohort, count2) %>% mutate(count_tier_2 = case_when(is.na(count_tier_2) ~ 0, TRUE ~ round(count_tier_2)))
cohort <- left_join(cohort, count3) %>%
  mutate(count_tier_3 = case_when(is.na(count_tier_3) ~ 0, TRUE ~ round(count_tier_3))) %>%
  compute()


cohortdxflags <- postgre_table(con, "c0_pool", "c0_cohortdxflags") %>%
  filter(age < 72)

# get infection episodes
infect <- cohortdx %>%
  filter(age < 24, dx_group %in% c(
    "Tier1_DX",
    "Tier2_DX", "Tier3_DX"
  )) %>%
  inner_join(cv_dates) %>%
  group_by(
    schema,
    patid
  ) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~ 0, tmp <= 14 ~ 0, TRUE ~
  date)) %>%
  dplyr::select(schema, patid, tmp) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_infect = n())

cohort <- left_join(cohort, infect) %>% mutate(count_infect = case_when(is.na(count_infect) ~
0, TRUE ~ round(count_infect)))

# get visits
visits <- cohortencounter24 %>%
  dplyr::select(
    schema, patid,
    age, age_order
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_visits = n()) %>%
  compute()

cohort <- left_join(cohort, visits) %>% mutate(count_visits = case_when(is.na(count_visits) ~
0, TRUE ~ round(count_visits)))

# get conditions
preterm <- cohortdx %>%
  filter(age < 24, dx_group == "Preterm") %>%
  dplyr::select(schema, patid) %>%
  distinct() %>%
  mutate(preterm = 1) %>%
  compute()

cohort <- cohort %>%
  left_join(preterm) %>%
  mutate(preterm = case_when(is.na(preterm) ~
  0, TRUE ~ preterm))

asthma <- postgre_table(con, "c0_pool", "c0_cohortdx") %>%
  filter(age < 72, dx_group == "Asthma") %>%
  dplyr::select(schema, patid, age, age_order, days_btw_rec) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(n = n()) %>%
  filter(n >= 2) %>%
  dplyr::select(-n) %>%
  mutate(has_asthma = 1) %>%
  compute()

cohort <- cohort %>%
  left_join(asthma) %>%
  mutate(has_asthma = case_when(is.na(has_asthma) ~
  0, TRUE ~ has_asthma))

growth_cond <- cohortdxflags %>%
  filter(age < 72) %>%
  dplyr::select(
    schema,
    patid, dx_group
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_growth_cond = n()) %>%
  compute()

cohort <- left_join(cohort, growth_cond) %>% mutate(count_growth_cond = case_when(is.na(count_growth_cond) ~
0, TRUE ~ round(count_growth_cond)))

# get reflux

reflux <- cohortpresc24 %>%
  filter(age <= 24, rxnorm_cui_group ==
    "Reflux_meds_all") %>%
  inner_join(cv_dates) %>%
  group_by(
    schema,
    patid
  ) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  compute()

reflux10 <- reflux %>%
  mutate(tmp = case_when(
    is.na(tmp) ~ 0,
    tmp <= 10 ~ 0, TRUE ~ date
  )) %>%
  dplyr::select(
    schema, patid,
    tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_reflux_10 = n())

cohort <- left_join(cohort, reflux10) %>% mutate(count_reflux_10 = case_when(is.na(count_reflux_10) ~
0, TRUE ~ round(count_reflux_10)))


reflux30 <- reflux %>%
  mutate(tmp = case_when(
    is.na(tmp) ~ 0,
    tmp <= 30 ~ 0, TRUE ~ date
  )) %>%
  dplyr::select(
    schema, patid,
    tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_reflux_30 = n())

cohort <- left_join(cohort, reflux30) %>% mutate(count_reflux_30 = case_when(is.na(count_reflux_30) ~
0, TRUE ~ round(count_reflux_30)))


# get Abx
Abx_Type <- read_csv("Broad vs Narrow.csv")
names(Abx_Type) <- tolower(names(Abx_Type))
Abx_Type$rxnorm_cui <- as.character(Abx_Type$rxnorm_cui)

abx <- cohortpresc24 %>%
  filter(rxnorm_cui_group == "Antibiotics_all") %>%
  inner_join(Abx_Type, copy = TRUE) %>%
  inner_join(cv_dates) %>%
  compute()

v1 <- abx %>%
  dplyr::select(schema, patid, abx_category, date) %>%
  distinct()
v2 <- v1 %>%
  group_by(schema, patid, date) %>%
  summarise(n = n()) %>%
  filter(n == 2)

abx_classification <- full_join(v1, v2) %>%
  mutate(abx_class = case_when(is.na(n) ~ abx_category, TRUE ~ "Broad")) %>%
  select(-n) %>%
  distinct() %>%
  compute()

# abx_both_types <- abx %>% dplyr::select(schema, patid, abx_category) %>%
#  distinct() %>% group_by(schema, patid) %>% summarise(n = n()) %>%
#  filter(n == 2) %>% mutate(has_both = 1) %>% dplyr::select(-n) %>%
#  inner_join(abx)

# abx_classification <- abx %>% dplyr::select(schema, patid, abx_category) %>%
#  group_by(schema, patid, abx_category) %>% summarise(n = n()) %>%
#  full_join(abx_both_types) %>% mutate(abx_class = case_when(has_both ==
#  1 ~ "Broad", TRUE ~ abx_category)) %>% dplyr::select(schema,
#  patid, abx_category, abx_class) %>% distinct() %>% compute()

abx_all <- abx %>%
  dplyr::select(schema, patid, date, age, abx_category) %>%
  inner_join(abx_classification) %>%
  compute()

abx_age <- abx_all %>%
  group_by(schema, patid) %>%
  summarise(age_first_abx = min(age), age_last_abx = max(age))
broad_age <- abx_all %>%
  filter(abx_class == "Broad") %>%
  group_by(schema, patid) %>%
  summarise(age_first_broad = min(age), age_last_broad = max(age))
narrow_age <- abx_all %>%
  filter(abx_class == "Narrow") %>%
  group_by(schema, patid) %>%
  summarise(age_first_narrow = min(age), age_last_narrow = max(age))

cohort <- left_join(cohort, abx_age) %>%
  left_join(broad_age) %>%
  left_join(narrow_age) %>%
  compute()

#_10 and _14 refers to implementation of 10-day and 14-day windows

abx_10 <- abx_all %>%
  dplyr::select(schema, patid, date) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~
  0, tmp <= 10 ~ 0, TRUE ~ date)) %>%
  dplyr::select(
    schema,
    patid, tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_abx_10 = n()) %>%
  compute()

cohort <- left_join(cohort, abx_10) %>% mutate(count_abx_10 = case_when(is.na(count_abx_10) ~
0, TRUE ~ round(count_abx_10)))


abx_14 <- abx_all %>%
  dplyr::select(schema, patid, date) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~
  0, tmp <= 14 ~ 0, TRUE ~ date)) %>%
  dplyr::select(
    schema,
    patid, tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_abx_14 = n()) %>%
  compute()

cohort <- left_join(cohort, abx_14) %>% mutate(count_abx_14 = case_when(is.na(count_abx_14) ~
0, TRUE ~ round(count_abx_14)))

broad_10 <- abx_all %>%
  filter(abx_class == "Broad") %>%
  dplyr::select(
    schema,
    patid, date
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~
  0, tmp <= 10 ~ 0, TRUE ~ date)) %>%
  dplyr::select(
    schema,
    patid, tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_broad_10 = n()) %>%
  compute()

cohort <- left_join(cohort, broad_10) %>% mutate(count_broad_10 = case_when(is.na(count_broad_10) ~
0, TRUE ~ round(count_broad_10)))

broad_14 <- abx_all %>%
  filter(abx_class == "Broad") %>%
  dplyr::select(
    schema,
    patid, date
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~
  0, tmp <= 14 ~ 0, TRUE ~ date)) %>%
  dplyr::select(
    schema,
    patid, tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_broad_14 = n()) %>%
  compute()

cohort <- left_join(cohort, broad_14) %>% mutate(count_broad_14 = case_when(is.na(count_broad_14) ~
0, TRUE ~ round(count_broad_14)))

narrow_10 <- abx_all %>%
  filter(abx_class == "Narrow") %>%
  dplyr::select(schema, patid, date) %>%
  distinct() %>%
  group_by(
    schema,
    patid
  ) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~ 0, tmp <= 10 ~ 0, TRUE ~
  date)) %>%
  dplyr::select(schema, patid, tmp) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_narrow_10 = n()) %>%
  compute()

cohort <- left_join(cohort, narrow_10) %>% mutate(count_narrow_10 = case_when(is.na(count_narrow_10) ~
0, TRUE ~ round(count_narrow_10)))

narrow_14 <- abx_all %>%
  filter(abx_class == "Narrow") %>%
  dplyr::select(schema, patid, date) %>%
  distinct() %>%
  group_by(
    schema,
    patid
  ) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~ 0, tmp <= 14 ~ 0, TRUE ~
  tmp)) %>%
  dplyr::select(schema, patid, tmp) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_narrow_14 = n()) %>%
  compute()

cohort <- left_join(cohort, narrow_14) %>% mutate(count_narrow_14 = case_when(is.na(count_narrow_14) ~
0, TRUE ~ round(count_narrow_14)))

# get steroids

Csteroid_Type <- read_csv("Csteroid Exclude.csv")
Csteroid_Type$rxnorm_cui <- as.character(Csteroid_Type$rxnorm_cuix)
OralSteroid <- filter(Csteroid_Type, exclude_steriods == 0)
InhaledSteroid <- filter(Csteroid_Type, exclude_steriods != 0)


csteroid <- cohortpresc24 %>%
  filter(rxnorm_cui_group == "Corticosteroids_all") %>%
  inner_join(Csteroid_Type, copy = TRUE) %>%
  inner_join(cv_dates) %>%
  compute()


cortico_all <- csteroid %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date - lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~
  0, tmp <= 10 ~ 0, TRUE ~ date)) %>%
  dplyr::select(
    schema,
    patid, tmp
  ) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_cortico_all = n())

cohort <- left_join(cohort, cortico_all) %>% mutate(count_cortico_all = case_when(is.na(count_cortico_all) ~
0, TRUE ~ round(count_cortico_all)))

cortico_oral <- csteroid %>%
  inner_join(OralSteroid, copy = TRUE) %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date -
    lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~ 0, tmp <=
    10 ~ 0, TRUE ~ date)) %>%
  dplyr::select(schema, patid, tmp) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_cortico_oral = n())

cohort <- left_join(cohort, cortico_oral) %>% mutate(count_cortico_oral = case_when(is.na(count_cortico_oral) ~
0, TRUE ~ round(count_cortico_oral)))
cortico_inhaled <- csteroid %>%
  inner_join(InhaledSteroid, copy = TRUE) %>%
  group_by(schema, patid) %>%
  arrange(date) %>%
  mutate(tmp = date -
    lag(date)) %>%
  mutate(tmp = case_when(is.na(tmp) ~ 0, tmp <=
    10 ~ 0, TRUE ~ date)) %>%
  dplyr::select(schema, patid, tmp) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(count_cortico_inhaled = n())

cohort <- left_join(cohort, cortico_inhaled) %>% mutate(count_cortico_inhaled = case_when(is.na(count_cortico_inhaled) ~
0, TRUE ~ round(count_cortico_inhaled)))


chronicdx <- postgre_table(con, "c0_pool", "c0_cohortdxflags") %>%
  filter(age >= 0, age < 72) %>%
  dplyr::select(schema, patid, dx_group, age, age_order, days_btw_rec) %>%
  distinct() %>%
  group_by(schema, patid) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  dplyr::select(schema, patid) %>%
  mutate(has_chronic = 1) %>%
  compute()

cohort <- left_join(cohort, chronicdx) %>%
  mutate(has_chronic = case_when(has_chronic == 1 ~ 1, TRUE ~ 0)) %>%
  compute()

cohort %>% compute(name = "c0_final_table_march12", temporary = FALSE)

check2 <- cohort %>%
  filter(matt_final == 1) %>%
  select(schema, patid, date) %>%
  ungroup() %>%
  group_by(schema, patid) %>%
  summarise(minDate = min(date), maxDate = max(date)) %>%
  filter(maxDate - minDate >= 180) %>%
  select(schema, patid) %>%
  distinct() %>%
  mutate(check2 = 1) %>%
  compute()

cohort2 <- cohort %>%
  rename(matt_cohort = matt_final) %>%
  mutate(race_hispanic = case_when(hispanic == 1 ~ "Hispanic", race == "Missing" ~ NA, TRUE ~ race)) %>%
  left_join(check2) %>%
  compute()

cohort2 %>% compute(name = "c0_finaltable_march13", temporary = FALSE)
