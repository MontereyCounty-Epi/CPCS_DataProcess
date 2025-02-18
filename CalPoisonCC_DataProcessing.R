# PCC file processing

library(tidyverse)
library(pdftools)

# set working directory to location of PDF files
setwd("K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files")
getwd()

# import all PDF files in the root (non-archive) folder
file <- list.files(pattern = "pdf$")
table <- lapply(paste0(getwd(), "/", file), pdf_text)
pcc <- lapply(paste0(getwd(), "/", file), pdf_text)%>% unlist() %>% as.data.frame() %>% mutate(across(everything(), as.character))
"text" -> names(pcc) 
#pcc <- pcc %>% subset(grepl("California Poison Control System", text) == TRUE)#%>% str_split("\n")
pcc <- paste(unlist(pcc), collapse = ' ; ')

#pcc2 <- paste(unlist(pcc2), collapse = "$$") %>% as.data.frame() %>% mutate(text = str_split(.,";", n = Inf)) %>% rowwise() %>% unnest(text)

# extract file metadata and split into sections
pcc2 <- pcc %>% as.data.frame() %>%
  mutate(text = str_split(.,";", n = Inf)) %>% select(-`.`) %>% 
  rowwise() %>%
  unnest(text) %>% as.data.frame() %>% 
  mutate(text2 = str_replace_all(text, "[[:cntrl:]]", " ; ")) %>% #str_split("\t") %>% # something is worng with file count -- 1 missing???
  rowwise() %>%
  separate(text2, c("date", "data"), "Reason for Exposure[[:punct:]]\\s") %>% #"\\(06053\\)\\s[[:punct:]]\\s") %>%
  mutate(reporting_date = str_extract(data, "\\d{1,2}/\\d{1,2}/\\d{4}"),
       #  data = str_replace_all(data, "\n", ""),
         unintentional = str_extract(data, "Unintentional (.+)Intentional."),
         unintentional = if_else(is.na(unintentional),str_extract(data, "Unintentional (.+)Other."),unintentional),
         unintentional = if_else(is.na(unintentional),str_extract(data, "Unintentional (.+)Call Volume."),unintentional),
         intentional = str_extract(data, "Intentional(.+)Call Volume.*"),
         intentional = if_else(is.na(intentional),str_extract(data, "Intentional (.+)Other.*"),intentional),
         Other = str_extract(data, "Other(.+)Call Volume.*")
    ) 


# extract cll details from each section
pcc3 <- pcc2 %>% mutate(unintentional = gsub("(?=;\\s[A-Z])", "|", unintentional, perl = TRUE),
                    intentional = gsub("(?=;\\s[A-Z])", "|", intentional, perl = TRUE),
                    Other = gsub("(?=;\\s[A-Z])", "|", Other, perl = TRUE)) %>% 
  separate_wider_delim(unintentional, delim = "|;", names_sep = "_", too_few = "align_start") %>% 
  separate_wider_delim(intentional, delim = "|;", names_sep = "_", too_few = "align_start") %>% 
  separate_wider_delim(Other, delim = "|;", names_sep = "_", too_few = "align_start") %>% select(-c(text, date)) %>% 
  group_by(reporting_date, data) %>% 
  pivot_longer(cols = -c(reporting_date, data), names_to = "type", values_to = "calls") %>% rowwise() %>% 
  mutate(calls = str_trim(calls, "both"),
         remove = if_else(grepl("^Unintentional|^Intentional|^Other|^Total|^MONTEREY|^Call Volume|^For additional",calls),1,NA),
         remove = if_else(is.na(remove) & is.na(calls),1,remove)
         ) %>% subset(is.na(remove))

# clean and format into separate variables
pcc4 <- pcc3 %>% mutate(calls = str_split(calls, "\\s{2,}")) %>% unnest_wider(calls, names_sep = "_") %>% 
  separate(calls_2, c("age", "age_unit", "sex"), "\\s") %>% 
  separate(calls_3, c("case_number", "datetime"), ";") %>% 
  mutate(date = str_extract(datetime, "\\d{1,2}/\\d{1,2}/\\d{4}"),
         time = str_extract(datetime, "\\d{1,2}:\\d{2}:\\d{2}"),
         time_units = str_extract(datetime, "[A-Z]M$")) %>% 
  select(data:datetime, date, time, time_units, calls_4:last_col()) %>% 
  separate(calls_4, c("zip", "city"), "-") %>% 
  separate(city, c("city", "remove_1"), ";") %>% 
  separate(calls_5, c("substance_1", "remove_2"), ";") %>% 
  separate(substance_1, c("remove_3", "substance_1"), ":") %>% select(-c(starts_with("remove"))) %>% rename(outcome = calls_1)

# clean substance variables
pcc5 <- pcc4 %>% 
  mutate_at(vars(calls_6:last_col()), ~ str_replace_all(., "^Grand|^[:digit:]", NA_character_)) %>% 
  mutate_at(vars(calls_6:last_col()), ~ gsub(".*[:]([^.]+)[;].*", "\\1", .)) %>% 
  mutate_at(vars(calls_6:last_col()), ~ gsub(".*[:]([^.]+).*", "\\1", .)) %>% 
  mutate(across(where(is.character), str_trim)) %>% select(-c(data)) %>% 
  mutate(type = str_replace(type, "[_][0-9]", "")) %>% distinct() %>% janitor::remove_empty("cols")

pcc5 <- pcc5[!duplicated(pcc5[, c("reporting_date", "case_number", "datetime")], fromLast=T),]
names(pcc5) <- c('reporting_date', 'type', 'outcome', 'age', 'age_unit', 'sex', 'case_number', 'datetime', 'date', 'time', 'time_units', 'zip', 'city', paste("substance", 1:100, sep="_"))

#colnames(pcc5) <- gsub("calls_", "substance_", colnames(pcc5))
pcc5 <- pcc5 %>% mutate(reporting_date = mdy(reporting_date)) %>% janitor::remove_empty("cols")
# combine with history file for complete dataset

#save(pcc5, file = "K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Datasets/2023.Rda")
# write 2023
historical <- read.csv("K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Datasets/Historical.csv")
historical <- historical %>% janitor::remove_empty("cols") %>% mutate(reporting_date = as.Date(reporting_date))
dataset <- bind_rows(pcc5, historical) %>% arrange(reporting_date) 
dataset <- dataset %>% distinct(case_number, .keep_all = TRUE)
dataset[dataset==""] <- NA
#save(dataset, "K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Data.Rda")

#write.csv(pcc5, "K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Datasets/2023.csv", row.names = FALSE, quote = FALSE, na= "")
write.csv(dataset, "K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Datasets/All.csv", row.names = FALSE, quote = FALSE, na= "")



write.csv(dataset, "K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/Datasets/Historical.csv", row.names = FALSE, quote = FALSE, na= "")
# old attempt -- didn't work
#str_extract_all(., ".*[:](.+)[;].*"))
# unintent2 <- x2 %>% select(reporting_date, unintentional) %>% separate(unintentional, c("meta", "calls"), ";", extra = "merge") %>% 
#   mutate(new = gsub("(?=;\\s[A-Z])", "|", calls, perl = TRUE)) %>%
#    separate_wider_delim(new, delim = "|;", names_sep = "caller_", too_few = "align_start") %>% select(-calls) %>% 
#   group_by(reporting_date, meta) %>%
#   pivot_longer(starts_with("newcaller_"), names_to = "caller", values_to = "data") %>% rowwise() %>%
# #  str_replace_all(data, "\\s{2,}", "$") ## delimiter not working...
#   separate_wider_regex(data, "\\s{2,}", names_sep = "var_", too_few = "align_start")
  
#   str_split("\n") 
# unintent2 <- tibble(unintent2) %>% unnest_wider(unintent2, names_sep = "_") %>%
#   pivot_longer(everything(), names_to = "id", values_to = "data") %>% 
#   separate_rows(data, sep = ";")
#   
#   paste(unlist(.), collapse = '|') %>% tibble(.) %>% tidyr::unnest_wider(.)
# #unlist() %>% as.data.frame()
# unintent2 <- paste(unlist(unintent2), collapse = '|') %>% as.data.frame() %>% 
#   mutate(text = str_split(.,"Unintentional", n = Inf)) %>% rowwise() %>% 
#   unnest(text) %>% separate_rows(everything())#separate_longer_delim(text, delim = ";")
# 
# unintent3 <- with(unintent2, strsplit(text, ';'))
# unintent32 <- data.frame(id = factor(rep(unintent2$id, times = lengths(unintent3)),
#                               levels = unintent2$id),
#                         caller = unlist(unintent3))
# 
# %>% rowwise() %>%
#    %>% 

#y <- pdf_text("K:/4000-PHB/COMMUNICABLE DISEASE/Epi and Surv/Cal-Poison/Data Files/CPCS-Monterey County-Daily DSAT Support Report_11.1.23.pdf") %>% str_split("\n")
