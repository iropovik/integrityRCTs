library(tidyverse)
wosData <- read_csv("data/wosExport.csv")

set.seed(1)
wosData <- wosData %>%
  select(4, 3, 5, 6, 19, 18, 16, 9, 11, 40, 12 ,21, 8, 36, 1) %>%
  mutate(paperID = sample(x = 1:nrow(wosData), size = nrow(wosData), replace = F)) %>%
  relocate(paperID, .before = 1) %>%
  arrange(paperID)

View(wosData)

write_csv(wosData, "data/wosData.csv")
