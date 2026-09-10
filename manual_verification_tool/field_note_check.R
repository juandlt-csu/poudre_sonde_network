library(tidyverse)
library(ross.wq.tools)
library(arrow)
sensor_notes <- load_mWater(creds = yaml::read_yaml(here::here("creds", "mWaterCreds.yml")))

parameter = "depth"
site_oi = "archery"
sensor_trim <- sensor_notes %>%
  filter(site == site_oi) %>%
  mutate(DT_round = with_tz(DT_round, tzone = "America/Denver"))%>%
  mutate(year = year(DT_round))%>%
  select(site, DT_round, crew, year, visit_type, visit_comments, contains(parameter))%>%
  filter(year == 2022)%>%
  # get column with the parameter in the name
  select(DT_round,crew, contains(parameter),  visit_comments)%>%
  arrange(DT_round)

view(sensor_trim)

sensor_data <- read_parquet(file = "data/collated/sensor/compiled_all_sensor_data_2026-01-29.parquet")%>%
  filter(site == "sfm" & parameter == "pH" & year(DT_round) == 2022)


plotly::ggplotly(ggplot(sensor_data%>%filter(!is.na(mean_cleaned)), aes(x = DT_round, y = mean_cleaned))+
  geom_point()+
    geom_line()+
  theme_bw())



grabs <- read_parquet(file = "data/collated/chem/ross_chem_data_2026-01-27.parquet")

plot<- grabs%>%
  filter(site_code %in% c("sfm", "pbr")  & year(DT_sample)== 2025)%>%
ggplot(aes(DT_sample, ChlA, color = site_code))+
  geom_point()+geom_line()

plotly::ggplotly(plot)
