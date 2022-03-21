# check relationship between LST, Ta, TS.5 


# Relationship in daily scale

# Ta vs LST
ggplot(data=BB_met, aes(x=LST, y=Ta)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 25) + 
  theme_bw()

# TS.5 vs LST
ggplot(data=BB_met, aes(x=LST, y= TS.5)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 17) + 
  theme_bw()

# Ta vs TS.5
ggplot(data=BB_met, aes(x=LST, y= Ta)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 17) + 
  theme_bw()

# TS.5 vs LST by year
ggplot(data=BB_met, aes(x= (LST), y= (TS.5), color = period)) +
  geom_smooth(method="loess", se = F) +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 12) + 
  theme_bw()

# Ta vs LST by year
ggplot(data=BB_met, aes(x= (LST), y= (Ta), color = period)) +
  geom_smooth(method="loess", se = F) +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 12) + 
  theme_bw()





# Relationship in Monthly scale

# Ta vs LST
ggplot(data=met_monthly, aes(x=LST, y=Ta)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 15) + 
  theme_bw()

# TS.5 vs LST
ggplot(data=met_monthly, aes(x=LST, y= TS.5)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 17) + 
  theme_bw()

# Ta vs TS.5
ggplot(data=met_monthly, aes(x= Ta, y= TS.5)) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y = 17) + 
  theme_bw()



# check relationship between TS.5 & wtd

BB_met %>% filter(season == "growing season") %>% 
  mutate(doy = yday(date)) %>% 
  ggplot() +
  # geom_smooth(method="loess", se = F) +
  geom_line(aes(x= (doy), y= (wtd), color = period),size = 0.9) +
  # geom_line(aes(x = doy, y = TS.5, color = period, size = 0.9))
  # stat_regline_equation() +
  # stat_cor(aes(label=..rr.label..)) + 
  theme_bw()


BB_met %>% filter(season == "growing season") %>% 
  mutate(doy = yday(date)) %>% 
  ggplot() +
  # geom_smooth(method="loess", se = F) +
  # geom_line(aes(x= (doy), y= (wtd), color = period),size = 0.9) +
  geom_line(aes(x = doy, y = rollmean(Ta, 7,fill = NA) , color = period), size = 1) +
  # stat_regline_equation() +
  # stat_cor(aes(label=..rr.label..)) + 
  theme_bw()


BB_met %>% filter(season == "growing season") %>% 
  mutate(doy = yday(date)) %>% 
  ggplot(aes(x = doy, y = TS.5 , color = period)) +
  geom_smooth(method="loess", se = F) +
  # geom_line(aes(x= (doy), y= (wtd), color = period),size = 0.9) +
  geom_point()+
  # stat_regline_equation() +
  # stat_cor(aes(label=..rr.label..)) + 
  theme_bw()

BB_met %>% filter(season == "growing season") %>% 
  mutate(doy = yday(date)) %>% 
  ggplot(aes(x = doy, y = LST , color = period)) +
  geom_smooth(method="loess", se = F) +
  # geom_line(aes(x= (doy), y= (wtd), color = period),size = 0.9) +
  geom_point()+
  # stat_regline_equation() +
  # stat_cor(aes(label=..rr.label..)) + 
  theme_bw()


BB_met_water <- BB_met %>% filter(season == "growing season") %>% 
  mutate(doy = yday(date)) %>% 
  group_by(period) %>% 
  mutate(cumsumET = cumsum(ET)) %>% 
  ungroup()

plot(BB_met_water$wtd, BB_met_water$cumsumET)

BB_met_water %>%
  ggplot(aes(x = doy, y = G , color = period)) +
  geom_smooth(method="loess", span = 0.2,  se = F) +
  # geom_line(aes(x= (doy), y= (wtd), color = period),size = 0.9) +
  geom_point()+
  # stat_regline_equation() +
  # stat_cor(aes(label=..rr.label..)) + 
  theme_bw()
