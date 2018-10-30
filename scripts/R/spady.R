data.plot <- PT31.all.cells[PT31.all.cells$phase=="G1", ]  %>% 
  dplyr::group_by(phase, time.1.1, stimulation.1.1) %>% 
  dplyr::summarise(mean = mean(log10(Mean_Alexa)), sd = sd(log10(Mean_Alexa)))

plots <- list()
title <- "PT31-image, cells"
plots[["both.phases"]] <- ggplot()+
  geom_line(data=model.all.cells, aes(x = time,
                                      y = log10(value),
                                      group=sigma), alpha = 0.1)+
  geom_errorbar(data=data.plot,
                aes(x = time.1.1, ymin = mean - sd,  ymax = mean + sd, group = time.1.1), colour = "green")+
  geom_point(data=data.plot,
             aes(x = time.1.1, y = mean, group = time.1.1), colour = "brown")+
  theme_jetka(base_size =7)+
  facet_grid(phase ~ stimulation.1.1)+
  ylim(1.5, 3)+
  ggtitle(title)
plots[["both.phases"]]