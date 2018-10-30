start <- Sys.time()
model.G1.cells <- prepare_model_par(input.file = file.to.G1.model,
                                    dane=PT31.all.cells[PT31.all.cells$phase=="G1", ], 
                                    factor=factor.for.model,
                                    sigma= -1,
                                    no.cores = 8,
                                    phase= "G1")
end <- Sys.time()
parallel2 <- end-start
# head(model.G1.cells)
model.all.cells <- model.G1.cells

title <- paste("PT31 cells, CV=", CV, "_", abundance, sep= '')
plots[["box_lines"]][[toString(CV)]] <- ggplot()+
  geom_line(data=model.all.cells[model.all.cells$stimulation.1.1 !=0, ], 
            aes(x = time,
                y = log10(value),
                group=sigma), 
            size = 0.8,
            alpha = 0.05)+
  geom_violin(data = PT31.all.cells[PT31.all.cells$phase=="G1" & 
                                      PT31.all.cells$stimulation.1.1 != 0 , ],
              aes(x = time.1.1, y = log10(Mean_Alexa), group = time.1.1),
              fill = "darkgoldenrod1", color= "goldenrod4")+
  theme_jetka(base_size =7)+
  facet_grid(phase ~ stimulation.1.1)+
  ylim(1.2, 3)+
  ggtitle(title)
plots[["box_lines"]][[toString(CV)]]