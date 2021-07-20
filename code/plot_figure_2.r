
red <- rgb(222, 39, 47, 255, maxColorValue = 255)
yellow <-  rgb(254, 229, 157, 255, maxColorValue = 255)
lblue <- rgb(168, 206, 225, 255, maxColorValue = 255)
dblue <- rgb(43, 122, 176, 255, maxColorValue = 255)
                
cols <- c(red, yellow, lblue, dblue, dblue, dblue, "black")
names(cols) <- c("NNRTI", "PI", "NRTI", "M184VI", "3TC", "3TC/FTC", "total")

#g.a is fit in fix_exponential_decay.r

g.B <- Bplot %>% 
    ggplot() +
        geom_bar(aes(x = Reg, y = y, fill = factor(numMuts)),
                 col = "white", stat = "identity", lwd = 0.05) + 
        theme_minimal() + 
        scale_fill_grey(start = 0.95, end = 0) + 
        labs(x = "", y = "Proportion", fill = "Number\nof DRMs") +
        facet_grid(regCats~., scales = "free_y", space = "free_y" ) + 
        coord_flip() +
            theme(
                strip.background = element_blank(),
                strip.text.y = element_blank()
            ) +
                theme(legend.position="bottom")

g.C <- Cplot %>% 
    ggplot() + 
    geom_bar(aes(x = Reg, y =c, fill = drugtype), col = "white",
             stat = "identity", width = 1, lwd = .2) + 
    theme_minimal() +#
    facet_grid(regCats~., scales = "free_y", space = "free_y" ) +
        coord_flip() +
    theme(
                strip.background = element_blank(),
                strip.text.y = element_blank()
            ) + 
    labs(y = "Proportion", fill = "Mutation type") +
    theme(axis.title.y = element_blank()) + 
    scale_fill_manual(values = cols)  +
    theme(legend.position="bottom") + theme(legend.title = element_blank())


g <- plot_grid(g.a, g.B, g.C, nrow = 1, rel_widths = c(0.65, 1, 1), labels = c("A", "B", "C"))
jpeg("../figures/F2_data.jpg", width = 12, height = 3.5, units = "in", res = 300)
print(g)
dev.off()
