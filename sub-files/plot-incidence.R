#prep dataset and rank by event A time
rk <- results %>% filter(grepl("^A_SC",resource) | resource=="time_in_model") %>% 
  group_by(name,preemptive,reactive) %>%
  mutate(rank=dense_rank(start_time)-1) #rank Xth event A

#Plot: 3N as denominator of cumulative incidence
#each run simulates the same cohort under 3 strategies, 
#experience of A can be different because of downstream B death competing with A risks
plot_incidence <- function(raw,threshold,force=FALSE) {
  tt <- raw %>% ungroup() %>% filter(rank==threshold) %>% 
    select(time=start_time) %>% 
    group_by(time) %>% summarise(ct=n()) %>%
    rbind(data.frame(time=c(0,inputs$vHorizon*365),ct=c(0,0))) %>% 
    arrange(time) %>% mutate(cum=cumsum(ct),incidence=cum/(3*inputs$vN)) 
  
  p <- ggplot(tt,aes_string(x="time",y="incidence")) +
    geom_step() + theme_bw() +
    ggtitle(paste0("Cumulative Incidence of Experiencing >= ",threshold," Event A")) +
    labs(x="Days",y="Probability") + 
    theme(legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  
  #plot 4 can be hard to read if y scale forced to 0-1
  if(force==TRUE) {
    p+scale_y_continuous(limits = c(0, 1)) 
  } else { p }
}

nm <- nrow(scenario.mapping) 
destdir <- getwd() ###need to define destination directory
for(i in 1:nm) {
  png(paste0(destdir,"/plot_",i,".png"))
  p <- plot_incidence(rk,i,FALSE) 
  print(p)
  dev.off()    
}
