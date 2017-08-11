all_genotyped <- function(attrs) {
	 attrs[['aGenotyped_SC_A']] ==1 && 
	 attrs[['aGenotyped_SC_B']] ==1 && 
	 attrs[['aGenotyped_SC_C']] ==1 && 
	 attrs[['aGenotyped_SC_D']] ==1
}
any_genotyped <- function(attrs) {
	 attrs[['aGenotyped_SC_A']] ==1 || 
	 attrs[['aGenotyped_SC_B']] ==1 || 
	 attrs[['aGenotyped_SC_C']] ==1 || 
	 attrs[['aGenotyped_SC_D']] ==1
}
panel_test <- function(traj,inputs) {
	traj %>% 
	 set_attribute('aGenotyped_SC_A', 1) %>% 
	 set_attribute('aGenotyped_SC_B', 1) %>% 
	 set_attribute('aGenotyped_SC_C', 1) %>% 
	 set_attribute('aGenotyped_SC_D', 1) %>%
	mark("panel_test")
}
