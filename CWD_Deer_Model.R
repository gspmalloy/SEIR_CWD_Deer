require('deSolve')
require('ggplot2')
#Original Data
df_1999_2019_SouthernFarmland <- read.csv('C:/Users/Giovanni Malloy.000/Desktop/Work_103020/USB 031420/Research_Chronic_Wasting_Disease/Data Wisconsin Wild White Tailed Deer/CWD Year 1999-2019 Southern Farmland Population Only.csv')

##################
# Section 1. Model
##################
CWD_mod <- function(Time, State, Pars)
{
  with(as.list(c(State,Pars)),{
    N_Deer <- S_Deer + E_Deer + I_Deer
    births_in <- death_rate_S*S_Deer + hunt_kill_rate*S_Deer*N_Hunters + death_rate_E*E_Deer + hunt_kill_rate*E_Deer*N_Hunters + death_rate_I*I_Deer + hunt_kill_rate*I_Deer*N_Hunters
    d_S_Deer <- births_in - death_rate_S*S_Deer - beta_CWD_Deer*S_Deer*(I_Deer)/N_Deer - beta_CWD_Natural*S_Deer/N_Deer - hunt_kill_rate*S_Deer*N_Hunters
    d_E_Deer <- beta_CWD_Deer*S_Deer*(I_Deer)/N_Deer + beta_CWD_Natural*S_Deer/N_Deer - epsilon_Deer*E_Deer - death_rate_E*E_Deer - hunt_kill_rate*E_Deer*N_Hunters
    d_I_Deer <- epsilon_Deer*E_Deer - death_rate_I*I_Deer - hunt_kill_rate*I_Deer*N_Hunters
    return(list(c(d_S_Deer, d_E_Deer, d_I_Deer)))
  })
}

##########################
# Section 2. Parameters
##########################
S_Deer <- 1377100/4
E_Deer <- 0
I_Deer <- 3
births_in <- 100
life_expectancy_deer <- 4.5 # years https://www.jstor.org/stable/3803059?seq=7#metadata_info_tab_contents
life_expectancy_deer <- life_expectancy_deer * 12
death_rate_S <- 1/life_expectancy_deer
epsilon_Deer <- 1/15 #https://pdfs.semanticscholar.org/75eb/8b27d8cd507c23f2a74bfc7f4391505a7b4a.pdf
death_rate_E <- 1/(life_expectancy_deer)
risk_ratio_CWD_Deer <- 4.5 # https://www.dec.ny.gov/docs/wildlife_pdf/cwdfactsheet.pdf
death_rate_I <- death_rate_S * risk_ratio_CWD_Deer

hunter_harvest_2017 <- 322054 #https://dnr.wi.gov/wideermetrics/DeerStats.aspx
Avg_days_hunted_2017 <- 4.3 # https://dnr.wi.gov/topic/WildlifeHabitat/documents/reports/gundeer2.pdf
hunter_days_2017 <- 2695768 # https://dnr.wi.gov/topic/WildlifeHabitat/documents/reports/gundeer2.pdf
num_hunters_2017 <- 2695768/4.3
num_deer_2017 <- 1377100
num_deers_per_hunter_2017 <- hunter_harvest_2017/num_hunters_2017
hunter_kill_rate_year <- num_deers_per_hunter_2017/num_deer_2017
hunter_kill_rate_month <- hunter_kill_rate_year/12
beta_CWD_Natural <- .000001

##########################
# Section 3. Calibration
##########################
beta_CWD_Deer <- 0.01

# Define list of possible values for beta_CWD_Deer
list_beta_CWD_Deer <- seq(0.01,0.5,by = 0.001)
# This is a placeholder value that also serves as a sanity check
best_beta_CWD_Deer <- 100
# Set the best score to something absurdly high as a santiy check
best_beta_CWD_Deer_Score <- 10000000000000000000000000

#Track for plot
list_all_betas <- c()
list_all_mse_scores <- c()

for(i in 1:length(list_beta_CWD_Deer))
{
  # Test the possible value of beta_CWD_Deer
  beta_CWD_Deer <- list_beta_CWD_Deer[i]
  #Initialize the model population
  ini_cwd_mod <- c(S_Deer = S_Deer, E_Deer = E_Deer, I_Deer = I_Deer)
  #Run the model for 252 months
  times_cwd_mod <- seq(1,252,by = 1) #By Month
  #Initialize the parameters in the model
  pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                    beta_CWD_Deer = beta_CWD_Deer, 
                    beta_CWD_Natural = beta_CWD_Natural, 
                    hunt_kill_rate = hunter_kill_rate_month,
                    epsilon_Deer = epsilon_Deer, 
                    death_rate_E = death_rate_E,
                    death_rate_I = death_rate_I, 
                    N_Hunters = num_hunters_2017)
  
  #Output using our previously defined function
  CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod, pars_cwd_mod)
  
  #create a data frame of the results
  yS_Deer <- CWD_mod_Out[,"S_Deer"]
  yE_Deer <- CWD_mod_Out[,"E_Deer"]
  yI_Deer <- CWD_mod_Out[,"I_Deer"]
  CWD_Results <- data.frame(yS_Deer, yE_Deer, yI_Deer)
  #Calculate the model predicted prevalence each year
  CWD_Results$Prevalence <- CWD_Results$yI_Deer / (CWD_Results$yS_Deer +
                                                     CWD_Results$yE_Deer +
                                                     CWD_Results$yI_Deer)
  CWD_Results_Year <- CWD_Results[c(12,24,36,48,60,72,84,96,108,120,
                                    132,144,156,168,180,192,204,
                                    216,228,240,252),]
  # Compare to the data
  CWD_Results_Year$Prevalence_Data <- df_1999_2019_SouthernFarmland$Incidence
  CWD_Results_Year$Year <- df_1999_2019_SouthernFarmland$Year
  
  CWD_Results_Year$Prevalence_Data <- as.numeric(as.character(CWD_Results_Year$Prevalence_Data))
  CWD_Results_Year$Prevalence <- as.numeric(as.character(CWD_Results_Year$Prevalence))
  # Calculate the mean squared error
  curr_mean_sq_err <- mean((CWD_Results_Year$Prevalence - CWD_Results_Year$Prevalence_Data)^2)
  
  list_all_betas <- c(list_all_betas, beta_CWD_Deer)
  list_all_mse_scores <- c(list_all_mse_scores, curr_mean_sq_err)
  
  #If it is the lowest MSE, save the beta value and the score as the best
  if(curr_mean_sq_err < best_beta_CWD_Deer_Score)
  {
    best_beta_CWD_Deer <- beta_CWD_Deer
    best_beta_CWD_Deer_Score <- curr_mean_sq_err
    print(curr_mean_sq_err)
  }
  print(i)
}

plot_df <- as.data.frame(matrix(c(list_all_betas, list_all_mse_scores), ncol = 2, byrow = FALSE))
colnames(plot_df) <- c('beta_CWD_Deer', 'Model_MSE')

ggplot()+
  geom_line(data = plot_df, aes(x = beta_CWD_Deer, y = Model_MSE, color = 'Model MSE'), size = 3)+
  xlab('beta_CWD_Deer Value')+
  ylab('MSE')+
  ggtitle('Calibration Results')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        title = element_text(size = 16),
        legend.text = element_text(size = 24))+
  scale_color_brewer(palette="Pastel1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(color=' ', shape = '')

###############################################
# Section 4. Plot Results
###############################################

beta_CWD_Deer <- best_beta_CWD_Deer

ini_cwd_mod <- c(S_Deer = S_Deer, E_Deer = E_Deer, I_Deer = I_Deer)

times_cwd_mod <- seq(1,252,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate = hunter_kill_rate_month,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results$Prevalence <- CWD_Results$yI_Deer / (CWD_Results$yS_Deer +
                                                  CWD_Results$yE_Deer +
                                                  CWD_Results$yI_Deer)
CWD_Results_Year <- CWD_Results[c(12,24,36,48,60,72,84,96,108,120,
                                  132,144,156,168,180,192,204,
                                  216,228,240,252),]
CWD_Results_Year$Prevalence_Data <- df_1999_2019_SouthernFarmland$Incidence
CWD_Results_Year$Year <- df_1999_2019_SouthernFarmland$Year

CWD_Results$Year <- seq((1998+(1/12)),2019,by=(1/12))

ggplot()+
  geom_point(data = df_1999_2019_SouthernFarmland, aes(x = Year, y = Incidence, color = 'Southern Farmland Zone Data'), size = 3)+
  geom_line(data = CWD_Results, aes(x = Year, y = Prevalence, color = 'Model'), size = 2)+
  xlab('Year')+
  ylab('Prevalence')+
  ggtitle('Chronic Wasting Disease, Whitetail Deer, Wisconsin')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        title = element_text(size = 16),
        legend.text = element_text(size = 24),
        legend.position = c(0.3,0.9))+
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15),
                     labels = c('0%','5%','10%','15%'))+
  scale_x_continuous(breaks = c(1999,2005,2011,2017),
                     labels = c('1999', '2005', '2011', '2017'))+
  scale_color_brewer(palette="Pastel1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(color=' ', shape = '')

###############################################
# Section 5. Plot Predicted prevalence 10 years forward
###############################################

beta_CWD_Deer <- best_beta_CWD_Deer

ini_cwd_mod <- c(S_Deer = S_Deer, E_Deer = E_Deer, I_Deer = I_Deer)

times_cwd_mod <- seq(1,492,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate = hunter_kill_rate_month,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results$Prevalence <- CWD_Results$yI_Deer / (CWD_Results$yS_Deer +
                                                   CWD_Results$yE_Deer +
                                                   CWD_Results$yI_Deer)
CWD_Results_Year <- CWD_Results[c(12,24,36,48,60,72,84,96,108,120,
                                  132,144,156,168,180,192,204,
                                  216,228,240,252,
                                  264,276,288,300,
                                  312,324,336,348,360,372,
                                  384,396,408,420,432,
                                  444,456,468,480,492),]

CWD_Results$Year <- seq((1998+(1/12)),2039,by=(1/12))

ggplot()+
  geom_point(data = df_1999_2019_SouthernFarmland, aes(x = Year, y = Incidence, color = 'Southern Farmland Zone Data'), size = 3)+
  geom_line(data = CWD_Results, aes(x = Year, y = Prevalence, color = 'Model'), size = 2)+
  xlab('Year')+
  ylab('Prevalence')+
  ggtitle('Chronic Wasting Disease, Whitetail Deer, Wisconsin')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        title = element_text(size = 16),
        legend.text = element_text(size = 24),
        legend.position = c(0.3,0.9))+
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15),
                     labels = c('0%','5%','10%','15%'))+
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020,2025,2030,2035),
                     labels = c('2000','2005','2010','2015','2020','2025','2030','2035'))+
  scale_color_brewer(palette="Pastel1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(color=' ', shape = '')

###################################
# Section 6. Policy Evaluation
###################################

# Run period with data
beta_CWD_Deer <- best_beta_CWD_Deer

ini_cwd_mod <- c(S_Deer = S_Deer, E_Deer = E_Deer, I_Deer = I_Deer)

times_cwd_mod <- seq(1,252,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate = hunter_kill_rate_month,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results_data <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results_data$Prevalence <- CWD_Results_data$yI_Deer / (CWD_Results_data$yS_Deer +
                                                             CWD_Results_data$yE_Deer +
                                                             CWD_Results_data$yI_Deer)
CWD_Results_data$Year <- seq((1998+(1/12)),2019,by=(1/12))

# Baseline model continues 20 years
beta_CWD_Deer <- best_beta_CWD_Deer

ini_cwd_mod <- c(S_Deer = CWD_Results_data$yS_Deer[252], E_Deer = CWD_Results_data$yE_Deer[252], I_Deer = CWD_Results_data$yI_Deer[252])

times_cwd_mod <- seq(253,492,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate = hunter_kill_rate_month,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results_baseline_projection <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results_baseline_projection$Prevalence <- CWD_Results_baseline_projection$yI_Deer / (CWD_Results_baseline_projection$yS_Deer +
                                                                                           CWD_Results_baseline_projection$yE_Deer +
                                                                                           CWD_Results_baseline_projection$yI_Deer)
CWD_Results_baseline_projection$Year <- seq((2019+(1/12)),2039,by=(1/12))

# Scenario 1 model continues 20 years
CWD_mod_diff_hunt_rates <- function(Time, State, Pars)
{
  with(as.list(c(State,Pars)),{
    N_Deer <- S_Deer + E_Deer + I_Deer
    births_in <- death_rate_S*S_Deer + hunt_kill_rate_S*S_Deer*N_Hunters + death_rate_E*E_Deer + hunt_kill_rate_EI*E_Deer*N_Hunters + death_rate_I*I_Deer + hunt_kill_rate_EI*I_Deer*N_Hunters
    d_S_Deer <- births_in - death_rate_S*S_Deer - beta_CWD_Deer*S_Deer*(I_Deer)/N_Deer - beta_CWD_Natural*S_Deer/N_Deer - hunt_kill_rate_S*S_Deer*N_Hunters
    d_E_Deer <- beta_CWD_Deer*S_Deer*(I_Deer)/N_Deer + beta_CWD_Natural*S_Deer/N_Deer - epsilon_Deer*E_Deer - death_rate_E*E_Deer - hunt_kill_rate_EI*E_Deer*N_Hunters
    d_I_Deer <- epsilon_Deer*E_Deer - death_rate_I*I_Deer - hunt_kill_rate_EI*I_Deer*N_Hunters
    return(list(c(d_S_Deer, d_E_Deer, d_I_Deer)))
  })
}

beta_CWD_Deer <- best_beta_CWD_Deer

ini_cwd_mod <- c(S_Deer = CWD_Results_data$yS_Deer[252], E_Deer = CWD_Results_data$yE_Deer[252], I_Deer = CWD_Results_data$yI_Deer[252])

times_cwd_mod <- seq(253,492,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate_S = hunter_kill_rate_month*1.15,
                  hunt_kill_rate_EI = hunter_kill_rate_month*1.35,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod_diff_hunt_rates, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results_scenario_1_projection <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results_scenario_1_projection$Prevalence <- CWD_Results_scenario_1_projection$yI_Deer / (CWD_Results_scenario_1_projection$yS_Deer +
                                                                                               CWD_Results_scenario_1_projection$yE_Deer +
                                                                                               CWD_Results_scenario_1_projection$yI_Deer)
CWD_Results_scenario_1_projection$Year <- seq((2019+(1/12)),2039,by=(1/12))

# Scenario 2 model continues 20 years
beta_CWD_Deer <- best_beta_CWD_Deer*1.3

ini_cwd_mod <- c(S_Deer = CWD_Results_data$yS_Deer[252], E_Deer = CWD_Results_data$yE_Deer[252], I_Deer = CWD_Results_data$yI_Deer[252])

times_cwd_mod <- seq(253,492,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate_S = hunter_kill_rate_month*1.15,
                  hunt_kill_rate_EI = hunter_kill_rate_month*1.35,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod_diff_hunt_rates, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results_scenario_2_projection <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results_scenario_2_projection$Prevalence <- CWD_Results_scenario_2_projection$yI_Deer / (CWD_Results_scenario_2_projection$yS_Deer +
                                                                                               CWD_Results_scenario_2_projection$yE_Deer +
                                                                                               CWD_Results_scenario_2_projection$yI_Deer)
CWD_Results_scenario_2_projection$Year <- seq((2019+(1/12)),2039,by=(1/12))

# Scenario 3 model continues 20 years
beta_CWD_Deer <- best_beta_CWD_Deer*1.3

ini_cwd_mod <- c(S_Deer = CWD_Results_data$yS_Deer[252], E_Deer = CWD_Results_data$yE_Deer[252], I_Deer = CWD_Results_data$yI_Deer[252])

times_cwd_mod <- seq(253,492,by = 1) #By Month

pars_cwd_mod <- c(births_in = births_in, death_rate_S = death_rate_S, 
                  beta_CWD_Deer = beta_CWD_Deer, 
                  beta_CWD_Natural = beta_CWD_Natural, 
                  hunt_kill_rate_S = hunter_kill_rate_month,
                  hunt_kill_rate_EI = hunter_kill_rate_month*1.7,
                  epsilon_Deer = epsilon_Deer, 
                  death_rate_E = death_rate_E,
                  death_rate_I = death_rate_I, 
                  N_Hunters = num_hunters_2017)

#Output
CWD_mod_Out <- ode(ini_cwd_mod, times_cwd_mod, CWD_mod_diff_hunt_rates, pars_cwd_mod)

yS_Deer <- CWD_mod_Out[,"S_Deer"]
yE_Deer <- CWD_mod_Out[,"E_Deer"]
yI_Deer <- CWD_mod_Out[,"I_Deer"]
CWD_Results_scenario_3_projection <- data.frame(yS_Deer, yE_Deer, yI_Deer)
CWD_Results_scenario_3_projection$Prevalence <- CWD_Results_scenario_3_projection$yI_Deer / (CWD_Results_scenario_3_projection$yS_Deer +
                                                                                               CWD_Results_scenario_3_projection$yE_Deer +
                                                                                               CWD_Results_scenario_3_projection$yI_Deer)
CWD_Results_scenario_3_projection$Year <- seq((2019+(1/12)),2039,by=(1/12))

ggplot()+
  geom_point(data = df_1999_2019_SouthernFarmland, aes(x = Year, y = Incidence, color = 'Southern Farmland Zone Data'), size = 3)+
  geom_line(data = CWD_Results_data, aes(x = Year, y = Prevalence, color = 'Calibrated Model'), size = 2)+
  geom_line(data = CWD_Results_baseline_projection, aes(x = Year, y = Prevalence, color = 'Baseline Projection'), size = 2)+
  geom_line(data = CWD_Results_scenario_1_projection, aes(x = Year, y = Prevalence, color = 'Scenario 1 Projection'), size = 2)+
  geom_line(data = CWD_Results_scenario_2_projection, aes(x = Year, y = Prevalence, color = 'Scenario 2 Projection'), size = 2)+
  geom_line(data = CWD_Results_scenario_3_projection, aes(x = Year, y = Prevalence, color = 'Scenario 3 Projection'), size = 2)+
  xlab('Year')+
  ylab('Prevalence')+
  ggtitle('Chronic Wasting Disease, Whitetail Deer, Wisconsin')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        title = element_text(size = 16),
        legend.text = element_text(size = 24),
        legend.position = c(0.2,0.8))+
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15),
                     labels = c('0%','5%','10%','15%'))+
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020,2025,2030,2035),
                     labels = c('2000','2005','2010','2015','2020','2025','2030','2035'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(color=' ', shape = '')
