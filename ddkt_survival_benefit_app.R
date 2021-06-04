# Shiny app for calculating the survival benefit of kidney transplantation

library(shiny)
library(shinythemes)
library(tidyverse)
library(coxme)
library(knitr)

to_model <- read_csv("clean_data.csv")
load("full_model_results.Rdata")

recips <- to_model %>% 
    filter(tx == 1)

i_time_cuts <- quantile(recips$ischemic_time, na.rm = TRUE, probs = c(0.25, 0.75))

centers <- c("mean center", unique(arrange(to_model, center)$center)) %>% na.omit()


kdpi_mapping_table <- read_csv("kdpi_mapping_table_2019.csv") %>%
    mutate(midpoint = (max-min)/2 + min,
           kdpi = row_number()-1) %>%
    filter(kdpi>0 & kdpi<100)


min_tx_kdri <- min(to_model$tx_kdri, na.rm = TRUE)

shift_kdri <- function(x, min_kdri =min_tx_kdri ) ifelse(x ==0, 0 ,x + abs(min_kdri) + 1)

KDRI_raw <- function(KDPI, scaling_factor = 1.28991045343964, mapping_table = kdpi_mapping_table){
    
    KDRI_median <- filter(mapping_table, kdpi == KDPI)$midpoint
    
    KDRI_rao <- KDRI_median*scaling_factor
    
    log(KDRI_rao)
}

# Define key functions to generate counterfactuals

log_hazard_w_tx <- function(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, model= full_me_model, icuts = i_time_cuts, ...){
    # 
    # if(between(kdpi, 1, 99) != TRUE){
    #   return("please enter an integer KDPI value between 1-99")
    # }
    
    no_dial <- ifelse(d_time == 0, 1,0)
    age_floor_25 <- ifelse(Age < 25, 0, Age - 25)
    log_dial_time <- log(d_time + 1)
    tx_kdri <- KDRI_raw(kdpi)
    #print(tx_kdri)
    shifted_tx_kdri <- shift_kdri(tx_kdri)
    
    
    i_time_low <- ifelse(ischemic_time <= icuts[[1]], 1, 0)
    i_time_high <- ifelse(ischemic_time > icuts[[2]], 1, 0)
    #print(shifted_tx_kdri)
    
    user_patient_params <- c(age_floor_25, log_dial_time, no_dial, diabetes, CAN_PREV_TX, i_time_low, i_time_high, shifted_tx_kdri)
    
    names(user_patient_params) <- c("age_floor_25", "log_dial_time", "no_dial", "diabetes", "CAN_PREV_TX",
                                    "i_time_low", "i_time_high", "shifted_tx_kdri")
    
    model <- full_me_model
    
    fixed_effects <- fixef(model)
    x_vars <- names(fixed_effects)
    ref <- data.frame(ranef(model)[[1]])
    ref$center <- row.names(ref)
    tx_ref <- filter(ref, center == ctr_id)$tx
    int_ref <- filter(ref, center == ctr_id)$Intercept
    
    #print(paste0("center intercept =", int_ref))
    #print(paste0("center tx effect  =", tx_ref))
    
    h_out <- tx_ref + unname(fixed_effects['tx']) +int_ref + 
        i_time_low*unname(fixed_effects['i_time_low'])+ 
        i_time_high*unname(fixed_effects['i_time_high'])+ 
        shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri']) +
        i_time_low*shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri:i_time_low'])+ 
        i_time_high*shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri:i_time_high'])
    
    #print(paste0("log hazard after ischemic time", h_out))
    
    
    for (x in names(user_patient_params)) {
        value <- unname(user_patient_params[x])
        #print(x)
        #print(value)
        if (!x %in% c("shifted_tx_kdri", "i_time_low", "i_time_high")){
            
            #print("added")
            int_var <- paste0(x, ":tx")
            kdri_int <- paste0(x, ":shifted_tx_kdri")
            
            
            # print(fixed_effects[x])
            # print(fixed_effects[int_var])
            # print(unname(fixed_effects[kdri_int]))
            
            add_hz <- value*(unname(fixed_effects[x]) + 
                                 unname(fixed_effects[int_var]) + 
                                 unname(fixed_effects[kdri_int])*shifted_tx_kdri)
            
            #print(add_hz)
            
            h_out <- h_out + add_hz
        } 
    }
    
    
    return(h_out)
    
}

log_hazard_waitlist <- function(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, model= full_me_model, ...){
    no_dial <- ifelse(d_time == 0, 1,0)
    age_floor_25 <- ifelse(Age < 25, 0, Age - 25)
    log_dial_time <- log(d_time + 1)
    
    user_patient_params <- c(age_floor_25, log_dial_time, no_dial, diabetes, CAN_PREV_TX)
    
    names(user_patient_params) <- c("age_floor_25", "log_dial_time", "no_dial", "diabetes", "CAN_PREV_TX")
    
    fixed_effects <- fixef(model)
    x_vars <- names(fixed_effects)
    ref <- data.frame(ranef(model)[[1]])
    ref$center <- row.names(ref)
    int_ref <- filter(ref, center == ctr_id)$Intercept
    
    h_out <- int_ref
    
    for (x in names(user_patient_params)) {
        value <- unname(user_patient_params[x])
        add_hz <- value*(unname(fixed_effects[x]))
        h_out <- h_out + add_hz
    }
    
    return(h_out)
    
}


percent <- function(x, digits = 1, format = "f", ...) {      
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}


survival_benefit <-  function(t_transplant, ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi,
                              output = "plot",  me_model = full_me_model, isch_cuts = i_time_cuts, b_hazard = full_base_hz){
    
    
    cum_hz <- filter(b_hazard, time == t_transplant)$hazard
    
    df <- b_hazard %>%
        mutate(log_hr_tx = log_hazard_w_tx(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                                           model = me_model, icuts = isch_cuts),
               log_hr_wait = log_hazard_waitlist(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, model = me_model),
               hazard = hazard - cum_hz) %>%
        filter(time >= t_transplant) %>%
        mutate(
            H_tx = hazard*exp(log_hr_tx),
            H_wait = hazard*exp(log_hr_wait),
            S_tx = exp(-H_tx),
            S_wait = exp(-H_wait)) %>%
        select(time, S_tx, S_wait) 
    
    
    if (output == "plot"){
        return(df %>%
                   filter(time < t_transplant + 365*5) %>%
                   mutate(years_from_transplant = (time - t_transplant)/365,
                          tx_dot = ifelse(years_from_transplant==0, 1, NA)) %>%
                   ggplot(aes(x = years_from_transplant)) +
                   geom_ribbon(aes(ymin = S_wait, ymax = S_tx, fill = "Survival Benefit of DDKT"), alpha = 0.75) + 
                   geom_step(aes(y = S_tx, color = "Post-Transplant"), size = 1.2) +
                   geom_step(aes(y = S_wait, color = "Waitlist"), size = 1.2) +
                   geom_point(aes(y= tx_dot, shape= "Transplant"), size = 5, colour = "darkgreen") +
                   scale_x_continuous(breaks = seq(0, 5)) +
                   scale_color_manual(name = "Survival function", values=c("blue", "red")) +
                   scale_fill_manual(name = NULL, values=c("skyblue")) +
                   ggthemes::theme_gdocs() +
                   #theme_minimal() +
                   lims(y = c(0,1)) + labs(y = "Survival", x = "Time (years after transplant)", shape = "") +
                   guides(shape = guide_legend(order = 1), color = guide_legend(order = 2), fill = guide_legend(order = 3))
        )
    }
    
    if (output == "rmst"){
        rmst_df <-df %>%
            mutate(days = (S_tx - S_wait)*(time-lag(time)),
                   days = ifelse(row_number() == 1, 0, days),
                   rmst = cumsum(days)) 
        return(
            paste0("DDKT would improve survival by ", round(filter(rmst_df, time == t_transplant + 5*365)$rmst), " days on average within the first 5 years following transplant")
        )
    }
    if (output == "table"){
        benefit_df <- df %>%
            filter(time %in% c(t_transplant + 365, t_transplant + 3*365,t_transplant + 5*365))
        
        
        return(
            benefit_df %>%
                mutate(S_tx = S_tx,
                       S_wait = S_wait,
                       surv_benefit = S_tx - S_wait,
                       years_from_transplant = (time - t_transplant)/365) %>%
                select("Years from transplant" = years_from_transplant, 
                       "Survival (transplant)" = S_tx, 
                       "Survival (waitlist)" = S_wait, 
                       "Absolute Survival Benefit" = surv_benefit) %>%
                mutate(across(c(2,3,4), percent)) 
        )
        
    }
    
}

log_hazard_w_tx_no_re <- function(Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, model= full_me_model, icuts = i_time_cuts, ...){
    # 
    # if(between(kdpi, 1, 99) != TRUE){
    #   return("please enter an integer KDPI value between 1-99")
    # }
    
    no_dial <- ifelse(d_time == 0, 1,0)
    age_floor_25 <- ifelse(Age < 25, 0, Age - 25)
    log_dial_time <- log(d_time + 1)
    tx_kdri <- KDRI_raw(kdpi)
    #print(tx_kdri)
    shifted_tx_kdri <- shift_kdri(tx_kdri)
    
    
    i_time_low <- ifelse(ischemic_time <= icuts[[1]], 1, 0)
    i_time_high <- ifelse(ischemic_time > icuts[[2]], 1, 0)
    #print(shifted_tx_kdri)
    
    user_patient_params <- c(age_floor_25, log_dial_time, no_dial, diabetes, CAN_PREV_TX, i_time_low, i_time_high, shifted_tx_kdri)
    
    names(user_patient_params) <- c("age_floor_25", "log_dial_time", "no_dial", "diabetes", "CAN_PREV_TX",
                                    "i_time_low", "i_time_high", "shifted_tx_kdri")
    
    model <- full_me_model
    
    fixed_effects <- fixef(model)
    x_vars <- names(fixed_effects)
    
    #print(paste0("center intercept =", int_ref))
    #print(paste0("center tx effect  =", tx_ref))
    
    h_out <- unname(fixed_effects['tx']) +
        i_time_low*unname(fixed_effects['i_time_low'])+ 
        i_time_high*unname(fixed_effects['i_time_high'])+ 
        shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri']) +
        i_time_low*shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri:i_time_low'])+ 
        i_time_high*shifted_tx_kdri*unname(fixed_effects['shifted_tx_kdri:i_time_high'])
    
    #print(paste0("log hazard after ischemic time", h_out))
    
    
    for (x in names(user_patient_params)) {
        value <- unname(user_patient_params[x])
        #print(x)
        #print(value)
        if (!x %in% c("shifted_tx_kdri", "i_time_low", "i_time_high")){
            
            #print("added")
            int_var <- paste0(x, ":tx")
            kdri_int <- paste0(x, ":shifted_tx_kdri")
            
            
            # print(fixed_effects[x])
            # print(fixed_effects[int_var])
            # print(unname(fixed_effects[kdri_int]))
            
            add_hz <- value*(unname(fixed_effects[x]) + 
                                 unname(fixed_effects[int_var]) + 
                                 unname(fixed_effects[kdri_int])*shifted_tx_kdri)
            
            #print(add_hz)
            
            h_out <- h_out + add_hz
        } 
    }
    
    
    return(h_out)
    
}





log_hazard_waitlist_no_re <- function(Age, d_time, diabetes, CAN_PREV_TX, model= full_me_model, ...){
    no_dial <- ifelse(d_time == 0, 1,0)
    age_floor_25 <- ifelse(Age < 25, 0, Age - 25)
    log_dial_time <- log(d_time + 1)
    
    user_patient_params <- c(age_floor_25, log_dial_time, no_dial, diabetes, CAN_PREV_TX)
    
    names(user_patient_params) <- c("age_floor_25", "log_dial_time", "no_dial", "diabetes", "CAN_PREV_TX")
    
    fixed_effects <- fixef(model)
    x_vars <- names(fixed_effects)
    
    h_out <- 0
    
    for (x in names(user_patient_params)) {
        value <- unname(user_patient_params[x])
        add_hz <- value*(unname(fixed_effects[x]))
        h_out <- h_out + add_hz
    }
    
    return(h_out)
    
}

survival_benefit_no_re <-  function(t_transplant, Age, d_time, diabetes, CAN_PREV_TX, 
                                    ischemic_time, kdpi,
                                    output = "plot", me_model = full_me_model, isch_cuts =  i_time_cuts, b_hazard = full_base_hz){
    
    
    cum_hz <- filter(b_hazard, time == t_transplant)$hazard
    
    df <- b_hazard %>%
        mutate(log_hr_tx = log_hazard_w_tx_no_re(Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                                                 model = me_model, icuts = isch_cuts),
               log_hr_wait = log_hazard_waitlist_no_re(Age, d_time, diabetes, CAN_PREV_TX, model = me_model),
               hazard = hazard - cum_hz) %>%
        filter(time >= t_transplant) %>%
        mutate(
            H_tx = hazard*exp(log_hr_tx),
            H_wait = hazard*exp(log_hr_wait),
            S_tx = exp(-H_tx),
            S_wait = exp(-H_wait)) %>%
        select(time, S_tx, S_wait) 
    
    if (output == "plot"){
        return(df %>%
                   filter(time < t_transplant + 365*5) %>%
                   mutate(years_from_transplant = (time - t_transplant)/365,
                          tx_dot = ifelse(years_from_transplant==0, 1, NA)) %>%
                   ggplot(aes(x = years_from_transplant)) +
                   geom_ribbon(aes(ymin = S_wait, ymax = S_tx, fill = "Survival Benefit of DDKT"), alpha = 0.75) + 
                   geom_step(aes(y = S_tx, color = "Post-Transplant"), size = 1.2) +
                   geom_step(aes(y = S_wait, color = "Waitlist"), size = 1.2) +
                   geom_point(aes(y= tx_dot, shape= "Transplant"), size = 5, colour = "darkgreen") +
                   scale_x_continuous(breaks = seq(0, 5)) +
                   scale_color_manual(name = "Survival function", values=c("blue", "red")) +
                   scale_fill_manual(name = NULL, values=c("skyblue")) +
                   ggthemes::theme_gdocs() +
                   #theme_minimal() +
                   lims(y = c(0,1)) + labs(y = "Survival", x = "Time (years after transplant)", shape = "") +
                   guides(shape = guide_legend(order = 1), color = guide_legend(order = 2), fill = guide_legend(order = 3))
        )
    }
    
    if (output == "rmst"){
        rmst_df <-df %>%
            mutate(days = (S_tx - S_wait)*(time-lag(time)),
                   days = ifelse(row_number() == 1, 0, days),
                   rmst = cumsum(days)) 
        return(
            paste0("DDKT would improve survival by ", round(filter(rmst_df, time == t_transplant + 5*365)$rmst), " days on average within the first 5 years following transplant")
        )
    }
    
    if (output == "table"){
        benefit_df <- df %>%
            filter(time %in% c(t_transplant + 365, t_transplant + 3*365,t_transplant + 5*365))
        
        
        return(
            benefit_df %>%
                mutate(S_tx = S_tx,
                       S_wait = S_wait,
                       surv_benefit = S_tx - S_wait,
                       years_from_transplant = (time - t_transplant)/365) %>%
                select("Years from transplant" = years_from_transplant, 
                       "Survival (transplant)" = S_tx, 
                       "Survival (waitlist)" = S_wait, 
                       "Absolute Survival Benefit" = surv_benefit) %>%
                mutate(across(c(2,3,4), percent))
        )
        
    }
    
    
    
}


#survival_benefit(t_transplant = 500, 782, 79, 2, 1, 0, 15, 25, plot = TRUE)

#survival_benefit_no_re(t_transplant = 500, 79, 2, 1, 0, 15, 25, plot = TRUE)



# Define UI 
ui <- fluidPage(
    theme = shinytheme("cerulean"),
    # Application title
    titlePanel("Survival Benefit of Deceased Donor Kidney Transplantation"),

    # Sidebar for inputs into survival benefit function
    sidebarLayout(
        sidebarPanel(
            sliderInput("Age",
                        "Candidate Age:",
                        min = 18,
                        max = 80,
                        value = 55),
            sliderInput("t_transplant", "Days on waitlist", 800, min = 1, max = 3000),
            sliderInput("d_time", "Dialysis Time (years)", 3.5, min = 0, max = 10, step = 0.5),
            selectInput("diabetes", "Diabetes", c("No", "Yes")),
            selectInput("CAN_PREV_TX", "History of previous solid organ transplant", c("No", "Yes")),
            sliderInput("kdpi",
                        "KDPI (%):",
                        min = 1,
                        max = 99,
                        value = 50),
            selectInput("ischemic_time",
                        "Cold ischemic time (hours):",
                        c("<12", "12-22", ">22"),
                        selected = "12-22"),
            selectInput("ctr_id", "Transplant Center", centers, selected = "mean center")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("survPlot"),
           tableOutput("survTable"),
           textOutput("rmstText"),
           h4("Survival with and without transplant estimated from all adult candidates listed for deceased donor kidney-alone transplanation, 2005-2010")
        )
    ),
    h4("William F Parker, MD, MS, Yolanda Becker MD, Robert Gibbons, PhD"),
    h5("A commissioned paper for National Academies of Sciences, Engineering, and Medicine study committee: A Fairer and More Equitable Cost-Effective and Transparent System of Donor Organ Procurement, Allocation, and Distribution")
)

# Define server logic required to calculate survival benefit
server <- function(input, output, session) {
    

    observeEvent(input$t_transplant,{
        
        cur_d_time <- isolate(input$d_time)
        
        update <- ifelse(cur_d_time > input$t_transplant/365, cur_d_time, input$t_transplant/365)
        updateSliderInput(session, "d_time", value = (update))
        
    })
    
    output$survPlot <- renderPlot({
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time <- case_when(
            input$ischemic_time == "<12" ~ 1,
            input$ischemic_time == ">22" ~ 30,
            TRUE ~ 16
        )
        
        if (input$ctr_id == "mean center"){
            survival_benefit_no_re(input$t_transplant, 
                                   input$Age, input$d_time, diabetes, CAN_PREV_TX,
                                   i_time,
                                   input$kdpi)
            
        } else {
        cur_ctr <- as.numeric(input$ctr_id)
        # model, b_hazard, t_transplant, ctr_id, Age, d_time, diabetes, CAN_PREV_TX, kdpi
        survival_benefit(input$t_transplant, cur_ctr, 
                         input$Age, input$d_time, diabetes, CAN_PREV_TX,
                         i_time,
                         input$kdpi)
        }
    })
    
    
    
    output$survTable <- renderTable({
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time <- case_when(
            input$ischemic_time == "<12" ~ 1,
            input$ischemic_time == ">22" ~ 30,
            TRUE ~ 16
        )
        
        # model, b_hazard, t_transplant, ctr_id, Age, d_time, diabetes, CAN_PREV_TX, kdpi
        if (input$ctr_id == "mean center"){
        survival_benefit_no_re(input$t_transplant,
                         input$Age, input$d_time, diabetes, CAN_PREV_TX,
                         i_time,
                         input$kdpi,
                         output = "table")
        } else {
            cur_ctr <- as.numeric(input$ctr_id)
            survival_benefit(input$t_transplant, cur_ctr,  
                             input$Age, input$d_time, diabetes, CAN_PREV_TX,
                             i_time,
                             input$kdpi,
                             output = "table")
            
        }
    }, digits =  0, align = "c")
    
    output$rmstText <- renderText({
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time <- case_when(
            input$ischemic_time == "<12" ~ 1,
            input$ischemic_time == ">22" ~ 30,
            TRUE ~ 16
        )
        
        # model, b_hazard, t_transplant, ctr_id, Age, d_time, diabetes, CAN_PREV_TX, kdpi
        if (input$ctr_id == "mean center"){
            survival_benefit_no_re(input$t_transplant,
                                   input$Age, input$d_time, diabetes, CAN_PREV_TX,
                                   i_time,
                                   input$kdpi,
                                   output = "rmst")
        } else {
            cur_ctr <- as.numeric(input$ctr_id)
            survival_benefit(input$t_transplant, cur_ctr,  
                             input$Age, input$d_time, diabetes, CAN_PREV_TX,
                             i_time,
                             input$kdpi,
                             output = "rmst")
            
        }
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
