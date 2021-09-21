# Shiny app for calculating the survival benefit of kidney transplantation

library(shiny)
library(shinythemes)
library(tidyverse)
library(coxme)
library(knitr)

load("app_data.Rdata")

kdpi_mapping_table <- read_csv("kdpi_mapping_table_2019.csv") %>%
    mutate(midpoint = (max-min)/2 + min,
           kdpi = row_number()-1) %>%
    filter(kdpi>0 & kdpi<100)

KDRI_raw <- function(KDPI, scaling_factor = 1.28991045343964, mapping_table = kdpi_mapping_table){
    
    KDRI_median <- filter(mapping_table, kdpi == KDPI)$midpoint
    
    KDRI_rao <- KDRI_median*scaling_factor
    
    log(KDRI_rao)
}

# Define key functions to generate counterfactuals
sum_hazard <- function(param_list,start_hazard, covariate_values, model_coef){
    h_out <- start_hazard
    
    for (wp in param_list ) {
        
        params <- str_split(wp, ":")
        value <- prod(sapply(params[[1]], function(x) unname(covariate_values[x])))
        #print(value)
        beta <- unname(model_coef[wp])
        #print(beta)
        add_hz <- value*beta
        #print(add_hz)
        h_out <- h_out + add_hz
        
    }
    return(h_out)
}


log_hazard <- function(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, 
                       ischemic_time, kdpi,
                       period = "wait",
                       model= full_me_model,
                       icuts = i_time_cuts,
                       ...){
    
    no_dial <- ifelse(d_time == 0, 1,0)
    age_floor_25 <- ifelse(Age < 25, 0, Age - 25)
    log_dial_time <- log(d_time + 1)
    tx_kdri <- KDRI_raw(kdpi)
    
    i_time_low <- ifelse(ischemic_time <= icuts[[1]], 1, 0)
    i_time_high <- ifelse(ischemic_time > icuts[[2]], 1, 0)
    
    user_patient_params <- c(age_floor_25, log_dial_time, no_dial, diabetes, CAN_PREV_TX, 
                             i_time_low, i_time_high, tx_kdri)
    
    names(user_patient_params) <- c("age_floor_25", "log_dial_time", "no_dial", "diabetes", "CAN_PREV_TX",
                                    "i_time_low", "i_time_high", "tx_kdri")
    
    user_patient_params[["i_time_unknown"]] <- 0
    
    fixed_effects <- fixef(model)
    x_vars <- names(fixed_effects)
    ref <- data.frame(ranef(model)[[1]])
    ref$center <- row.names(ref)
    
    if (ctr_id == "mean center"){
        tx_ref <- 0
        int_ref <- 0
    } else {
        tx_ref <- filter(ref, center == ctr_id)$tx
        int_ref <- filter(ref, center == ctr_id)$Intercept
    }
    
    if (period == "wait"){
        entered_params <- c(user_patient_params, "tx" =0, "immediate_post_trans"=0, "second_post_trans_int" =0)
        entered_params[["tx_kdri"]] <- 0
        entered_params[["i_time_low"]] <- 0
        entered_params[["i_time_high"]] <- 0
        st_hazard <- int_ref
    }
    
    if (period == "post_1"){
        entered_params <- c(user_patient_params, "tx" =1, "immediate_post_trans"=1, "second_post_trans_int" =0)
        st_hazard <- int_ref + tx_ref
    }
    
    if (period == "post_2"){
        entered_params <- c(user_patient_params, "tx" =1, "immediate_post_trans"=0, "second_post_trans_int" =1)
        st_hazard <- int_ref + tx_ref
    }
    
    if (period == "post_3"){
        entered_params <- c(user_patient_params, "tx" =1, "immediate_post_trans"=0, "second_post_trans_int" =0)
        st_hazard <- int_ref + tx_ref
    }
    # print(x_vars)
    # print(st_hazard)
    # print(entered_params)
    # print(fixed_effects)
    h_out <- sum_hazard(x_vars, st_hazard, entered_params, fixed_effects)
    return(h_out)
    
}


percent <- function(x, digits = 1, format = "f", ...) {      
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

survival_benefit <-  function(t_transplant, 
                              ctr_id, 
                              Age, d_time, diabetes, CAN_PREV_TX, 
                              ischemic_time, kdpi,
                              output = "plot",  
                              me_model = full_me_model, 
                              isch_cuts = i_time_cuts,
                              b_hazard = full_base_hz,
                              post_tx_int_1 = 30,
                              post_tx_int_2 = 180){
    
    log_hr_w <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                           model = me_model, icuts = isch_cuts, period = "wait")
    
    log_hr_post_1 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                                model = me_model, icuts = isch_cuts, period = "post_1") 
    
    log_hr_post_2 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                                model = me_model, icuts = isch_cuts, period = "post_2")              
    
    log_hr_post_3 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, ischemic_time, kdpi, 
                                model = me_model, icuts = isch_cuts, period = "post_3")    
    
    cum_hz <- filter(b_hazard, time == t_transplant)$hazard
    
    df <- b_hazard %>%
        mutate(hazard = hazard- lag(hazard),
               hazard = ifelse(is.na(hazard), 0, hazard)) %>%
        filter(time >= t_transplant) %>%
        mutate(
            time = time - t_transplant,
            hazard_wait = hazard*exp(log_hr_w),
            hazard_tx = case_when(
                time < post_tx_int_1 ~ hazard*exp(log_hr_post_1),
                time < post_tx_int_2 ~ hazard*exp(log_hr_post_2),
                TRUE ~ hazard*exp(log_hr_post_3)
            ),
            H_tx = cumsum(hazard_tx),
            H_wait = cumsum(hazard_wait),
            S_tx = exp(-H_tx),
            S_wait = exp(-H_wait),
            S_tx_better = ifelse(S_tx >= S_wait, S_tx, NA),
            S_tx_worse = ifelse(S_tx <= S_wait, S_tx, NA),
            S_wait_better = ifelse(S_wait >= S_tx, S_wait, NA),
            S_wait_worse = ifelse(S_wait <= S_tx, S_wait, NA))
    
    if (output == "dataframe"){
        return(df)
    }
    
    if (output == "plot"){
        return(df %>%
                   filter(time < 365*5) %>%
                   mutate(years_from_transplant = time/365,
                          tx_dot = ifelse(years_from_transplant==0, 1, NA)) %>%
                   ggplot(aes(x = years_from_transplant)) +
                   geom_ribbon(aes(ymin = S_wait_worse, ymax = S_tx_better, fill = "Survival Benefit"), alpha = 0.75) +
                   geom_step(aes(y = S_tx, color = "Post-Transplant"), size = 1.1) +
                   geom_step(aes(y = S_wait, color = "Waitlist"), size = 1.1) +
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
            paste0("DDKT would improve survival by ", round(filter(rmst_df, time == 5*365)$rmst/365, digits = 1), " years on average within the first 5 years following transplant")
        )
    }
    if (output == "table"){
        benefit_df <- df %>%
            filter(time %in% c( 365, 3*365, 5*365))
        
        
        return(
            benefit_df %>%
                mutate(S_tx = S_tx,
                       S_wait = S_wait,
                       surv_benefit = S_tx - S_wait,
                       years_from_transplant = time/365) %>%
                select("Years from transplant" = years_from_transplant,
                       "Survival (transplant)" = S_tx,
                       "Survival (waitlist)" = S_wait,
                       "Absolute Survival Benefit" = surv_benefit) %>%
                mutate(across(c(2,3,4), percent)) %>%
                kable()
        )
        
    }
    
}


wait_for_better <- function(t_a, kdpi_a, i_time_a,
                            t_b, kdpi_b, i_time_b,
                            Age, d_time, diabetes, CAN_PREV_TX, ctr_id = "mean center", 
                            output = "plot",  me_model = full_me_model, 
                            isch_cuts = i_time_cuts, b_hazard = full_base_hz,
                            post_tx_int_1 = 30, post_tx_int_2 = 180){
    
    df <- survival_benefit(t_a, ctr_id, Age, d_time, diabetes, CAN_PREV_TX, i_time_a, kdpi_a,
                           output = "dataframe") %>%
        select(time, hazard, hazard_wait, S_a = S_tx, S_wait)
    
    
    log_hr_post_1 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, i_time_b, kdpi_b, 
                                model = me_model, icuts = isch_cuts, period = "post_1") 
    
    log_hr_post_2 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, i_time_b, kdpi_b, 
                                model = me_model, icuts = isch_cuts, period = "post_2")              
    
    log_hr_post_3 <- log_hazard(ctr_id, Age, d_time, diabetes, CAN_PREV_TX, i_time_b, kdpi_b, 
                                model = me_model, icuts = isch_cuts, period = "post_3") 
    
    t_b <- t_b - t_a
    
    df <- df %>%
        mutate(hazard_b = case_when(
            time < t_b ~ hazard_wait,
            time < t_b + post_tx_int_1 ~ hazard*exp(log_hr_post_1),
            time < t_b + post_tx_int_2 ~ hazard*exp(log_hr_post_2),
            TRUE ~ hazard*exp(log_hr_post_3)
        ),
        H_b = cumsum(hazard_b),
        S_b = exp(-H_b),
        years_from_first_offer = time/365,
        kidney_a = case_when(years_from_first_offer == 0 ~ S_a),
        kidney_b = case_when(time == t_b ~ S_b),
        transplant = case_when(
            #years_from_first_offer == 0 ~ "First Offer",
            time == t_b ~ "Receive Kidney B"
        ),
        transplant_surv = case_when(
            #years_from_first_offer == 0 ~ S_a,
            time == t_b ~ S_b
        ),
        s_a_better = ifelse(S_a >= S_b, S_a, NA),
        s_a_worse = ifelse(S_b >= S_a, S_a, NA),
        s_b_better = ifelse(S_b >= S_a, S_b, NA),
        s_b_worse = ifelse(S_a >= S_b, S_b, NA)
        ) %>%
        filter(time < t_b + 365*5)
    
    if (all(is.na(df$s_a_worse))) {
        p <- df %>% 
            ggplot(aes(x =years_from_first_offer)) +
            geom_ribbon(aes(ymax = s_a_better, ymin = s_b_worse, fill = "Benefit of A"), alpha = 0.75) 
        
    } else if (all(is.na(df$s_b_worse))) {
        p <- df %>% 
            ggplot(aes(x =years_from_first_offer)) +
            geom_ribbon(aes(ymax = s_b_better, ymin = s_a_worse, fill = "Benefit of B"), alpha = 0.75)
        
    } else {
        p <- df %>% 
            ggplot(aes(x =years_from_first_offer)) +
            geom_ribbon(aes(ymax = s_a_better, ymin = s_b_worse, fill = "Benefit of A"), alpha = 0.75) + 
            geom_ribbon(aes(ymax = s_b_better, ymin = s_a_worse, fill = "Benefit of B"), alpha = 0.75)
        
    }    
    
    if (all(is.na(df$s_b_worse))) {
        p <- p +            
            lims(y = c(0, 1)) +
            ggthemes::theme_gdocs() +
            scale_x_continuous(breaks = seq(0,10, 1)) + 
            scale_color_manual(values = c("darkgreen", "blue", "red")) + 
            scale_fill_manual(values = c("skyblue")) + 
            labs(x = "Years from first offer", 
                 y = "Survival",
                 shape = "Transplant",
                 color = "Scenario",
                 fill = "")
    } else {
        p <- p +            
            lims(y = c(0, 1)) +
            ggthemes::theme_gdocs() +
            scale_x_continuous(breaks = seq(0,10, 1)) + 
            scale_color_manual(values = c("darkgreen", "blue", "red")) + 
            scale_fill_manual(values = c("forestgreen", "skyblue")) + 
            labs(x = "Years from first offer", 
                 y = "Survival",
                 shape = "Transplant",
                 color = "Scenario",
                 fill = "")
    }

    
    if (output == "plot") {
        
        out <- p +        
            geom_step(aes(y = S_a, color = "Take A now"), size = 0.7) +
            geom_step(aes(y = S_b, color = "Wait for B"), size = 0.7) +
            geom_point(aes(y = kidney_a, shape ="Kidney A"), size = 3) + 
            geom_point(aes(y = kidney_b, shape ="Kidney B"), size = 3) +
            ggrepel::geom_label_repel(aes(y = transplant_surv, label = transplant),
                                      min.segment.length = 0,
                                      nudge_y = -0.15)
    } 
    
    if (output == "wait_plot"){
        
        out <- p + 
            geom_step(aes(y = S_a, color = "Take kidney A now"), size = 0.7) +
            geom_step(aes(y = S_wait, color = "Waitlist"), size = 0.7)+
            geom_step(aes(y = S_b, color = "Wait for kidney B"), size = 0.7) +
            geom_point(aes(y = kidney_a, shape ="Kidney A"), size = 3) + 
            geom_point(aes(y = kidney_b, shape ="Kidney B"), size = 3) +
            ggrepel::geom_label_repel(aes(y = transplant_surv, label = transplant),
                                      min.segment.length = 0,
                                      nudge_y = -0.15)
        
    }
    
    
    if (output == "rmst"){
        
        df <- df %>%
            mutate(surv_diff = S_a - S_b)
        
        prob_death_pre_b <- (1- filter(df, time == t_b)$S_b)
        
        
        rmst_df <-df %>%
            mutate(days = (S_a - S_b)*(time-lag(time)),
                   days = ifelse(row_number() == 1, 0, days),
                   rmst = cumsum(days)) 
        
        days_a <- round(filter(rmst_df, time == 5*365)$rmst, digits = 1)
        
        if (days_a < 0){
            out <-  paste0("Waiting ", 
                           (t_b), " days for Kidney B would increase survival by ", 
                           (-days_a), " days on average within the first 5 years following transplant. The probability the patient dies while waiting is for Kidney B is ", 
                           round(100*prob_death_pre_b), "%.")
        } else {
            out <-  paste0("Accepting Kidney A now would increase survival by ", 
                           days_a, " days on average within the first 5 years following transplant compared to waiting ",(t_b), " days for Kidney B. The probability the patient dies while waiting for Kidney B is ", round(100*prob_death_pre_b), "%.")
            
        }
        
    }
    
    if (output == "table"){
        benefit_df <- df %>%
            filter(time %in% c(365, 3*365, 5*365))
        
        
        out <- benefit_df %>%
            mutate(surv_benefit = S_a - S_b,
                   sb_a = S_a - S_wait,
                   sb_b = S_b - S_wait,
                   years_from_first_offer = time/365) %>%
            select("Years from first offer" = years_from_first_offer , 
                   "Survival (accept A)" = S_a,
                   "Survival (wait for B)" = S_b,
                   "Survival (waitlist)" = S_wait, 
                   "Survival Benefit of Accepting A compared to waiting forever" = sb_a,
                   "Survival Benefit of wating for B compared to waiting forever" = sb_b,
                   "Survival Benefit of Accepting A over waiting for B" = surv_benefit) %>%
            mutate(across(c(2,3,4,5,6,7), percent)) 
        
    }
    
    
    
    return(out)
}

# Define UI 
ui <- fluidPage(
    theme = shinytheme("cerulean"),
    # Application title
    titlePanel("Take the offer or wait for a better one?"),
    h4("Comparing survival under different accept-reject strategies in Deceased Donor Kidney Transplantation"),

    # Sidebar for inputs into survival benefit function
    sidebarLayout(
        sidebarPanel(
            h4("Candidate Characterstics (at time of first offer)"),
            sliderInput("Age",
                        "Age:",
                        min = 18,
                        max = 80,
                        value = 55),
            sliderInput("t_transplant", "Days on waitlist", 800, min = 1, max = 3000),
            sliderInput("d_time", "Dialysis Time (years)", 3.5, min = 0, max = 10, step = 0.5),
            selectInput("diabetes", "Diabetes", c("No", "Yes")),
            selectInput("CAN_PREV_TX", "History of previous solid organ transplant", c("No", "Yes")),
            selectInput("ctr_id", "Transplant Center", centers, selected = "mean center"),
            checkboxInput("wait_plot", "Show waitlist survival?", value = TRUE),
            #h4("Survival with accepting kidney A now compared to waiting for kidney B"),
            h6("Model estimated from all adult candidates listed for deceased donor kidney-alone transplanation, 2005-2010")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("survPlot"),
           tags$head(tags$style("#survText{color: black;
                                 font-size: 16px;
                                 font-style: bold;
                                 }"
           )
           ),
           fluidRow(column(4,textOutput("survText")),
                   column(8,tableOutput("survTable"))),
           fluidRow(
               column(6,
                      wellPanel(
                          h5("Kidney A (current offer)"),
                          sliderInput("kdpi_a",
                                      "KDPI (%)",
                                      min = 1,
                                      max = 99,
                                      value = 80),
                          selectInput("ischemic_time_a",
                                      "Cold ischemic time (hours)",
                                      c("<12", "12-22", ">22"),
                                      selected = "12-22")
                      )),
               column(6,
                      wellPanel(
                          h5("Kidney B (future offer)"),
                          sliderInput("t_delta", "Days waiting for Kidney B", 150, min = 1, max = 500),
                          sliderInput("kdpi_b",
                                      "KDPI (%)",
                                      min = 1,
                                      max = 99,
                                      value = 40),
                          selectInput("ischemic_time_b",
                                      "Cold ischemic time (hours)",
                                      c("<12", "12-22", ">22"),
                                      selected = "12-22")
                      ))
           ),
           # tableOutput("survTable"),
           # textOutput("rmstText"),
        )
    ),
    h4("William F Parker, MD, MS, Yolanda Becker MD, Robert Gibbons, PhD"),
    h5("A commissioned paper for National Academies of Sciences, Engineering, and Medicine study committee: A Fairer and More Equitable Cost-Effective and Transparent System of Donor Organ Procurement, Allocation, and Distribution")
)

# Define server logic required to calculate survival benefit
server <- function(input, output, session) {
    
    prev_t_transplant <- isolate(input$t_transplant)
    
    observeEvent(input$t_transplant,{
        
        cur_d_time <- isolate(input$d_time)
        
        cur_age <- isolate(input$Age)
        update_years <- (input$t_transplant - prev_t_transplant)/365
        
        
        #update <- ifelse(cur_d_time > input$t_transplant/365, cur_d_time, input$t_transplant/365)
        updateSliderInput(session, "d_time", value = cur_d_time + update_years)
        
        updateSliderInput(session, "Age", value = cur_age + update_years)
    })
    

    
    output$survPlot <- renderPlot({
        
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time_a <- case_when(
            input$ischemic_time_a == "<12" ~ 1,
            input$ischemic_time_a == ">22" ~ 30,
            TRUE ~ 16
        )
        
        i_time_b <- case_when(
            input$ischemic_time_b == "<12" ~ 1,
            input$ischemic_time_b == ">22" ~ 30,
            TRUE ~ 16
        )
        
        t_b <- input$t_transplant + input$t_delta
        
        if (input$wait_plot == FALSE){
            wait_for_better(input$t_transplant, input$kdpi_a, i_time_a,
                            t_b, input$kdpi_b, i_time_b,
                            input$Age, input$d_time, diabetes, CAN_PREV_TX, input$ctr_id)
        } else {
            wait_for_better(input$t_transplant, input$kdpi_a, i_time_a,
                            t_b, input$kdpi_b, i_time_b,
                            input$Age, input$d_time, diabetes, CAN_PREV_TX,input$ctr_id,
                            output = "wait_plot")
        }
        
    })
    
    output$survText <- renderText({
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time_a <- case_when(
            input$ischemic_time_a == "<12" ~ 1,
            input$ischemic_time_a == ">22" ~ 30,
            TRUE ~ 16
        )
        
        i_time_b <- case_when(
            input$ischemic_time_b == "<12" ~ 1,
            input$ischemic_time_b == ">22" ~ 30,
            TRUE ~ 16
        )
        
        t_b <- input$t_transplant + input$t_delta
        
        wait_for_better(input$t_transplant, input$kdpi_a, i_time_a,
                        t_b, input$kdpi_b, i_time_b,
                        input$Age, input$d_time, diabetes, CAN_PREV_TX, input$ctr_id,
                        output = "rmst")
    })
    
    
    
    output$survTable <- renderTable({
        diabetes <- ifelse(input$diabetes == "Yes", 1, 0)
        CAN_PREV_TX <- ifelse(input$CAN_PREV_TX == "Yes", 1, 0)
        
        i_time_a <- case_when(
            input$ischemic_time_a == "<12" ~ 1,
            input$ischemic_time_a == ">22" ~ 30,
            TRUE ~ 16
        )
        
        i_time_b <- case_when(
            input$ischemic_time_b == "<12" ~ 1,
            input$ischemic_time_b == ">22" ~ 30,
            TRUE ~ 16
        )
        
        t_b <- input$t_transplant + input$t_delta
        
        wait_for_better(input$t_transplant, input$kdpi_a, i_time_a,
                        t_b, input$kdpi_b, i_time_b,
                        input$Age, input$d_time, diabetes, CAN_PREV_TX, input$ctr_id,
                        output = "table") %>% 
            select("Years from first offer","Survival (accept A)",
                   "Survival (wait for B)", 
                   "Survival Benefit of Accepting A over waiting for B")

    }, digits =  0, align = "c")
}

# Run the application 
shinyApp(ui = ui, server = server)
