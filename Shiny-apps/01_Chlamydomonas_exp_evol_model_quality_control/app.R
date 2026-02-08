# J R Laurich

# November 13 3, 2025

# Let's look at my models and raw data together

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(nls.multstart)

# Functions -------------------------------------------------

pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

empty_na_plot <- function(xlim = c(0, 100), ylim = c(-0.1, 2)) {
  ggplot() +
    theme_classic() +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    labs(x = "Days", y = "RFU") +
    annotate("text",
             x = mean(xlim), y = mean(ylim),
             label = "N/A", size = 6, fontface = "bold")
}

# Load & explore the data --------------------------------------------------------

###### Summary data ######

df.sum <- read.csv('data/27_summary_table.csv') # summary data

pop_order <- unique(df.sum$population) # Unique values for each population

pop_order<- c(
  as.character(sort(as.numeric(pop_order[grepl("^[0-9]+$", pop_order)]))),
  sort(pop_order[!grepl("^[0-9]+$", pop_order)])
) # re-arrange them to increase numerically, then alphabetically. 

df.sum <- df.sum %>%
  arrange(factor(population, levels = pop_order), rep.ID) # Now this is ordered by population and well ID

df.sum$rep <- rep (1:4, 37) # Add replicate numbers

rep_levels <- 1:4

###### µ estimates ######

df.l.mu <- read.csv('data/07_µ_estimates_light.csv') # Growth data, light
head(df.l.mu)

df.l.mu <- df.l.mu %>%
  arrange(factor(population, levels = pop_order), light, well.ID) # Now this is ordered by population, light level, and well ID
  
df.l.mu$rep <- rep (1:4, 370) # Add replicate numbers

df.n.mu <- read.csv('data/11_µ_estimates_nitrogen.csv') # Growth data, nitrogen
head(df.n.mu)

df.n.mu <- df.n.mu %>%
  arrange(factor(population, levels = pop_order), nit, well.ID) # Now this is ordered by population, nit level, and well ID

df.n.mu$rep <- rep (1:4, 370) # Add replicate numbers

df.p.mu <- read.csv('data/15_µ_estimates_phosphorous.csv') # Growth data, phosphorous
head(df.p.mu)

df.p.mu <- df.p.mu %>%
  arrange(factor(population, levels = pop_order), phos, well.ID) # Now this is ordered by population, phos level, and well ID

df.p.mu$rep <- rep (1:4, 370) # Add replicate numbers

df.s.mu <- read.csv('data/19_µ_estimates_salt.csv') # Growth data, salt
head(df.s.mu)

df.s.mu <- df.s.mu %>%
  mutate(population = case_when(
    population == "Anc 2" ~ "anc2",
    population == "Anc 3" ~ "anc3",
    population == "Anc 4" ~ "anc4",
    population == "Anc 5" ~ "anc5",
    TRUE ~ population
  ))

df.s.mu <- df.s.mu %>%
  arrange(factor(population, levels = pop_order), salt) # Now this is ordered by population, salt level, and well ID

df.t.mu <- read.csv('data/02_µ_estimates_temp.csv') # Growth data, temperature
head(df.t.mu)

df.t.mu <- df.t.mu %>%
  arrange(factor(population, levels = pop_order), temp, well.ID) # Now this is ordered by population, temp level, and well ID

df.t.mu$rep <- rep (1:4, 222) # Add replicate numbers

###### Raw growth rate time series data ######

df.l <- read.csv('data/06_light_rfus_time.csv') # Raw growth data, light
head(df.l)

df.l$percentage <- as.numeric(df.l$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df.l$percentage[is.na(df.l$percentage)] <- 0.6 # For now, let's set this to 0.6, but I need to talk with Joey about this. 

df.l <- df.l %>% 
  rename(well.ID = well_plate) %>% 
  mutate(level = percentage * 2.5)

df.l <- df.l %>% 
  select(well.ID, RFU, population, days, level) %>% 
  arrange(factor(population, levels = pop_order), well.ID)

rep <- numeric()

for (i in unique(df.l$population)){
  
  n <- nrow(df.l[df.l$population == i,])/4
  reps <- c(rep(1,n), rep(2,n), rep(3,n), rep(4,n))
  rep<- c (rep, reps)
  
}

df.l <- cbind(df.l, rep)
df.l$grad <- "Light"
df.l$unit <- "µM photons m^-2s^-1"

df.n <- read.csv('data/10_nitrogen_rfus_time.csv') # Raw growth data, nitrogen
head(df.n)

df.n <- df.n %>% 
  rename(well.ID = well_plate,
         level = nitrate_concentration) %>% 
  select(well.ID, RFU, population, days, level) %>% 
  arrange(factor(population, levels = pop_order), well.ID)

rep <- numeric()

for (i in unique(df.n$population)){
  
  n <- nrow(df.n[df.n$population == i,])/4
  reps <- c(rep(1,n), rep(2,n), rep(3,n), rep(4,n))
  rep<- c (rep, reps)
  
}

df.n <- cbind(df.n, rep)
df.n$grad <- "Nitrogen"
df.n$unit <- "µM nitrate"

df.p <- read.csv('data/14_phosphorous_rfus_time.csv') # Raw growth data, phosphorous
head(df.p)

df.p <- df.p %>% 
  rename(well.ID = well_plate,
         level = phosphate_concentration) %>% 
  select(well.ID, RFU, population, days, level) %>% 
  arrange(factor(population, levels = pop_order), well.ID)

rep <- numeric()

for (i in unique(df.p$population)){
  
  n <- nrow(df.p[df.p$population == i,])/4
  reps <- c(rep(1,n), rep(2,n), rep(3,n), rep(4,n))
  rep<- c (rep, reps)
  
}

df.p <- cbind(df.p, rep)
df.p$grad <- "Phosphorous"
df.p$unit <- "µM phosphate"

df.s <- read.csv('data/18_salt_rfus_time.csv') # Raw growth data, salt
head(df.s)

df.s$salt_level <- factor(df.s$treatment, levels = sort(unique(df.s$treatment)), ordered = TRUE) # Keep the numerical sorting.
df.s$salt <- as.numeric(df.s$salt_level) - 1# Subtract 1 to make it run from 0 to 9

df.s <- df.s %>%
  mutate(population = case_when(
    population == "Anc 2" ~ "anc2",
    population == "Anc 3" ~ "anc3",
    population == "Anc 4" ~ "anc4",
    population == "Anc 5" ~ "anc5",
    TRUE ~ population
  ))

df.s <- df.s %>% 
  rename(level = salt,
         days = time,
         well.ID = well) %>% 
  select(well.ID, RFU, population, days, level) %>%
  arrange(factor(population, levels = pop_order))

df.s$rep <- 1
df.s$grad <- "Salt"
df.s$unit <- "g/L NaCl"

df.t <- read.csv('data/01_temp_rfus_time.csv') # Raw growth data, temperature
head(df.t)

df.t <- df.t %>% 
  filter(population != "cc1629",
         temperature != 20) %>% 
  rename(well.ID = well_plate,
         level = temperature) %>% 
  select(well.ID, RFU, population, days, level) %>% 
  arrange(factor(population, levels = pop_order), well.ID)

rep <- numeric()

for (i in unique(df.t$population)){
  
  n <- nrow(df.t[df.t$population == i,])/4
  reps <- c(rep(1,n), rep(2,n), rep(3,n), rep(4,n))
  rep<- c (rep, reps)
  
}

df.t <- cbind(df.t, rep)
df.t$grad <- "Temperature"
df.t$unit <- "°C"

df.raw <- rbind (df.l, df.n, df.p, df.s, df.t)

df.raw$logRFU <- log(df.raw$RFU + 0.001)
df.raw$grad <- as.factor(df.raw$grad)

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Pick a population and replicate"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop", "Select a population:",
        choices  = pop_order,        
        selected = pop_order[1],
        width    = "100%"
      ),
      actionButton("prev_pop", "Previous population"),
      actionButton("next_pop", "Next population")
    )
  ),
  
  fluidRow(
    column(
      3,
      selectInput(
        "rep", "Select a replicate:",
        choices  = rep_levels,
        selected = rep_levels[1],
        width    = "100%"
      ),
      actionButton("prev_rep", "Previous replicate"),
      actionButton("next_rep", "Next replicate")
    )
  ),
  
  fluidRow(
    column(4, textOutput("anc.text")),
    column(4, ""),
    column(4, textOutput("evol.text"))
  ),
  
  hr(),
  
  fluidRow(
    column(4, plotOutput("p.light")),
    column(4, plotOutput("p.nit")),
    column(4, plotOutput("p.phos"))
  ),
  
  fluidRow(
    column(4, plotOutput("p.salt")),
    column(4, plotOutput("p.temp")),
    column(4, "")
  ),
  
  hr(),
  hr(),
  
  fluidRow(
    
    column(3,
           selectInput("grad", "Select a gradient:",
                       choices = sort(unique(df.raw$grad)),
                       width = '100%')
  ),
  
  fluidRow(column(12, "")),
  
  fluidRow(
    column(3, plotOutput("p.1")),
    column(3, plotOutput("p.2")),
    column(3, plotOutput("p.3")),
    column(3, plotOutput("p.4"))
  ),
  
  fluidRow(
    column(3, plotOutput("p.5")),
    column(3, plotOutput("p.6")),
    column(3, plotOutput("p.7")),
    column(3, plotOutput("p.8"))
  ),
  
  fluidRow(
    column(3, plotOutput("p.9")),
    column(3, plotOutput("p.10")),
    column(3, ""),
    column(3, "")
  )
  
  )
  
)


server <- function(input, output, session) {
  
  ###### Load the information ######
  
  observeEvent(input$next_pop, { # Next population
    req(input$pop)
    i <- match(input$pop, pop_order)
    next_i <- ifelse(i == length(pop_order), 1, i + 1)
    updateSelectInput(session, "pop", selected = pop_order[next_i])
  })
  
  observeEvent(input$prev_pop, {   # Previous population
    req(input$pop)
    i <- match(input$pop, pop_order)
    prev_i <- ifelse(i == 1, length(pop_order), i - 1)
    updateSelectInput(session, "pop", selected = pop_order[prev_i])
  })
  
  observeEvent(input$next_rep, {   # Next replicate
    req(input$rep)
    j <- match(as.numeric(input$rep), rep_levels)
    next_j <- ifelse(j == length(rep_levels), 1, j + 1)
    updateSelectInput(session, "rep", selected = rep_levels[next_j])
  })
  
  observeEvent(input$prev_rep, {   # Previous replicate
    req(input$rep)
    j <- match(as.numeric(input$rep), rep_levels)
    prev_j <- ifelse(j == 1, length(rep_levels), j - 1)
    updateSelectInput(session, "rep", selected = rep_levels[prev_j])
  })
  
  sel.row <- reactive({                               # Which row was selected based on pop
    df.sum %>% filter(population == input$pop)
  })
  
  output$anc.text <- renderText({                     # Left column: ancestor info
    req(sel.row())
    paste("Descended from population:", sel.row()$Anc[1])
  })
  
  output$evol.text <- renderText({                    # Right column: evolutionary condition
    req(sel.row())
    paste("Evolutionary condition:", sel.row()$Evol[1])
  })
  
  ###### Generate curve plots ######
  
  output$p.light <- renderPlot({
    
    l.mu <- df.l.mu %>% filter(population == input$pop, rep == input$rep)
    l.monod <- df.sum %>% filter(population == input$pop, rep == input$rep)
    
    curve.l <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = l.monod$I.µ.max, k.s = l.monod$I.K.s)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = l.mu,
                 aes(x = light, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (µM photons m^-2s^-1)",
        y = "Exponential growth rate",
        title = "Monod curve — Light"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = l.monod$I.µ.max, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 1/l.monod$I.comp, linetype = "dashed", colour = "darkred")
  })
  
  output$p.nit <- renderPlot({
    
    n.mu <- df.n.mu %>% filter(population == input$pop, rep == input$rep)
    n.monod <- df.sum %>% filter(population == input$pop, rep == input$rep)
    
    curve.n <- tibble::tibble(
      res  = seq(0, 1000, length.out = 200),
      rate = pred_mon(res, r.max = n.monod$N.µ.max, k.s = n.monod$N.K.s)
    )
    
    ggplot(curve.n, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = n.mu,
                 aes(x = nit, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Monod curve — Nitrogen"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = n.monod$N.µ.max, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 1/n.monod$N.comp, linetype = "dashed", colour = "darkred")
  })
  
  output$p.phos <- renderPlot({
    
    p.mu <- df.p.mu %>% filter(population == input$pop, rep == input$rep)
    p.monod <- df.sum %>% filter(population == input$pop, rep == input$rep)
    
    curve.p <- tibble::tibble(
      res  = seq(0, 50, length.out = 200),
      rate = pred_mon(res, r.max = p.monod$P.µ.max, k.s = p.monod$P.K.s)
    )
    
    ggplot(curve.p, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "darkorange2") +
      geom_point(data = p.mu,
                 aes(x = phos, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Phosphate concentration (µM)",
        y = "Exponential growth rate",
        title = "Monod curve — Phosphorous"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = p.monod$P.µ.max, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 1/p.monod$P.comp, linetype = "dashed", colour = "darkred")
    })
  
  output$p.salt <- renderPlot({
    
    s.mu <- df.s.mu %>% filter(population == input$pop)
    s.tol <- df.sum %>% filter(population == input$pop)
    
    curve.s <- tibble::tibble(
      salt  = seq(0, 10, length.out = 200),
      rate = pred_salt(salt, a = s.tol$S.µ.max, b = s.tol$S.b, c = s.tol$S.c)
    )
    
    ggplot(curve.s, aes(x = salt, y = rate)) +
      geom_line(size = 1.5, colour = "orchid2") +
      geom_point(data = s.mu,
                 aes(x = salt, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Salt concentration (g/L)",
        y = "Exponential growth rate",
        title = "Salt tolerance curve"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = s.tol$S.µ.max, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = s.tol$S.c, linetype = "dashed", colour = "darkred")
  })
  
  output$p.temp <- renderPlot({
    
    t.mu <- df.t.mu %>% filter(population == input$pop, rep == input$rep)
    t.tpc <- df.sum %>% filter(population == input$pop, rep == input$rep)
    
    curve.t <- tibble::tibble(
     res  = seq(0, 42, length.out = 200),
     rate = pred_lact(res, a = t.tpc$T.a, b = t.tpc$T.b, delta_t = t.tpc$T.d.t, tmax = t.tpc$T.tmax)
    )
    
    ggplot(curve.t, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = t.mu,
                 aes(x = temp, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Thermal performance curve (Lactin II)"
      ) +
      theme_classic() +
      ylim(-0.1, 6) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_vline(xintercept = max(-1.8, t.tpc$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = t.tpc$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = t.tpc$T.opt, linetype = "dashed", colour = 'forestgreen') +
      
      geom_segment(
        aes(x = t.tpc$T.br.min, xend = t.tpc$T.br.max,
            y = t.tpc$T.µ.max / 2, yend = t.tpc$T.µ.max / 2)
      )
    
  })
  
  ###### Generate µ estimation plots ######
  
  output$p.1 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p1 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p1 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p1$level)))
    df.p1 <- df.p1 %>% filter(level == levels[1])
      
    df.p1 <- df.p1[order(df.p1$days), ]
      
    t.series <- unique(df.p1$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p1.sl <- df.p1[df.p1$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p1.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p1.th <- df.p1[df.p1$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
        
    df.p1.th$N0 <- df.p1$RFU[1]
        
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p1.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
        
    smt.days <- seq(min(df.p1.th$days, na.rm = TRUE), max(df.p1.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
        
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
        
    ggplot() + 
      geom_point(data=df.p1, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p1$level, df.p1$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.2 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p2 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p2 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p2$level)))
    df.p2 <- df.p2 %>% filter(level == levels[2])
    
    df.p2 <- df.p2[order(df.p2$days), ]
    
    t.series <- unique(df.p2$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p2.sl <- df.p2[df.p2$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p2.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p2.th <- df.p2[df.p2$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p2.th$N0 <- df.p2$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p2.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p2.th$days, na.rm = TRUE), max(df.p2.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p2, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p2$level, df.p2$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.3 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p3 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p3 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p3$level)))
    df.p3 <- df.p3 %>% filter(level == levels[3])
    
    df.p3 <- df.p3[order(df.p3$days), ]
    
    t.series <- unique(df.p3$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p3.sl <- df.p3[df.p3$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p3.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p3.th <- df.p3[df.p3$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p3.th$N0 <- df.p3$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p3.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p3.th$days, na.rm = TRUE), max(df.p3.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p3, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p3$level, df.p3$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.4 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p4 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p4 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p4$level)))
    df.p4 <- df.p4 %>% filter(level == levels[4])
    
    df.p4 <- df.p4[order(df.p4$days), ]
    
    t.series <- unique(df.p4$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p4.sl <- df.p4[df.p4$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p4.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p4.th <- df.p4[df.p4$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p4.th$N0 <- df.p4$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p4.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p4.th$days, na.rm = TRUE), max(df.p4.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p4, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p4$level, df.p4$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.5 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p5 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p5 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p5$level)))
    df.p5 <- df.p5 %>% filter(level == levels[5])
    
    df.p5 <- df.p5[order(df.p5$days), ]
    
    t.series <- unique(df.p5$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p5.sl <- df.p5[df.p5$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p5.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p5.th <- df.p5[df.p5$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p5.th$N0 <- df.p5$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p5.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p5.th$days, na.rm = TRUE), max(df.p5.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p5, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p5$level, df.p5$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.6 <- renderPlot({
    
    if(input$grad %in% c("Light", "Nitrogen", "Phosphorous", "Temperature")){
      
      df.p6 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p6 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    }
    
    levels <- sort(unique(as.numeric(df.p6$level)))
    df.p6 <- df.p6 %>% filter(level == levels[6])
    
    df.p6 <- df.p6[order(df.p6$days), ]
    
    t.series <- unique(df.p6$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p6.sl <- df.p6[df.p6$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p6.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p6.th <- df.p6[df.p6$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p6.th$N0 <- df.p6$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p6.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p6.th$days, na.rm = TRUE), max(df.p6.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p6, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p6$level, df.p6$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.7 <- renderPlot({
    
    if(input$grad == "Temperature"){
      
      empty_na_plot()
      
    } else if(input$grad %in% c("Light", "Nitrogen", "Phosphorous")){
      
      df.p7 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p7 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    } 
    
    levels <- sort(unique(as.numeric(df.p7$level)))
    df.p7 <- df.p7 %>% filter(level == levels[7])
    
    df.p7 <- df.p7[order(df.p7$days), ]
    
    t.series <- unique(df.p7$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p7.sl <- df.p7[df.p7$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p7.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p7.th <- df.p7[df.p7$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p7.th$N0 <- df.p7$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p7.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p7.th$days, na.rm = TRUE), max(df.p7.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p7, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p7$level, df.p7$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.8 <- renderPlot({
    
    if(input$grad == "Temperature"){
      
      empty_na_plot()
      
    } else if(input$grad %in% c("Light", "Nitrogen", "Phosphorous")){
      
      df.p8 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p8 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    } 
    
    levels <- sort(unique(as.numeric(df.p8$level)))
    df.p8 <- df.p8 %>% filter(level == levels[8])
    
    df.p8 <- df.p8[order(df.p8$days), ]
    
    t.series <- unique(df.p8$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p8.sl <- df.p8[df.p8$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p8.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p8.th <- df.p8[df.p8$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p8.th$N0 <- df.p8$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p8.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p8.th$days, na.rm = TRUE), max(df.p8.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p8, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p8$level, df.p8$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.9 <- renderPlot({
    
    if(input$grad == "Temperature"){
      
      empty_na_plot()
      
    } else if(input$grad %in% c("Light", "Nitrogen", "Phosphorous")){
      
      df.p9 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p9 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    } 
    
    levels <- sort(unique(as.numeric(df.p9$level)))
    df.p9 <- df.p9 %>% filter(level == levels[9])
    
    df.p9 <- df.p9[order(df.p9$days), ]
    
    t.series <- unique(df.p9$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p9.sl <- df.p9[df.p9$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p9.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p9.th <- df.p9[df.p9$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p9.th$N0 <- df.p9$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p9.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p9.th$days, na.rm = TRUE), max(df.p9.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p9, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p9$level, df.p9$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.10 <- renderPlot({
    
    if(input$grad == "Temperature"){
      
      empty_na_plot()
      
    } else if(input$grad %in% c("Light", "Nitrogen", "Phosphorous")){
      
      df.p10 <- df.raw %>% filter(population == input$pop, rep == input$rep, grad == input$grad)
      
    } else if(input$grad == "Salt"){
      
      df.p10 <- df.raw %>% filter(population == input$pop, grad == input$grad)
      
    } 
    
    levels <- sort(unique(as.numeric(df.p10$level)))
    df.p10 <- df.p10 %>% filter(level == levels[10])
    
    df.p10 <- df.p10[order(df.p10$days), ]
    
    t.series <- unique(df.p10$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including days 0 
      
      df.p10.sl <- df.p10[df.p10$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(logRFU~days, data = df.p10.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p10.th <- df.p10[df.p10$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
    
    df.p10.th$N0 <- df.p10$RFU[1]
    
    r_exp <- nls_multstart(RFU ~ N0 * exp(r*days), # Exponential growth model (N0 is in our dataframe)
                           data = df.p10.th,
                           start_lower = c(r = -4.5), 
                           start_upper = c(r = 4.5),   
                           iter = 500,
                           supp_errors = 'Y',
                           control = nls.control(maxiter = 200))
    
    smt.days <- seq(min(df.p10.th$days, na.rm = TRUE), max(df.p10.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p10, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Level", df.p10$level, df.p10$unit), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
}

shinyApp(ui, server)

