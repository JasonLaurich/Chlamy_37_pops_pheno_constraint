# J R Laurich

# December 15, 2025

# We're going to take a look at my interspecific TPC models and data to check fits.

# That's the Thomas 2012, Bestion 2018, Edwards 2016, and Lewington-Pearce 2019 datasets

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(nls.multstart)

# Functions -------------------------------------------------

pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

# Load & explore the data --------------------------------------------------------

###### Thomas 2012 ######

df.t.raw <- read.csv('data/28_Thomas_2012_raw_data.csv') # Thomas raw data
head(df.t.raw)

df.t.raw <-df.t.raw %>% 
  rename(temp = Temperature,
         mu = Growth.rate)

df.t.fit <- read.csv('data/29_Thomas_2012_TPCs.csv') # Thomas TPCs
head(df.t.fit)

df.t.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter((T.max.max-T.max.min)>= 0.25*(T.max-T.min)) %>% 
  nrow() # We lose one population at 25% (7 at 15%), leaving us with 168 or 162/169 fitted curves

df.t.fit <- df.t.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter(!(T.max.max-T.max.min)>= 0.25*(T.max-T.min)) 

pop_order.t <- unique(df.t.fit$Sp.id) # Unique numbers for each species. 

###### Bestion 2018 ######

df.b.raw <- read.csv('data/31_Bestion_2018_raw_data.csv') # Bestion raw data
head(df.b.raw)

df.b.raw <- df.b.raw[df.b.raw$Phosphate_c == 30,]

df.b.fit <- read.csv('data/32_Bestion_2018_TPCs.csv') # Bestion TPCs
head(df.b.fit)

df.b.fit <- df.b.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter(!(T.max.max-T.max.min)>= 0.25*(T.max-T.min)) # All 5 species pass this test. 

pop_order.b <- unique(df.b.fit$Sp.id) # Unique numbers for each species.

###### Lewington-Pearce 2019 ######

df.l.raw <- read.csv('data/35_Lewington-Pearce_2019_µ_estimates_temp.csv') # Lewington-Pearce raw data
head(df.l.raw)

df.l.raw$Sp.num <- c(rep(1,12), rep(2,12), rep(3,12), rep(4,12), rep(5,12), rep(6,12), 7)

df.l.fit <- read.csv('data/36_Lewington_2019_TPCs.csv') # Lewington TPC fits
head(df.l.fit)

df.l.fit <- df.l.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter(!(T.max.max-T.max.min)>= 0.25*(T.max-T.min)) # All 6 species pass this test. 

pop_order.l <- unique(df.l.fit$Sp.id) # Unique numbers for each species.

###### Edwards 2016 ######

df.e.raw <- read.csv('data/38_Edwards_2016_raw_data.csv') # Edwards raw data
head(df.e.raw)

df.e.raw <- df.e.raw %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = paste(species, reference, sep = "_")) %>% # So that each entry is treated seperately!
  filter(irradiance <= 250) # Eliminate the insanely high light data. 

df.e.raw$sp.num <- as.integer(factor(df.e.raw$unique.id))

df.e.fit <- read.csv('data/39_Edwards_2016_TPCs.csv') # Edwards TPC fits
head(df.e.fit)

df.e.fit <- df.e.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter(!(T.max.max-T.max.min)>= 0.25*(T.max-T.min)) # All 29 species pass this test. 

pop_order.e <- unique(df.e.fit$Sp.id) # Unique numbers for each species.

###### Levasseur 2025 ######

df.lv.raw <- read.csv('data/44_Levasseur_2025_µ_estimates_temp.csv') # Levasseur raw data
head(df.lv.raw)

df.lv.fit <- read.csv('data/45_Levasseur_2025_TPCs.csv') # Levasseur TPC fits
head(df.lv.fit)

df.lv.fit <- df.lv.fit %>% 
  mutate(T.min = if_else(is.na(T.min), -1.8, T.min)) %>% 
  filter(!(T.max.max-T.max.min)>= 0.25*(T.max-T.min)) # All 20 species pass this test. 

pop_order.lv <- unique(df.lv.fit$Sp.id) # Unique numbers for each species.

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  
  ###### Thomas 2012 ######
  
  fluidRow(
    column(12, align = "center",
           tags$h3(tags$b("Thomas 2012 dataset"))
    )
  ),
  
  titlePanel("Pick a species number"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop.t", "Select a number:",
        choices  = pop_order.t,        
        selected = pop_order.t[1],
        width    = "100%"
      ),
      actionButton("prev_pop.t", "Previous species number"),
      actionButton("next_pop.t", "Next species number")
    )
  ),
  
  fluidRow(
    column(4, textOutput("sp.text.t")),
    column(4, ""),
    column(4, textOutput("study.text.t"))
  ),
  
  fluidRow(
    column(5, plotOutput("p.t.tpc1")),
    column(2, ""),
    column(5, plotOutput("p.t.tpc2"))
  ),
  
  ###### Bestion 2018 ######
  
  fluidRow(
    column(12, align = "center",
           tags$h3(tags$b("Bestion 2018 dataset"))
    )
  ),
  
  titlePanel("Pick a species number"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop.b", "Select a number:",
        choices  = pop_order.b,        
        selected = pop_order.b[1],
        width    = "100%"
      ),
      actionButton("prev_pop.b", "Previous species number"),
      actionButton("next_pop.b", "Next species number")
    )
  ),
  
  fluidRow(
    column(4, textOutput("sp.text.b")),
    column(4, ""),
    column(4, "")
  ),
  
  fluidRow(
    column(5, plotOutput("p.b.tpc1")),
    column(2, ""),
    column(5, plotOutput("p.b.tpc2"))
  ),
  
  ###### Lewington-Pearce 2019 ######
  
  fluidRow(
    column(12, align = "center",
           tags$h3(tags$b("Lewington-Pearce 2019 dataset"))
    )
  ),
  
  titlePanel("Pick a species number"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop.l", "Select a number:",
        choices  = pop_order.l,        
        selected = pop_order.l[1],
        width    = "100%"
      ),
      actionButton("prev_pop.l", "Previous species number"),
      actionButton("next_pop.l", "Next species number")
    )
  ),
  
  fluidRow(
    column(4, textOutput("sp.text.l")),
    column(4, ""),
    column(4, "")
  ),
  
  fluidRow(
    column(5, plotOutput("p.l.tpc1")),
    column(2, ""),
    column(5, plotOutput("p.l.tpc2"))
  ),
  
  ###### Edwards 2016 ######
  
  fluidRow(
    column(12, align = "center",
           tags$h3(tags$b("Edwards 2016 dataset"))
    )
  ),
  
  titlePanel("Pick a species number"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop.e", "Select a number:",
        choices  = pop_order.e,        
        selected = pop_order.e[1],
        width    = "100%"
      ),
      actionButton("prev_pop.e", "Previous species number"),
      actionButton("next_pop.e", "Next species number")
    )
  ),
  
  fluidRow(
    column(4, textOutput("sp.text.e")),
    column(4, ""),
    column(4, textOutput("study.text.e"))
  ),
  
  fluidRow(
    column(5, plotOutput("p.e.tpc1")),
    column(2, ""),
    column(5, plotOutput("p.e.tpc2"))
  ),
  
  ###### Levasseur 2025 ######
  
  fluidRow(
    column(12, align = "center",
           tags$h3(tags$b("Levasseur 2025 dataset"))
    )
  ),
  
  titlePanel("Pick a species number"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "pop.lv", "Select a number:",
        choices  = pop_order.lv,        
        selected = pop_order.lv[1],
        width    = "100%"
      ),
      actionButton("prev_pop.lv", "Previous species number"),
      actionButton("next_pop.lv", "Next species number")
    )
  ),
  
  fluidRow(
    column(4, textOutput("sp.text.lv")),
    column(4, ""),
    column(4, "")
  ),
  
  fluidRow(
    column(5, plotOutput("p.lv.tpc1")),
    column(2, ""),
    column(5, plotOutput("p.lv.tpc2"))
  )
  
)

server <- function(input, output, session) {
  
  ###### Thomas 2012 ######
  
  observeEvent(input$next_pop.t, { # Next population
    req(input$pop.t)
    i <- match(input$pop.t, pop_order.t)
    next_i.t <- ifelse(i == length(pop_order.t), 1, i + 1)
    updateSelectInput(session, "pop.t", selected = pop_order.t[next_i.t])
  })
  
  observeEvent(input$prev_pop.t, {   # Previous population
    req(input$pop.t)
    i <- match(input$pop.t, pop_order.t)
    prev_i.t <- ifelse(i == 1, length(pop_order.t), i - 1)
    updateSelectInput(session, "pop.t", selected = pop_order.t[prev_i.t])
  })
  
  sel.row.t <- reactive({
    req(input$pop.t)
    df.t.raw %>% filter(id.number == as.numeric(input$pop.t))
  })
  
  output$sp.text.t <- renderText({                     # Left column: species
    req(sel.row.t())
    paste("Species:", sel.row.t()$Species.name[1])
  })
  
  output$study.text.t <- renderText({                    # Right column: study
    req(sel.row.t())
    paste("Study:", sel.row.t()$Study[1])
  })
  
  output$p.t.tpc1 <- renderPlot({
    
    mu.t <- df.t.raw %>% filter(id.number == input$pop.t)
    tpc.t <- df.t.fit %>% filter(Sp.id == input$pop.t)
    
    curve.t <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.t$a.mod, b = tpc.t$b.mod, delta_t = tpc.t$d.t.mod, tmax = tpc.t$tmax.mod)
    )
    
    ggplot(curve.t, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.t,
                 aes(x = temp, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Thomas 2012: Lactin II TPC (nls.LM model fits)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.t$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.t$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.t$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.t.tpc2 <- renderPlot({
    
    mu.t <- df.t.raw %>% filter(id.number == input$pop.t)
    tpc.t <- df.t.fit %>% filter(Sp.id == input$pop.t)
    
    curve.t <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.t$a, b = tpc.t$b, delta_t = tpc.t$d.t, tmax = tpc.t$tmax)
    )
    
    ggplot(curve.t, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.t,
                 aes(x = temp, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Thomas 2012: Lactin II TPC (bootstrapped medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.t$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.t$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.t$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  ###### Bestion 2018 ######
  
  observeEvent(input$next_pop.b, { # Next population
    req(input$pop.b)
    i <- match(input$pop.b, pop_order.b)
    next_i.b <- ifelse(i == length(pop_order.b), 1, i + 1)
    updateSelectInput(session, "pop.b", selected = pop_order.b[next_i.b])
  })
  
  observeEvent(input$prev_pop.b, {   # Previous population
    req(input$pop.b)
    i <- match(input$pop.b, pop_order.b)
    prev_i.b <- ifelse(i == 1, length(pop_order.b), i - 1)
    updateSelectInput(session, "pop.b", selected = pop_order.b[prev_i.b])
  })
  
  sel.row.b <- reactive({
    req(input$pop.b)
    df.b.raw %>% filter(SpeciesNb == as.numeric(input$pop.b))
  })
  
  output$sp.text.b <- renderText({                     # Left column: species
    req(sel.row.b())
    paste("Species:", sel.row.b()$SpeciesName[1])
  })
  
  output$p.b.tpc1 <- renderPlot({
    
    mu.b <- df.b.raw %>% filter(SpeciesNb == input$pop.b)
    tpc.b <- df.b.fit %>% filter(Sp.id == input$pop.b)
    
    curve.b <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.b$a.mod, b = tpc.b$b.mod, delta_t = tpc.b$d.t.mod, tmax = tpc.b$tmax.mod)
    )
    
    ggplot(curve.b, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.b,
                 aes(x = Temperature_c, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Bestion 2018: Lactin II TPC (nls.LM model fits)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.b$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.b$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.b$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.b.tpc2 <- renderPlot({
    
    mu.b <- df.b.raw %>% filter(SpeciesNb == input$pop.b)
    tpc.b <- df.b.fit %>% filter(Sp.id == input$pop.b)
    
    curve.b <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.b$a, b = tpc.b$b, delta_t = tpc.b$d.t, tmax = tpc.b$tmax)
    )
    
    ggplot(curve.b, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.b,
                 aes(x = Temperature_c, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Bestion 2018: Lactin II TPC (bootstrapped medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.b$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.b$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.b$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  ###### Lewington-Pearce 2019 ######
  
  observeEvent(input$next_pop.l, { # Next population
    req(input$pop.l)
    i <- match(input$pop.l, pop_order.l)
    next_i.l <- ifelse(i == length(pop_order.l), 1, i + 1)
    updateSelectInput(session, "pop.l", selected = pop_order.l[next_i.l])
  })
  
  observeEvent(input$prev_pop.l, {   # Previous population
    req(input$pop.l)
    i <- match(input$pop.l, pop_order.l)
    prev_i.l <- ifelse(i == 1, length(pop_order.l), i - 1)
    updateSelectInput(session, "pop.l", selected = pop_order.l[prev_i.l])
  })
  
  sel.row.l <- reactive({
    req(input$pop.l)
    df.l.raw %>% filter(Sp.num == as.numeric(input$pop.l))
  })
  
  output$sp.text.l <- renderText({                     # Left column: species
    req(sel.row.l())
    paste("Species:", sel.row.l()$Sp.id[1])
  })
  
  output$p.l.tpc1 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(Sp.num == input$pop.l)
    tpc.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    curve.l <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.l$a.mod, b = tpc.l$b.mod, delta_t = tpc.l$d.t.mod, tmax = tpc.l$tmax.mod)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.l,
                 aes(x = temperature, y = r.exp),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Lactin II TPC (nls.LM model fits)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.l$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.l$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.l$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.l.tpc2 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(Sp.num == input$pop.l)
    tpc.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    curve.l <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.l$a, b = tpc.l$b, delta_t = tpc.l$d.t, tmax = tpc.l$tmax)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.l,
                 aes(x = temperature, y = r.exp),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Lactin II TPC (bootstrapped medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.l$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.l$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.l$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  ###### Edwards 2016 ######
  
  observeEvent(input$next_pop.e, { # Next population
    req(input$pop.e)
    i <- match(input$pop.e, pop_order.e)
    next_i.e <- ifelse(i == length(pop_order.e), 1, i + 1)
    updateSelectInput(session, "pop.e", selected = pop_order.e[next_i.e])
  })
  
  observeEvent(input$prev_pop.e, {   # Previous population
    req(input$pop.e)
    i <- match(input$pop.e, pop_order.e)
    prev_i.e <- ifelse(i == 1, length(pop_order.e), i - 1)
    updateSelectInput(session, "pop.e", selected = pop_order.e[prev_i.e])
  })
  
  sel.row.e <- reactive({
    req(input$pop.e)
    df.e.raw %>% filter(sp.num == as.numeric(input$pop.e))
  })
  
  output$sp.text.e <- renderText({                     # Left column: species
    req(sel.row.e())
    paste("Species:", sel.row.e()$species[1])
  })
  
  output$study.text.e <- renderText({                    # Right column: study
    req(sel.row.e())
    paste("Study:", sel.row.e()$reference[1])
  })
  
  output$p.e.tpc1 <- renderPlot({
    
    mu.e <- df.e.raw %>% filter(sp.num == input$pop.e) %>% filter(irradiance == max(irradiance))
    tpc.e <- df.e.fit %>% filter(Sp.id == input$pop.e)
    
    curve.e <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.e$a.mod, b = tpc.e$b.mod, delta_t = tpc.e$d.t.mod, tmax = tpc.e$tmax.mod)
    )
    
    ggplot(curve.e, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.e,
                 aes(x = temperature, y = growth.rate),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Edwards 2016: Lactin II TPC (nls.LM model fits)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.e$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.e$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.e$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.e.tpc2 <- renderPlot({
    
    mu.e <- df.e.raw %>% filter(sp.num == input$pop.e) %>% filter(irradiance == max(irradiance))
    tpc.e <- df.e.fit %>% filter(Sp.id == input$pop.e)
    
    curve.e <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.e$a, b = tpc.e$b, delta_t = tpc.e$d.t, tmax = tpc.e$tmax)
    )
    
    ggplot(curve.e, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.e,
                 aes(x = temperature, y = growth.rate),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Edwards 2016: Lactin II TPC (bootstrapped medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.e$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.e$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.e$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  ###### Levasseur 2025 ######
  
  observeEvent(input$next_pop.lv, { # Next population
    req(input$pop.lv)
    i <- match(input$pop.lv, pop_order.lv)
    next_i.lv <- ifelse(i == length(pop_order.lv), 1, i + 1)
    updateSelectInput(session, "pop.lv", selected = pop_order.lv[next_i.lv])
  })
  
  observeEvent(input$prev_pop.lv, {   # Previous population
    req(input$pop.lv)
    i <- match(input$pop.lv, pop_order.lv)
    prev_i.lv <- ifelse(i == 1, length(pop_order.lv), i - 1)
    updateSelectInput(session, "pop.lv", selected = pop_order.lv[prev_i.lv])
  })
  
  sel.row.lv <- reactive({
    req(input$pop.lv)
    df.lv.raw %>% filter(id == as.numeric(input$pop.lv))
  })
  
  output$sp.text.lv <- renderText({                     # Left column: species
    req(sel.row.lv())
    paste("Species:", sel.row.lv()$Sp.id[1])
  })
  
  output$p.lv.tpc1 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id == input$pop.lv)
    tpc.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    curve.lv <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.lv$a.mod, b = tpc.lv$b.mod, delta_t = tpc.lv$d.t.mod, tmax = tpc.lv$tmax.mod)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.lv,
                 aes(x = temperature, y = r.exp),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Lactin II TPC (nls.LM model fits)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.lv$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.lv$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.lv$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.lv.tpc2 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id == input$pop.lv)
    tpc.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    curve.lv <- tibble::tibble(
      res  = seq(-10, 50, length.out = 200),
      rate = pred_lact(res, a = tpc.lv$a, b = tpc.lv$b, delta_t = tpc.lv$d.t, tmax = tpc.lv$tmax)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.lv,
                 aes(x = temperature, y = r.exp),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Lactin II TPC (bootstrapped medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 4) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = max(-1.8, tpc.lv$T.min), linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.lv$T.max, linetype = "dashed", colour = 'red3') +
      geom_vline(xintercept = tpc.lv$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
}

shinyApp(ui, server)

