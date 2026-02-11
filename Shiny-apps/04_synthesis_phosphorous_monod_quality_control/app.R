# J R Laurich

# January 3rd, 2026

# We're going to take a look at my interspecific phosphorous Monod models and data to check fits.

# That's the Bestion 2018, Narwani 2015, Edwards 2015, and Levasseur 2025 datasets
# Edwards 2015 and Narwani 2015 only report estimates, not raw data. 

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)

# Functions -------------------------------------------------

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

# Load & explore the data --------------------------------------------------------

###### Bestion 2018 ######

df.b.raw <- read.csv('data/61_Bestion_2018_phos_monods.csv') # Bestion raw data
head(df.b.raw)

df.b.raw <-df.b.raw %>% 
  rename(temp = Temperature_c,
         phos = Phosphate_c,
         Species.name = SpeciesName) %>% 
  mutate(id.number = as.integer(factor(df.b.raw$SpeciesNb)))

df.b.fit <- read.csv('data/62_Bestion_2018_phos_monods_fits.csv') # Bestion phos monod fits
head(df.b.fit)

pop_order.b <- unique(df.b.fit$Sp.id)

###### Levasseur 2025 ######

df.lv.raw <- read.csv('data/63_Levasseur_2025_µ_estimates_phosphorous.csv') # Levasseur growth estimates
head(df.lv.raw)

df.lv.raw <-df.lv.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.lv.raw$Sp.id)))

df.lv.fit <- read.csv('data/64_Levasseur_2025_phos_monods.csvv') # Levasseur phos monod fits
head(df.lv.fit)

pop_order.lv <- unique(df.lv.fit$Sp.id)

# Build the app -----------------------------------------------------------

ui <- fluidPage(
  
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
    column(5, plotOutput("p.b.phos1")),
    column(2, ""),
    column(5, plotOutput("p.b.phos2"))
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
    column(5, plotOutput("p.lv.phos1")),
    column(2, ""),
    column(5, plotOutput("p.lv.phos2"))
  )
  
)

server <- function(input, output, session) {
  
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
    paste("Species:", sel.row.b()$Species.name[1])
  })
  
  output$p.b.phos1 <- renderPlot({
    
    mu.b <- df.b.raw %>% filter(id.number == input$pop.b) 
    phos.b <- df.b.fit %>% filter(Sp.id == input$pop.b)
    
    t.sel <- ifelse(is.na(phos.b$T.opt.TPC), phos.b$T.max.raw, phos.b$T.opt.TPC)
    
    t.sel2 <- mu.b$temp[which.min(abs(mu.b$temp - t.sel))]
    
    mu.b <- mu.b %>% 
      filter(temp == t.sel2)
    
    curve.b <- tibble::tibble(
      res  = seq(0, 50, length.out = 200),
      rate = pred_mon(res, r.max = phos.b$r.max.mod, k.s = phos.b$K.s.mod)
    )
    
    ggplot(curve.b, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "darkorange2") +
      geom_point(data = mu.b,
                 aes(x = phos, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Phosphate concentration (µM)",
        y = "Exponential growth rate",
        title = "Bestion 2018: Phosphorous Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = phos.b$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*phos.b$K.s.mod/(phos.b$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.b.phos2 <- renderPlot({
    
    mu.b <- df.b.raw %>% filter(id.number == input$pop.b) 
    phos.b <- df.b.fit %>% filter(Sp.id == input$pop.b)
    
    t.sel <- ifelse(is.na(phos.b$T.opt.TPC), phos.b$T.max.raw, phos.b$T.opt.TPC)
    
    t.sel2 <- mu.b$temp[which.min(abs(mu.b$temp - t.sel))]
    
    mu.b <- mu.b %>% 
      filter(temp == t.sel2)
    
    curve.b <- tibble::tibble(
      res  = seq(0, 50, length.out = 200),
      rate = pred_mon(res, r.max = phos.b$r.max.post, k.s = phos.b$K.s.post)
    )
    
    ggplot(curve.b, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "darkorange2") +
      geom_point(data = mu.b,
                 aes(x = phos, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Phosphate concentration (µM)",
        y = "Exponential growth rate",
        title = "Bestion 2018: Phosphorous Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = phos.b$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*phos.b$K.s.post/(phos.b$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
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
    df.lv.raw %>% filter(id.number == as.numeric(input$pop.lv))
  })
  
  output$sp.text.lv <- renderText({                     # Left column: species
    req(sel.row.lv())
    paste("Species:", sel.row.lv()$Species.name[1])
  })
  
  output$p.lv.phos1 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    phos.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(phos.lv$T.opt.TPC), phos.lv$T.max.raw, phos.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 50, length.out = 200),
      rate = pred_mon(res, r.max = phos.lv$r.max.mod, k.s = phos.lv$K.s.mod)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "darkorange2") +
      geom_point(data = mu.lv,
                 aes(x = phos, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Phosphate concentration (µM)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Phosphorous Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = phos.lv$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*phos.lv$K.s.mod/(phos.lv$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.lv.phos2 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    phos.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(phos.lv$T.opt.TPC), phos.lv$T.max.raw, phos.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 50, length.out = 200),
      rate = pred_mon(res, r.max = phos.lv$r.max.post, k.s = phos.lv$K.s.post)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "darkorange2") +
      geom_point(data = mu.lv,
                 aes(x = phos, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Phosphate concentration (µM)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Phosphorous Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = phos.lv$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*phos.lv$K.s.post/(phos.lv$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
  })
  
}

shinyApp(ui, server)
