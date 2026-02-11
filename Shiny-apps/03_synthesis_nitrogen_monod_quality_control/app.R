# J R Laurich

# January 3rd, 2026

# We're going to take a look at my interspecific nitrogen Monod models and data to check fits.

# That's the Edwards 2015, Lewington-Pearce 2019, and Levasseur 2025 datasets
# Edwards 2015 only reports estimates, not raw data. 

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)

# Functions -------------------------------------------------

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

# Load & explore the data --------------------------------------------------------

###### Lewington-Pearce 2019 ######

df.l.raw <- read.csv('data/55_Lewington-Pearce_2019_µ_estimates_nitrogen.csv') # Lewington-Pearce raw data
head(df.l.raw)

df.l.raw <-df.l.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.l.raw$Sp.id)))

df.l.fit <- read.csv('data/56_Lewington_2019_nit_monods.csv') # Lewington nit monod fits
head(df.l.fit)

pop_order.l <- unique(df.l.fit$Sp.id)

###### Levasseur 2025 ######

df.lv.raw <- read.csv('data/58_Levasseur_2025_µ_estimates_nitrogen.csv') # Levasseur raw data
head(df.lv.raw)

df.lv.raw <-df.lv.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.lv.raw$Sp.id)))

df.lv.fit <- read.csv('data/59_Levasseur_2025_nit_monods.csv') # Levasseur nit monod fits
head(df.lv.fit)

pop_order.lv <- unique(df.lv.fit$Sp.id)

# Build the app -----------------------------------------------------------

ui <- fluidPage(
  
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
    column(5, plotOutput("p.l.nit1")),
    column(2, ""),
    column(5, plotOutput("p.l.nit2"))
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
    column(5, plotOutput("p.lv.nit1")),
    column(2, ""),
    column(5, plotOutput("p.lv.nit2"))
  )
  
)

server <- function(input, output, session) {
  
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
    df.l.raw %>% filter(id.number == as.numeric(input$pop.l))
  })
  
  output$sp.text.l <- renderText({                     # Left column: species
    req(sel.row.l())
    paste("Species:", sel.row.l()$Species.name[1])
  })
  
  output$p.l.nit1 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(id.number == input$pop.l) 
    nit.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    t.sel <- ifelse(is.na(nit.l$T.opt.TPC), nit.l$T.max.raw, nit.l$T.opt.TPC)
    
    t.sel2 <- mu.l$temp[which.min(abs(mu.l$temp - t.sel))]
    
    mu.l <- mu.l %>% 
      filter(temp == t.sel2)
    
    curve.l <- tibble::tibble(
      res  = seq(0, 1000, length.out = 200),
      rate = pred_mon(res, r.max = nit.l$r.max.mod, k.s = nit.l$K.s.mod)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.l,
                 aes(x = nit, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Nitrogen Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = nit.l$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*nit.l$K.s.mod/(nit.l$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.l.nit2 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(id.number == input$pop.l) 
    nit.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    t.sel <- ifelse(is.na(nit.l$T.opt.TPC), nit.l$T.max.raw, nit.l$T.opt.TPC)
    
    t.sel2 <- mu.l$temp[which.min(abs(mu.l$temp - t.sel))]
    
    mu.l <- mu.l %>% 
      filter(temp == t.sel2)
    
    curve.l <- tibble::tibble(
      res  = seq(0, 1000, length.out = 200),
      rate = pred_mon(res, r.max = nit.l$r.max.post, k.s = nit.l$K.s.post)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.l,
                 aes(x = nit, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Nitrogen Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = nit.l$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*nit.l$K.s.post/(nit.l$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
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
  
  output$p.lv.nit1 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    nit.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(nit.lv$T.opt.TPC), nit.lv$T.max.raw, nit.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 1000, length.out = 200),
      rate = pred_mon(res, r.max = nit.lv$r.max.mod, k.s = nit.lv$K.s.mod)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.lv,
                 aes(x = nit, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Nitrogen Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = nit.lv$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*nit.lv$K.s.mod/(nit.lv$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.lv.nit2 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    nit.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(nit.lv$T.opt.TPC), nit.lv$T.max.raw, nit.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 1000, length.out = 200),
      rate = pred_mon(res, r.max = nit.lv$r.max.post, k.s = nit.lv$K.s.post)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.lv,
                 aes(x = nit, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Nitrogen Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = nit.lv$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*nit.lv$K.s.post/(nit.lv$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
  })
  
}

shinyApp(ui, server)
