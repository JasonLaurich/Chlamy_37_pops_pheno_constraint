# J R Laurich

# January 3rd, 2026

# We're going to take a look at my interspecific light Monod models and data to check fits.

# That's the Edwards 2016, Lewington-Pearce 2019, Narwani 2015, and Levasseur 2025 datasets
# Narwani 2015 only reports the actual model parameters, no raw data. 

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)

# Functions -------------------------------------------------

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

# Load & explore the data --------------------------------------------------------

###### Edwards 2016 ######

df.e.raw <- read.csv('data/38_Edwards_2016_raw_data.csv') # Edwards raw data
head(df.e.raw)

df.e.raw <- df.e.raw %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = paste(species, reference, sep = "_")) # So that each entry is treated seperately!

df.e.raw <- df.e.raw %>% 
  filter(irradiance<= 250) # After looking at the light models, insanely high lights are driving down growth rates. We'll cap this out at 250 to improve model fits.

df.e.raw <-df.e.raw %>% 
  rename(mu = growth.rate,
         Species.name = species,
         temp = temperature,
         light = irradiance) %>% 
  mutate(id.number = as.integer(factor(df.e.raw$unique.id)))

df.e.fit <- read.csv('data/47_Edwards_2016_light_monods.csv') # Edwards light Monods
head(df.e.fit)

df.e.fit <- df.e.fit %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = as.integer(factor(Sp.id, levels = sort(unique(Sp.id)))))

pop_order.e <- unique(df.e.fit$Sp.id) # Unique numbers for each species. 

###### Lewington-Pearce 2019 ######

df.l.raw <- read.csv('data/49_Lewington-Pearce_2019_µ_estimates_light.csv') # Lewington-Pearce raw data
head(df.l.raw)

df.l.raw <-df.l.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.l.raw$Sp.id)))

df.l.fit <- read.csv('data/50_Lewington_2019_light_monods.csv') # Lewington light monod fits
head(df.l.fit)

pop_order.l <- unique(df.l.fit$Sp.id)

###### Levasseur 2025 ######

df.lv.raw <- read.csv('data/52_Levasseur_2025_µ_estimates_light.csv') # Levasseur raw data
head(df.lv.raw)

df.lv.raw <-df.lv.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.lv.raw$Sp.id)))

df.lv.fit <- read.csv('data/53_Levasseur_2025_light_monods.csvv') # Levasseur light monod fits
head(df.lv.fit)

pop_order.lv <- unique(df.lv.fit$Sp.id)

# Build the app -----------------------------------------------------------

ui <- fluidPage(
  
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
    column(5, plotOutput("p.e.light1")),
    column(2, ""),
    column(5, plotOutput("p.e.light2"))
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
    column(5, plotOutput("p.l.light1")),
    column(2, ""),
    column(5, plotOutput("p.l.light2"))
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
    column(5, plotOutput("p.lv.light1")),
    column(2, ""),
    column(5, plotOutput("p.lv.light2"))
  )
  
)

server <- function(input, output, session) {
  
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
    df.e.raw %>% filter(id.number == as.numeric(input$pop.e))
  })
  
  output$sp.text.e <- renderText({                     # Left column: species
    req(sel.row.e())
    paste("Species:", sel.row.e()$Species.name[1])
  })
  
  output$study.text.e <- renderText({                    # Right column: study
    req(sel.row.e())
    paste("Study:", sel.row.e()$reference[1])
  })
  
  output$p.e.light1 <- renderPlot({
    
    mu.e <- df.e.raw %>% filter(id.number == input$pop.e) 
    light.e <- df.e.fit %>% filter(unique.id == input$pop.e)
    
    t.sel <- ifelse(is.na(light.e$T.opt.TPC), light.e$T.max.raw, light.e$T.opt.TPC)
    
    t.sel2 <- mu.e$temp[which.min(abs(mu.e$temp - t.sel))]
    
    mu.e <- mu.e %>% 
      filter(temp == t.sel2)
    
    curve.e <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.e$r.max.mod, k.s = light.e$K.s.mod)
    )
    
    ggplot(curve.e, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.e,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Edward 2016: Light Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.e$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.e$K.s.mod/(light.e$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.e.light2 <- renderPlot({
    
    mu.e <- df.e.raw %>% filter(id.number == input$pop.e) 
    light.e <- df.e.fit %>% filter(unique.id == input$pop.e)
    
    t.sel <- ifelse(is.na(light.e$T.opt.TPC), light.e$T.max.raw, light.e$T.opt.TPC)
    
    t.sel2 <- mu.e$temp[which.min(abs(mu.e$temp - t.sel))]
    
    mu.e <- mu.e %>% 
      filter(temp == t.sel2)
    
    curve.e <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.e$r.max.post, k.s = light.e$K.s.post)
    )
    
    ggplot(curve.e, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.e,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Edward 2016: Light Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.e$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.e$K.s.post/(light.e$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
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
    df.l.raw %>% filter(id.number == as.numeric(input$pop.l))
  })
  
  output$sp.text.l <- renderText({                     # Left column: species
    req(sel.row.l())
    paste("Species:", sel.row.l()$Species.name[1])
  })
  
  output$p.l.light1 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(id.number == input$pop.l) 
    light.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    t.sel <- ifelse(is.na(light.l$T.opt.TPC), light.l$T.max.raw, light.l$T.opt.TPC)
    
    t.sel2 <- mu.l$temp[which.min(abs(mu.l$temp - t.sel))]
    
    mu.l <- mu.l %>% 
      filter(temp == t.sel2)
    
    curve.l <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.l$r.max.mod, k.s = light.l$K.s.mod)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.l,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Light Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.l$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.l$K.s.mod/(light.l$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.l.light2 <- renderPlot({
    
    mu.l <- df.l.raw %>% filter(id.number == input$pop.l) 
    light.l <- df.l.fit %>% filter(Sp.id == input$pop.l)
    
    t.sel <- ifelse(is.na(light.l$T.opt.TPC), light.l$T.max.raw, light.l$T.opt.TPC)
    
    t.sel2 <- mu.l$temp[which.min(abs(mu.l$temp - t.sel))]
    
    mu.l <- mu.l %>% 
      filter(temp == t.sel2)
    
    curve.l <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.l$r.max.post, k.s = light.l$K.s.post)
    )
    
    ggplot(curve.l, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.l,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Lewington-Pearce 2019: Light Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.l$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.l$K.s.post/(light.l$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
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
  
  output$p.lv.light1 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    light.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(light.lv$T.opt.TPC), light.lv$T.max.raw, light.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.lv$r.max.mod, k.s = light.lv$K.s.mod)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.lv,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Light Monod (model parameters)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.lv$r.max.mod, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.lv$K.s.mod/(light.lv$r.max.mod - 0.1), linetype = "dashed", colour = "darkred")
  })
  
  output$p.lv.light2 <- renderPlot({
    
    mu.lv <- df.lv.raw %>% filter(id.number == input$pop.lv) 
    light.lv <- df.lv.fit %>% filter(Sp.id == input$pop.lv)
    
    t.sel <- ifelse(is.na(light.lv$T.opt.TPC), light.lv$T.max.raw, light.lv$T.opt.TPC)
    
    t.sel2 <- mu.lv$temp[which.min(abs(mu.lv$temp - t.sel))]
    
    mu.lv <- mu.lv %>% 
      filter(temp == t.sel2)
    
    curve.lv <- tibble::tibble(
      res  = seq(0, 250, length.out = 200),
      rate = pred_mon(res, r.max = light.lv$r.max.post, k.s = light.lv$K.s.post)
    )
    
    ggplot(curve.lv, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "goldenrod1") +
      geom_point(data = mu.lv,
                 aes(x = light, y = mu),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Light level (irradiance)",
        y = "Exponential growth rate",
        title = "Levasseur 2025: Light Monod (posterior medians)"
      ) +
      theme_classic() +
      ylim(-0.1, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      geom_hline(yintercept = light.lv$r.max.post, linetype= "dashed", colour = "forestgreen") +
      geom_vline(xintercept = 0.1*light.lv$K.s.post/(light.lv$r.max.post - 0.1), linetype = "dashed", colour = "darkred")
  })
}

shinyApp(ui, server)
