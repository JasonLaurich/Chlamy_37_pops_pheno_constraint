# app.R — µ estimation comparison: JL (indirect) vs JB (direct Monod)
# Phosphorus gradient only
#
# Data folder needs:
#   14_phosphorous_rfus_time.csv   (from processed-data/)
#   15_µ_estimates_phosphorous.csv (from processed-data/)
#   jb_dir_params.csv              (write.csv(jb.dir.params %>% select(population, umax, ks), ...))
#   jb_exp_slopes.csv              (write.csv(jb.slopes, ...) — see prep code)

library(shiny)
library(tidyverse)
library(minpack.lm)
library(nls.multstart)

col.jb <- "#E8931F"   # orange — JB throughout
col.jl <- "#2166AC"   # blue   — JL throughout

# Load data ---------------------------------------------------------------

df.raw <- read.csv("data/14_phosphorous_rfus_time.csv") %>%
  rename(well.ID = well_plate, phos = phosphate_concentration) %>%
  filter(population != "COMBO", !is.na(RFU)) %>%
  mutate(population = as.character(population),
         well = sub("_.*", "", well.ID),
         logRFU = log(RFU + 0.001)) %>%
  group_by(population) %>%
  mutate(rep = as.integer(factor(well, levels = sort(unique(well))))) %>%
  ungroup()

df.mu <- read.csv("data/15_µ_estimates_phosphorous.csv")
names(df.mu)[names(df.mu) == "µ"] <- "mu"
df.mu <- df.mu %>%
  mutate(population = as.character(population),
         well = sub("_.*", "", well.ID)) %>%
  group_by(population) %>%
  mutate(rep = as.integer(factor(well, levels = sort(unique(well))))) %>%
  ungroup()

df.jb <- read.csv("data/jb_dir_params.csv") %>%
  mutate(population = as.character(population))

df.jb.slopes <- read.csv("data/jb_exp_slopes.csv") %>%
  mutate(population = as.character(population))

# N0_mean per population for JB timeseries curve
n0_means <- df.raw %>%
  group_by(population, well.ID) %>%
  summarise(N0 = RFU[which.min(days)], .groups = "drop") %>%
  group_by(population) %>%
  summarise(N0_mean = mean(N0), .groups = "drop")

df.jb <- df.jb %>% left_join(n0_means, by = "population")

# Population ordering (numeric first, then alphabetical)
pops <- unique(df.raw$population)
pop_order <- c(
  as.character(sort(as.numeric(pops[grepl("^[0-9]+$", pops)]))),
  sort(pops[!grepl("^[0-9]+$", pops)])
)

phos_levels <- sort(unique(df.raw$phos))   # 10 levels
p_seq       <- seq(0, 55, length.out = 200)

# UI ----------------------------------------------------------------------

ui <- fluidPage(

  titlePanel("µ estimation: JL indirect vs JB direct Monod — Phosphorus"),

  fluidRow(
    column(3,
           selectInput("pop", "Population:", choices = pop_order, selected = pop_order[1]),
           actionButton("prev_pop", "← Prev"),
           actionButton("next_pop", "Next →")
    ),
    column(3,
           selectInput("rep", "Replicate:", choices = 1:4, selected = 1),
           actionButton("prev_rep", "← Prev"),
           actionButton("next_rep", "Next →")
    ),
    column(6,
           tags$div(style = "padding-top: 25px; font-size: 14px;",
                    tags$span(style = paste0("color:", col.jl, "; font-weight:bold;"), "— JL (indirect, per well)   "),
                    tags$span(style = paste0("color:", col.jb, "; font-weight:bold;"), "— JB (direct Monod, pooled)   "),
                    tags$span(style = "color: grey30;", "● data points")
           )
    )
  ),

  hr(),

  # Top row: combined Monod curve comparison
  fluidRow(
    column(8, plotOutput("monod_compare", height = "320px")),
    column(4, tags$div(style = "padding-top: 140px; font-size: 13px; line-height: 2;",
                       tags$div(tags$span(style = paste0("color:", col.jl, "; font-weight:bold;"), "● "),
                                "JL µ estimates (this rep)"),
                       tags$div(tags$span(style = paste0("color:", col.jb, "; font-weight:bold;"), "● "),
                                "JB exponential-window slopes (this rep)"),
                       tags$div(tags$span(style = paste0("color:", col.jl, "; font-weight:bold;"), "— "),
                                "JL Monod fit (this rep)"),
                       tags$div(tags$span(style = paste0("color:", col.jb, "; font-weight:bold;"), "— "),
                                "JB Monod fit (pooled)")))
  ),

  hr(),

  # Bottom: RFU timeseries per phos level
  fluidRow(
    column(2, plotOutput("p1",  height = "200px")),
    column(2, plotOutput("p2",  height = "200px")),
    column(2, plotOutput("p3",  height = "200px")),
    column(2, plotOutput("p4",  height = "200px")),
    column(2, plotOutput("p5",  height = "200px")),
    column(2, plotOutput("p6",  height = "200px"))
  ),
  fluidRow(
    column(2, plotOutput("p7",  height = "200px")),
    column(2, plotOutput("p8",  height = "200px")),
    column(2, plotOutput("p9",  height = "200px")),
    column(2, plotOutput("p10", height = "200px"))
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  observeEvent(input$next_pop, {
    i <- match(input$pop, pop_order)
    updateSelectInput(session, "pop", selected = pop_order[ifelse(i == length(pop_order), 1, i + 1)])
  })
  observeEvent(input$prev_pop, {
    i <- match(input$pop, pop_order)
    updateSelectInput(session, "pop", selected = pop_order[ifelse(i == 1, length(pop_order), i - 1)])
  })
  observeEvent(input$next_rep, {
    j <- as.integer(input$rep)
    updateSelectInput(session, "rep", selected = ifelse(j == 4, 1, j + 1))
  })
  observeEvent(input$prev_rep, {
    j <- as.integer(input$rep)
    updateSelectInput(session, "rep", selected = ifelse(j == 1, 4, j - 1))
  })

  # Reactive: JL Monod fit for selected pop + rep (fitted to this rep's µ estimates)
  jl_fit <- reactive({
    d <- df.mu %>%
      filter(population == input$pop,
             rep        == as.integer(input$rep),
             mu > 0, is.finite(mu))
    if (nrow(d) < 3) return(NULL)
    tryCatch(
      nlsLM(mu ~ umax * (phos / (ks + phos)),
            data    = d,
            start   = c(umax = 1, ks = 1),
            lower   = c(umax = 0, ks = 0),
            upper   = c(umax = 3, ks = 100),
            control = nls.control(maxiter = 1024, minFactor = 1/204800000)),
      error = function(e) NULL
    )
  })

  # Top: combined Monod curve comparison
  output$monod_compare <- renderPlot({

    d.jb.pts <- df.jb.slopes %>%
      filter(population == input$pop, rep == as.integer(input$rep))

    d.jl.pts <- df.mu %>%
      filter(population == input$pop, rep == as.integer(input$rep))

    jb_row <- df.jb %>% filter(population == input$pop)
    fit     <- jl_fit()

    jb_curve <- if (nrow(jb_row) == 1)
      data.frame(phos = p_seq, mu = jb_row$umax * p_seq / (jb_row$ks + p_seq), method = "JB")
    else NULL

    jl_curve <- if (!is.null(fit)) {
      umax_jl <- coef(fit)["umax"]; ks_jl <- coef(fit)["ks"]
      data.frame(phos = p_seq, mu = umax_jl * p_seq / (ks_jl + p_seq), method = "JL")
    } else NULL

    curves <- bind_rows(jb_curve, jl_curve)

    p <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      labs(title = "Monod curves: JB (pooled direct) vs JL (indirect, this rep)",
           x = "Phosphate (µM P)", y = "Growth rate (day⁻¹)") +
      theme_minimal(base_size = 13)

    if (nrow(d.jb.pts) > 0)
      p <- p + geom_point(data = d.jb.pts, aes(x = phos, y = slope),
                          color = col.jb, size = 3, alpha = 0.8)

    if (nrow(d.jl.pts) > 0)
      p <- p + geom_point(data = d.jl.pts, aes(x = phos, y = mu),
                          color = col.jl, size = 3, alpha = 0.8)

    if (nrow(curves) > 0)
      p <- p + geom_line(data = curves, aes(x = phos, y = mu, color = method),
                         linewidth = 1.4) +
        scale_color_manual(values = c("JB" = col.jb, "JL" = col.jl), guide = "none")
    p
  })

  # Bottom: RFU timeseries per phos level with both curves
  lapply(1:10, function(idx) {
    local({
      i <- idx
      output[[paste0("p", i)]] <- renderPlot({

        p_level <- phos_levels[i]

        d <- df.raw %>%
          filter(population == input$pop,
                 rep        == as.integer(input$rep),
                 phos       == p_level) %>%
          arrange(days)

        if (nrow(d) == 0)
          return(ggplot() + labs(title = paste(p_level, "µM P")) + theme_minimal())

        # --- JL curve: re-run sliding window + nls_multstart (your exact approach) ---
        t.series  <- unique(d$days)
        ln.slopes <- c()
        for (z in t.series[2:length(t.series)]) {
          d.sl      <- d[d$days <= z, ]
          ln_slope  <- lm(logRFU ~ days, data = d.sl)
          ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2, 1])
        }
        s    <- max(2, which.max(ln.slopes))
        d.th <- d[d$days <= t.series[s + 1], ]
        d.th$N0 <- d$RFU[1]

        jl_fit_ts <- tryCatch(
          nls_multstart(RFU ~ N0 * exp(r * days),
                        data        = d.th,
                        start_lower = c(r = -4.5),
                        start_upper = c(r =  4.5),
                        iter        = 500,
                        supp_errors = "Y",
                        control     = nls.control(maxiter = 200)),
          error = function(e) NULL
        )

        jl_curve <- if (!is.null(jl_fit_ts)) {
          smt <- seq(min(d.th$days), max(d.th$days), length.out = 100)
          data.frame(days = smt,
                     RFU  = predict(jl_fit_ts, newdata = data.frame(days = smt)),
                     method = "JL")
        } else NULL

        # --- JB curve: Monod exponential over JB exponential window only ---
        jb_row    <- df.jb %>% filter(population == input$pop)
        jb_window <- df.jb.slopes %>%
          filter(population == input$pop,
                 rep        == as.integer(input$rep),
                 phos       == p_level)

        d.jb.th <- if (nrow(jb_window) == 1)
          d %>% arrange(days) %>% slice(1:jb_window$n_pts)
        else
          d %>% arrange(days) %>% slice(1:3)   # fallback

        jb_curve <- if (nrow(jb_row) == 1 && nrow(d.jb.th) > 0) {
          t_jb  <- seq(min(d.jb.th$days), max(d.jb.th$days), length.out = 100)
          mu_jb <- jb_row$umax * p_level / (jb_row$ks + p_level)
          data.frame(days = t_jb,
                     RFU  = jb_row$N0_mean * exp(mu_jb * t_jb),
                     method = "JB")
        } else NULL

        curves <- bind_rows(jl_curve, jb_curve)

        p <- ggplot(d, aes(x = days, y = RFU)) +
          geom_point(size = 2, color = "grey50") +
          labs(title = paste(p_level, "µM P"), x = "Days", y = "RFU") +
          theme_minimal(base_size = 10)

        # Highlight JB's exponential window points in orange
        if (nrow(d.jb.th) > 0)
          p <- p + geom_point(data = d.jb.th, aes(x = days, y = RFU),
                              color = col.jb, size = 2.5)

        # Highlight JL's fitted window in blue
        if (!is.null(jl_fit_ts))
          p <- p + geom_point(data = d.th, aes(x = days, y = RFU),
                              color = col.jl, size = 2.5)

        if (!is.null(curves) && nrow(curves) > 0)
          p <- p +
            geom_line(data = curves, aes(x = days, y = RFU, color = method),
                      linewidth = 1.2) +
            scale_color_manual(values = c("JL" = col.jl, "JB" = col.jb),
                               guide = "none")
        p
      })
    })
  })
}

shinyApp(ui, server)
