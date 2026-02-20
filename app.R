library(shiny)
library(bslib)
library(shinyvalidate)
library(survival)
library(dplyr)
library(ggplot2)

# Helper: Kaplan-Meier step plot -----------------------------------------

surv_at_time <- function(fit, t) {
  m <- if (is.matrix(fit$surv)) fit$surv else matrix(fit$surv, ncol = 1)
  sapply(seq_len(ncol(m)), function(i) {
    idx <- which(fit$time <= t)
    if (length(idx) == 0L) 1.0 else m[max(idx), i]
  })
}

surv_plot <- function(fit, vline, labels = c("Case", "Ref.")) {
  m <- if (is.matrix(fit$surv)) fit$surv else matrix(fit$surv, ncol = 1)
  surv_data <- do.call(
    rbind,
    lapply(seq_len(ncol(m)), function(i) {
      data.frame(
        time = fit$time,
        surv = m[, i],
        group = factor(labels[i], levels = labels)
      )
    })
  )

  # Dots at vline crossing
  pt_y <- surv_at_time(fit, vline)
  pt_data <- data.frame(
    time = vline,
    surv = pt_y,
    group = factor(labels, levels = labels)
  )

  ggplot(surv_data, aes(x = time, y = surv, color = group, linetype = group)) +
    geom_step(linewidth = 1.0) +
    geom_vline(
      xintercept = vline,
      linetype = "dashed",
      color = "gray55",
      linewidth = 0.5
    ) +
    geom_point(
      data = pt_data,
      aes(fill = group),
      shape = 21,
      size = 4,
      color = "white",
      stroke = 1.8
    ) +
    scale_color_manual(values = c("Case" = "#c0392b", "Ref." = "#2c3e50")) +
    scale_fill_manual(
      values = c("Case" = "#c0392b", "Ref." = "#2c3e50"),
      guide = "none"
    ) +
    scale_linetype_manual(values = c("Case" = "solid", "Ref." = "dotted")) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) {
      paste0(x * 100, "%")
    }) +
    scale_x_continuous(breaks = 0:5) +
    coord_cartesian(xlim = c(0, 5)) +
    labs(
      x = "Time (years)",
      y = "Survival probability",
      color = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      axis.line = element_line(color = "gray60", linewidth = 0.4),
      axis.ticks = element_line(color = "gray60"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 15, 5, 10)
    )
}

# Load models once at startup ---------------------------------------------

rfs_model <- local({
  e <- new.env()
  load("PFS_model.Rdata", envir = e)
  get(ls(e)[1], envir = e)
})

os_model <- local({
  e <- new.env()
  load("OS_model.Rdata", envir = e)
  get(ls(e)[1], envir = e)
})

# UI ----------------------------------------------------------------------

ui <- page_sidebar(
  title = div(
    strong("Adv. CC HALP Survival"),
    tags$small(
      " — Advanced Cervical Cancer",
      class = "text-white-50 fw-normal ms-1"
    )
  ),
  theme = bs_theme(version = 5, primary = "#0e4b8f"),

  tags$style(HTML(
    "
    /* Compact inputs */
    .control-label { font-size: 0.78rem; font-weight: 600; color: #374151; margin-bottom: 1px; }
    .form-control, .form-select { font-size: 0.82rem; padding: 0.2rem 0.5rem; height: auto; }
    .form-group, .shiny-input-container { margin-bottom: 0.35rem !important; }

    /* Compact slider */
    .irs--shiny .irs-bar, .irs--shiny .irs-line { height: 3px; }
    .irs--shiny .irs-handle { width: 16px; height: 16px; top: 22px; }
    .irs { height: 40px; }
    .irs-with-grid { height: 52px; }

    /* Value boxes — tightly hugged */
    .bslib-value-box { min-height: unset !important; }
    .bslib-value-box > .card-body { padding: 0 !important; }
    .bslib-value-box .value-box-grid { padding: 0.3rem 0.5rem !important; gap: 0.25rem !important; min-height: unset !important; }
    .bslib-value-box .value-box-area { padding: 0 !important; gap: 0 !important; justify-content: center; }
    .bslib-value-box .value-box-value { font-size: 1.3rem !important; font-weight: 700; line-height: 1.1; margin: 0 !important; }
    .bslib-value-box .value-box-title { font-size: 0.7rem !important; text-transform: uppercase; letter-spacing: 0.03em; margin: 0 !important; line-height: 1.2; }
    .bslib-value-box .value-box-showcase { font-size: 1.1rem !important; padding: 0.3rem 0.5rem !important; flex: 0 0 auto !important; align-self: center !important; }

    /* Cards */
    .card { box-shadow: 0 2px 10px rgba(0,0,0,0.07) !important; }
    .card-body { padding: 0.6rem !important; }
    .card-header { padding: 0.4rem 0.75rem !important; font-size: 0.88rem; }

    /* HALP badge */
    .halp-badge { border-radius: 8px; font-size: 0.78rem; }

    /* Remove sidebar top gap */
    .sidebar-content { padding-top: 0.25rem !important; }
    .bslib-sidebar-layout > .sidebar > .sidebar-content { padding-top: 0.25rem !important; }

    /* Divider */
    hr { margin: 0.4rem 0; opacity: 0.15; }

    /* Layout gap */
    .bslib-gap-spacing { gap: 0.5rem !important; }
  "
  )),

  sidebar = sidebar(
    width = 220,
    tags$small("HEMATOLOGIC / NUTRITIONAL", class = "text-muted fw-bold"),
    numericInput("hb", "Hemoglobin (g/dL)", 10, step = 0.1, min = 1, max = 15),
    numericInput("alb", "Albumin (g/dL)", 3.5, step = 0.1, min = 1, max = 5),
    numericInput(
      "lymph",
      "Lymphocyte (\u00d710\u2079/L)",
      3,
      step = 0.1,
      min = 0.1,
      max = 5
    ),
    numericInput(
      "plt",
      "Platelets (\u00d710\u2079/L)",
      150,
      step = 10,
      min = 50,
      max = 500
    ),
    uiOutput("halp_display"),
    hr(),
    tags$small("TUMOR CHARACTERISTICS", class = "text-muted fw-bold"),
    numericInput("tsize", "Tumor size (cm)", 4, step = 0.1, min = 1, max = 10),
    selectInput("stage", "Stage", choices = c("I", "II", "III", "IV")),
    selectInput(
      "cmt",
      "Treatment",
      choices = c("RT", "CCRT"),
      selected = "CCRT"
    ),
    selectInput("histo", "Histology", choices = c("SCC", "ADC", "ASC")),
    hr(),
    actionButton(
      "button",
      "Predict!",
      icon = icon("chart-line"),
      class = "btn-primary w-100"
    )
  ),

  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header(
        tagList(icon("rotate-left"), " Recurrence-Free Survival"),
        class = "bg-danger text-white fw-bold"
      ),
      plotOutput("RFS_plot", height = "230px"),
      sliderInput(
        "survyear_RFS",
        "Landmark year",
        value = 1,
        min = 0.1,
        max = 5,
        step = 0.1
      )
    ),
    card(
      card_header(
        tagList(icon("heart-pulse"), " Overall Survival"),
        class = "bg-primary text-white fw-bold"
      ),
      plotOutput("OS_plot", height = "230px"),
      sliderInput(
        "survyear_OS",
        "Landmark year",
        value = 1,
        min = 0.1,
        max = 5,
        step = 0.1
      )
    )
  ),

  layout_columns(
    col_widths = c(3, 3, 3, 3),
    value_box(
      "Case RFS",
      uiOutput("RFS_text_case"),
      showcase = icon("person"),
      theme = "danger"
    ),
    value_box(
      "Baseline RFS",
      uiOutput("RFS_text_ref"),
      showcase = icon("users"),
      theme = "danger"
    ),
    value_box(
      "Case OS",
      uiOutput("OS_text_case"),
      showcase = icon("person"),
      theme = "primary"
    ),
    value_box(
      "Baseline OS",
      uiOutput("OS_text_ref"),
      showcase = icon("users"),
      theme = "primary"
    )
  ),

  tags$p(
    class = "text-muted mt-1 mb-0",
    style = "font-size:0.78rem;",
    icon("circle-info"),
    " Ref. group: high HALP, <4 cm, Stage I, CCRT, SCC",
    tags$span(
      class = "ms-3",
      icon("triangle-exclamation"),
      " For research purpose only."
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  iv <- InputValidator$new()
  iv$add_rule("hb", sv_between(1, 15))
  iv$add_rule("alb", sv_between(1, 5))
  iv$add_rule("lymph", sv_between(0.1, 5))
  iv$add_rule("plt", sv_between(50, 500))
  iv$add_rule("tsize", sv_between(1, 10))
  iv$enable()

  HALP <- eventReactive(input$button, ignoreNULL = FALSE, {
    (input$hb * 0.1) * input$lymph * input$alb / (input$plt / 1e3)
  })

  output$halp_display <- renderUI({
    req(HALP())
    val <- round(HALP(), 2)
    is_high <- val > 22.2
    div(
      class = paste0(
        "alert alert-",
        if (is_high) "success" else "warning",
        " halp-badge py-1 px-2 mt-1 mb-0"
      ),
      div(class = "fw-bold", paste0("HALP: ", val)),
      tags$small(
        if (is_high) "\u2191 High (> 22.2)" else "\u2193 Low (\u2264 22.2)"
      )
    )
  })

  user_df <- eventReactive(input$button, ignoreNULL = FALSE, {
    tsize_cat <- as.character(
      cut(input$tsize, breaks = c(0, 4, 8, Inf), labels = c("<4", "4-8", ">8"))
    )
    histo_cat <- case_when(
      input$histo == "SCC" ~ "SCC",
      input$histo == "ADC" ~ "AD",
      input$histo == "ASC" ~ "ASC"
    )
    data.frame(
      HALP_cutoff = c(ifelse(HALP() <= 22.2, "<= 22.2", "> 22.2"), "> 22.2"),
      stage_cat = c(
        ifelse(input$stage %in% c("I", "II"), "I + II", input$stage),
        "I + II"
      ),
      cmt_cat = c(input$cmt, "CCRT"),
      histo_cat = c(histo_cat, "SCC"),
      tsize_cat = c(tsize_cat, "<4")
    )
  })

  RFS_fit <- reactive({
    req(user_df())
    survfit(rfs_model, newdata = user_df())
  })
  OS_fit <- reactive({
    req(user_df())
    survfit(os_model, newdata = user_df())
  })
  RFS_surv <- reactive({
    req(user_df())
    surv_at_time(RFS_fit(), input$survyear_RFS)
  })
  OS_surv <- reactive({
    req(user_df())
    surv_at_time(OS_fit(), input$survyear_OS)
  })

  output$RFS_plot <- renderPlot(
    {
      req(user_df())
      surv_plot(RFS_fit(), input$survyear_RFS)
    },
    bg = "white"
  )

  output$OS_plot <- renderPlot(
    {
      req(user_df())
      surv_plot(OS_fit(), input$survyear_OS)
    },
    bg = "white"
  )

  fmt_pct <- function(v, idx) {
    if (length(v) >= idx) paste0(round(v[idx] * 100, 2), "%") else "\u2014"
  }

  output$RFS_text_case <- renderUI({
    req(user_df())
    strong(fmt_pct(RFS_surv(), 1))
  })
  output$RFS_text_ref <- renderUI({
    req(user_df())
    strong(fmt_pct(RFS_surv(), 2))
  })
  output$OS_text_case <- renderUI({
    req(user_df())
    strong(fmt_pct(OS_surv(), 1))
  })
  output$OS_text_ref <- renderUI({
    req(user_df())
    strong(fmt_pct(OS_surv(), 2))
  })
}

# Run ---------------------------------------------------------------------

shinyApp(ui, server)
