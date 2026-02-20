library(shiny)
library(bslib)
library(shinyvalidate)
library(survival)
library(dplyr)
library(ggplot2)
library(plotly)

# Helper: Kaplan-Meier step plot (replaces ggsurvplot) --------------------

surv_plot <- function(fit, vline, labels = c("Case", "Ref.")) {
  group_var <- rep(labels, as.integer(fit$strata))
  surv_data <- data.frame(
    time  = fit$time,
    surv  = fit$surv,
    group = factor(group_var, levels = labels)
  )
  p <- ggplot(surv_data, aes(x = time, y = surv, color = group, linetype = group)) +
    geom_step(linewidth = 0.8) +
    geom_vline(xintercept = vline, linetype = "dashed", color = "gray40") +
    scale_y_continuous(limits = c(0, 1),
                       labels = function(x) paste0(x * 100, "%")) +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(values    = c("Case" = "#e74c3c", "Ref." = "#2c3e50")) +
    scale_linetype_manual(values = c("Case" = "solid",   "Ref." = "dotted")) +
    labs(x = "Time (years)", y = "Survival probability",
         color = NULL, linetype = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
  ggplotly(p) |> layout(legend = list(orientation = "h", x = 0.3, y = -0.2))
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
  title = tagList(icon("signal"), strong("Adv. CC HALP Survival")),
  theme = bs_theme(bootswatch = "flatly", primary = "#0e4b8f"),

  sidebar = sidebar(
    tags$small("For research purpose only.", class = "text-muted"),
    hr(),
    numericInput("hb",    "Hemoglobin (g/dL)",          10,  step = 0.1, min = 1,  max = 15),
    numericInput("alb",   "Albumin (g/dL)",               3.5, step = 0.1, min = 1,  max = 5),
    numericInput("lymph", "Lymphocyte (cells/\u00b5L)",   3,   step = 0.1, min = 1,  max = 10),
    numericInput("plt",   "Platelets (\u00d710\u2079/L)", 50,  step = 10,  min = 50, max = 500),
    numericInput("tsize", "Tumor size (cm)",               4,  step = 0.1, min = 1,  max = 10),
    selectInput("stage", "Stage",      choices = c("I", "II", "III", "IV")),
    selectInput("cmt",   "Modalities", choices = c("RT", "CCRT"), selected = "CCRT"),
    selectInput("histo", "Histology",
                choices = c("Squamous cell carcinoma",
                            "Adenocarcinoma",
                            "Adenosquamous carcinoma")),
    actionButton("button", "Predict!", icon = icon("table"),
                 class = "btn-primary w-100 mt-2")
  ),

  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header("Recurrence-Free Survival",
                  class = "bg-danger text-white fw-bold"),
      plotlyOutput("RFS_plot"),
      sliderInput("survyear_RFS", "Year", value = 1, min = 0, max = 5, step = 0.1)
    ),
    card(
      card_header("Overall Survival",
                  class = "bg-primary text-white fw-bold"),
      plotlyOutput("OS_plot"),
      sliderInput("survyear_OS", "Year", value = 1, min = 0, max = 5, step = 0.1)
    )
  ),

  layout_columns(
    col_widths = c(3, 3, 3, 3),
    value_box("Case RFS",     uiOutput("RFS_text_case"), theme = "danger"),
    value_box("Baseline RFS", uiOutput("RFS_text_ref"),  theme = "danger"),
    value_box("Case OS",      uiOutput("OS_text_case"),  theme = "primary"),
    value_box("Baseline OS",  uiOutput("OS_text_ref"),   theme = "primary")
  ),

  p(strong("Ref. group:"), "low HALP, <4 cm, Stage I, RT, SCC")
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  iv <- InputValidator$new()
  iv$add_rule("hb",    sv_between(1,  15))
  iv$add_rule("alb",   sv_between(1,  5))
  iv$add_rule("lymph", sv_between(1,  10))
  iv$add_rule("plt",   sv_between(50, 500))
  iv$add_rule("tsize", sv_between(1,  10))
  iv$enable()

  HALP <- eventReactive(input$button, {
    (input$hb * 0.1) * input$lymph * input$alb / (input$plt / 1e3)
  })

  user_df <- eventReactive(input$button, {
    tsize_cat <- as.character(
      cut(input$tsize, breaks = c(0, 4, 8, Inf), labels = c("<4", "4-8", ">8"))
    )
    histo_cat <- case_when(
      input$histo == "Squamous cell carcinoma"  ~ "SCC",
      input$histo == "Adenocarcinoma"            ~ "AD",
      input$histo == "Adenosquamous carcinoma"   ~ "ASC"
    )
    data.frame(
      HALP_cutoff = c(ifelse(HALP() <= 22.2, "<= 22.2", "> 22.2"), "<= 22.2"),
      stage_cat   = c(ifelse(input$stage %in% c("I", "II"), "I + II", input$stage), "I + II"),
      cmt_cat     = c(input$cmt, "RT"),
      histo_cat   = c(histo_cat, "SCC"),
      tsize_cat   = c(tsize_cat, "<4")
    )
  })

  RFS_fit      <- reactive({ survfit(rfs_model, newdata = user_df()) })
  OS_fit       <- reactive({ survfit(os_model,  newdata = user_df()) })
  RFS_fit_summ <- reactive({ summary(RFS_fit(), time = input$survyear_RFS) })
  OS_fit_summ  <- reactive({ summary(OS_fit(),  time = input$survyear_OS) })

  output$RFS_plot <- renderPlotly({
    req(user_df())
    surv_plot(RFS_fit(), input$survyear_RFS)
  })

  output$OS_plot <- renderPlotly({
    req(user_df())
    surv_plot(OS_fit(), input$survyear_OS)
  })

  surv_pct <- function(summ, idx) {
    s <- tryCatch(summ$surv, error = function(e) NULL)
    if (!is.null(s) && length(s) >= idx) paste0(round(s[idx] * 100, 2), "%") else "\u2014"
  }

  output$RFS_text_case <- renderUI({ strong(surv_pct(RFS_fit_summ(), 1)) })
  output$RFS_text_ref  <- renderUI({ strong(surv_pct(RFS_fit_summ(), 2)) })
  output$OS_text_case  <- renderUI({ strong(surv_pct(OS_fit_summ(),  1)) })
  output$OS_text_ref   <- renderUI({ strong(surv_pct(OS_fit_summ(),  2)) })
}

# Run ---------------------------------------------------------------------

shinyApp(ui, server)
