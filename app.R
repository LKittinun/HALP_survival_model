library(shinydashboard)
library(shinyWidgets)
library(shinyvalidate)
library(patchwork)
library(shiny)
library(survival)
library(survminer)
# UI ----------------------------------------------------------------------

header <- dashboardHeader(title = span(icon("signal", lib = "font-awesome"), 
                                       strong("Adv. CC HALP")),
                          tags$li(class = "dropdown", tags$a(("This application is intended for                                                                     research purpose only.")))
                          )


sidebar <- dashboardSidebar(
      numericInput("hb", "Hemoglobin (g/L)", 10 , step = 0.1, min = 1, max = 15, width = "70%"),
      numericInput("alb", "Albumin (g/L)", 3.5 , min = 1, max = 5, step = 0.1, width = "70%"),
      numericInput("lymph", "Lymphocyte (cells/L)", 3, min = 1, max = 10, step = 0.1, width = "70%"),
      numericInput("plt", "Platelets (cells/L)", 50, min = 50, max = 500, step = 10, width = "70%"),
      numericInput("tsize", "Tumor size (cm)", 4, min = 1, max = 10, step = 0.1, width = "70%"),
      selectInput("stage", "Stage", choices = c("I", "II", "III", "IV")),
      selectInput("cmt", "Modalities", selected = "CCRT", choices = c("RT", "CCRT")),
      selectInput("histo", "Histology",
                  choices = c("Squamous cell carcinoma", "Adenocarcinoma", "Adenosquamous carcinoma")),
      actionButton("button", " Predict! ",icon = icon("table"))
)

body <- dashboardBody(
  tags$style(".box-header h3.box-title {
               font-weight: bold;
               font-size: 20px; }"),
  tags$style(HTML(".box.box-solid.box-primary>.box-header {
                background:#0e4b8f;}
                
             .box.box-solid.box-primary{
                border-bottom-color:#0e4b8f;
                border-left-color:#0e4b8f;
                border-right-color:#0e4b8f;
                border-top-color:#0e4b8f;}
                
                     /* logo */
        .skin-blue .main-header .logo {
                              background-color: #0e4b8f;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #0e4b8f;
                              }        

             ")),
  tags$style(".small-box.bg-blue { background-color: #0e4b8f !important; color: #FFFFFF !important; }"),
  mainPanel(
    fluidRow(
      box(plotly::plotlyOutput("RFS_plot"), title = "Recurrent free survival", 
          status = "danger", solidHeader = TRUE),
      box(plotly::plotlyOutput("OS_plot"), title = "Overall Survival", 
          status = "primary", solidHeader = TRUE)
    ),
    fluidRow(
      column(width = 6,
        valueBoxOutput("RFS_text_case", width = 6),
        valueBoxOutput("RFS_text_ref", width = 6)
      ),
      column(width = 6,
        valueBoxOutput("OS_text_case", width = 6),
        valueBoxOutput("OS_text_ref", width = 6)
        
      )
    ),
    fluidRow(
      column(width = 6,
              box(width = 12, sliderInput("survyear_RFS", NULL, value = 1, min = 0, max = 5, step = 0.1),
                  status = "danger", title = "RFS year adjustment", 
                  solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE)
      ),
      column(width = 6,
             box(width = 12, sliderInput("survyear_OS", NULL, value = 1, min = 0, max = 5, step = 0.1),
                 status = "primary", title = "OS year adjustment", 
                 solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE)
      )

    ),
  h4(strong("Ref. group"), "is based on low HALP, <4 cm, stage I, RT, and, SCC ")
))

ui <- dashboardPage(header, sidebar, body)

# Server ------------------------------------------------------------------

## Preload data

train_RFS_model_shiny <- 
  reactive({get(load("RFS_model.Rdata"))})

train_OS_model_shiny <- 
  reactive({get(load("OS_model.Rdata"))})

server <- function(input, output) {

  #Validator
  iv <- InputValidator$new()
  iv$add_rule("hb", sv_between(1, 15))
  iv$add_rule("alb", sv_between(1, 5))
  iv$add_rule("lymph", sv_between(1, 10))
  iv$add_rule("plt", sv_between(50, 500))
  iv$add_rule("tsize", sv_between(1, 10))
  iv$enable()
  
   HALP <- reactive({
     input$button
     HALP <- isolate((input$hb * 0.1) * input$lymph * input$alb / (input$plt / 10^3))
     HALP
   })

   ## Generate new dataframe
   user_df <- reactive({
    input$button
    tsize <- isolate(as.character(cut(input$tsize, breaks = c(0, 4, 8, Inf), labels = c("<4", "4-8", ">8"))))
    df <- isolate(
      data.frame(
      HALP_cutoff = c(ifelse(HALP() <= 22.2, "<= 22.2", "> 22.2"), "<= 22.2"),
      stage_cat = c(ifelse(input$stage %in% c("I", "II"), "I + II", input$stage), "I + II"),
      cmt_cat = c(input$cmt, "RT"),
      histo_cat = c(dplyr::case_when(input$histo == "Squamous cell carcinoma" ~ "SCC",
                                   input$histo == "Adenocarcinoma" ~ "AD",
                                   input$histo == "Adenosquamous carcinoma" ~ "ASC"), "SCC"),
      tsize_cat = c(tsize, "<4")
    ))
    df
   })

   ## Generate fitted model

    RFS_fit <- reactive({ survfit(train_RFS_model_shiny(), newdata = user_df()) })
    RFS_fit_summ <- reactive({summary(RFS_fit(), time = input$survyear_RFS)})
    
    OS_fit <- reactive({ survfit(train_OS_model_shiny(), newdata = user_df()) })
    OS_fit_summ <- reactive({summary(OS_fit(), time = input$survyear_OS)})
    
  ## Plot recurrent free survival

      output$RFS_plot <- plotly::renderPlotly({
      
      validate(need(RFS_fit(), "Model is not found"))  
      
      RFS_p <- ggsurvplot(RFS_fit(), data = user_df(), 
                          censor.shape=".", censor.size = 0.2, break.time.by = 1,
                          xlim = c(0, 5), linetype = c("solid", "dotted"),
                          title = NULL , legend.labs = c("Case", "Ref"))$plot +
                          geom_vline(xintercept = input$survyear_RFS, linetype = "dashed")
      RFS_p <- plotly::ggplotly(RFS_p)

      }) 
      
      output$RFS_text_case <- renderValueBox({
        valueBox(
          paste0(round(RFS_fit_summ()$surv[1]*100, 2), " % "), 
          paste0("Case's RFS at ", input$survyear_RFS, " year"), color = "red"
        )
      })
      
      output$RFS_text_ref <- renderValueBox({
        valueBox(
          paste0(round(RFS_fit_summ()$surv[2]*100, 2), " % "), 
          paste0("Baseline RFS at ", input$survyear_RFS, " year"), color = "red"
        )
      })
      
  ## Plot overall survival

    output$OS_plot <- plotly::renderPlotly({
      
      validate(need(OS_fit(), "Model is not found"))  
      
      OS_p <- ggsurvplot(OS_fit(), data = user_df(), 
                          censor.shape=".", censor.size = 0.2, break.time.by = 1,
                          xlim = c(0,5), linetype = c("solid", "dotted"),
                          title = NULL, legend.labs = c("Case", "Ref"))$plot +
                          geom_vline(xintercept = input$survyear_OS, linetype = "dashed")
      OS_p <- plotly::ggplotly(OS_p)
      
    }) 
    
    output$OS_text_case <- renderValueBox({
      valueBox(
       paste0(round(OS_fit_summ()$surv[1]*100, 2), " % "), color = "blue",
       paste0("Case's OS at ", input$survyear_OS, " year")
      )
    })

    output$OS_text_ref <- renderValueBox({
      valueBox(
        paste0(round(OS_fit_summ()$surv[2]*100, 2), " % "), color = "blue", 
        paste0("Baseline OS at ", input$survyear_OS, " year")
      )
    })
}

# Run the application -----------------------------------------------------

shinyApp(ui = ui, server = server)
