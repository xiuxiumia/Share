# 1. ---- Load the Packages ----
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)
library(dda)
library(moments)
library(shinyjs)
library(rmarkdown)
library(shinyBS)

# 2. --- User Interface (UI) ----
ui <- navbarPage(
  # Set up the title for the app
  title = "Direction Dependence Analysis (DDA)",

  # First tab - Upload and Read Data
  tabPanel(
    title = "Upload Data", # Title for first tab
    sidebarLayout( # Side bar for first tab
      sidebarPanel(
        fileInput(
          inputId = "document",
          label = "Upload A CSV File",
          multiple = FALSE,
          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
          buttonLabel = "Browse",
          placeholder = "No file selected"
        )
      ),
      mainPanel( # Main Panel for first tab
        conditionalPanel(
          condition = "output.showMain1 == true",
          h2("Preview of Uploaded Data"),
          DTOutput("tablePreview")
        )
      )
    )
  ),

  # Second tab - Descriptive Statistics and Visualization
  tabPanel(
    title = "Descriptive Statistics", # Title for first tab
    sidebarLayout( # Side bar for second tab
      sidebarPanel(
        selectInput( # Select input cause
          inputId = "cause",
          label = "Select Cause",
          choices = NULL,
          multiple = FALSE
        ),
        selectInput( # Select input effect
          inputId = "effect",
          label = "Select Effect",
          choices = NULL,
          multiple = FALSE
        )
      ),
      mainPanel( # Main Panel for second tab
        tabsetPanel(
          tabPanel(title = "Summary Tables", dataTableOutput("descTable_mainvars")),
          tabPanel(title = "Visualiazation", plotOutput("hist"))
        )
      )
    )
  ),

  # Third tab - Run DDA
  tabPanel(
    title = "DDA Summary", # Title for third tab
    useShinyjs(), # To run .Rmd file
    sidebarLayout( # Side bar for third tab
      sidebarPanel(
        h2("Select Covariates"), # Select covariates section
        selectInput(
          inputId = "Covar_n_dda",
          label = "Select Numerical Covariates",
          choices = NULL,
          multiple = TRUE
        ),
        selectInput(
          inputId = "Covar_c_dda",
          label = "Select Categorical Covariates",
          choices = NULL,
          multiple = TRUE
        ),
        br(),
        h2("DDA - Variables"), # DDA variables conditions
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Number of Bootstrap Samples",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_boot1",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        numericInput(
          "BootN_v_dda",
          label = NULL,
          value = 200
        ),
        bsPopover(
          id = "help_boot1",
          title = NULL,
          content = "The bootstrap sample size should be greater than 1.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "BootType_v_dda",
          label = "Bootstrap Type",
          choices = c("perc", "bca"),
          selected = "perc"
        ),
        numericInput(
          inputId = "CI_v_dda",
          label = "Bootstrap Confidence Intervals",
          value = .95,
          min = 0,
          max = 1
        ),
        br(),
        h2("DDA - Residuals"), # DDA residuals conditions
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Number of Bootstrap Samples",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_boot2",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        numericInput(
          "BootN_r_dda",
          label = NULL,
          value = 200
        ),
        bsPopover(
          id = "help_boot2",
          title = NULL,
          content = "The bootstrap sample size should be greater than 1.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "BootType_r_dda",
          label = "Bootstrap Type",
          choices = c("perc", "bca"),
          selected = "perc"
        ),
        numericInput(
          inputId = "CI_r_dda",
          label = "Bootstrap Confidence Intervals",
          value = .95,
          min = 0,
          max = 1
        ),
        selectInput(
          inputId = "ProbTrans_r_dda",
          label = "Probability Integral Transformation",
          choices = c("TRUE", "FASLE"),
          selected = "FALSE"
        ),
        br(),
        h2("DDA - Independence"), # DDA independence conditions
        numericInput(
          inputId = "NlFun_i_dda",
          label = "Values Used for Power Transformation",
          value = NA
        ),
        selectInput(
          inputId = "Hetero_i_dda",
          label = "Separate Homoscedasticity Tests",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        selectInput(
          inputId = "HsicMethod_i_dda",
          label = "HSCI Inference Method",
          choices = c("gamma", "eigenvalue", "boot", "permutation"),
          selected = "gamma"
        ),
        selectInput(
          inputId = "DiffTest_i_dda",
          label = "Differences Tests",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        selectInput(
          inputId = "Paral_i_dda",
          label = "Multiple Cores",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        numericInput(
          inputId = "CoresN_i_dda",
          label = "Number of Cores",
          value = 1
        ),
        tags$div(
          id = "boot_section_i_dda",
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Number of Bootstrap Samples",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_boot3",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        numericInput(
          "BootN_i_dda",
          label = NULL,
          value = 200
        ),
        bsPopover(
          id = "help_boot3",
          title = NULL,
          content = "The bootstrap sample size should be greater than 1.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "BootType_i_dda",
          label = "Bootstrap Type",
          choices = c("perc", "bca"),
          selected = "perc"
        ),
        numericInput(
          inputId = "CI_i_dda",
          label = "Bootstrap Confidence Intervals",
          value = .95,
          min = 0,
          max = 1
        ),
        br(),
        downloadButton("downloadReport1", "Export Report", icon = shiny::icon("file-pdf")) # Download DDA Report
      ),
      mainPanel( # Main Panel for third tab
        tabsetPanel(
          tabPanel(title = "Variables", verbatimTextOutput("dda_var")),
          tabPanel(title = "Residuals", verbatimTextOutput("dda_res")),
          tabPanel(title = "Independence", verbatimTextOutput("dda_ind"))
        )
      )
    )
  ),

  # Fouth tab - Run CDDA
  tabPanel(
    title = "CDDA Summary", # Title for fourth tab
    useShinyjs(), # To run .Rmd file
    sidebarLayout( # Side bar for fourth tab
      sidebarPanel(
        h2("Select Moderator"), # Select moderator section
        selectInput(
          inputId = "mod_cdda",
          label = "",
          choices = NULL,
          multiple = FALSE
        ),
        checkboxInput(inputId = "Cat_Mod", label = "Categorical Moderator", value = FALSE), # Select categorical moderator
        br(),
        h2("Select Covariates"), # Select covariates section
        selectInput(
          inputId = "Covar_n_cdda",
          label = "Select Numerical Covariates",
          choices = NULL,
          multiple = TRUE
        ),
        selectInput(
          inputId = "Covar_c_cdda",
          label = "Select Categorical Covariates",
          choices = NULL,
          multiple = TRUE
        ),
        br(),
        h2("CDDA - Variables"), # CDDA variables conditions
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Enter A Moderator Value",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_mod1",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        textInput(
          inputId = "Modval_v_cdda",
          label = NULL,
          placeholder = NULL
        ),
        bsPopover(
          id = "help_mod1",
          title = NULL,
          content = "Characters or a numeric sequence specifying the moderator values used in post-hoc probing. Possible characters include mean, median, and JN.",
          placement = "right",
          trigger = "hover"
        ),
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Number of Bootstrap Samples",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_boot4",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        numericInput(
          "BootN_v_cdda",
          label = NULL,
          value = 200
        ),
        bsPopover(
          id = "help_boot4",
          title = NULL,
          content = "The bootstrap sample size should be greater than 1.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "BootType_v_cdda",
          label = "Bootstrap Type",
          choices = c("perc", "bca"),
          selected = "perc"
        ),
        numericInput(
          inputId = "CI_v_cdda",
          label = "Bootstrap Confidence Intervals",
          value = .95,
          min = 0,
          max = 1
        ),
        br(),
        h2("CDDA - Independence"), # CDDA independence conditions
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Enter A Moderator Value",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_mod2",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        textInput(
          inputId = "Modval_i_cdda",
          label = NULL,
          placeholder = NULL
        ),
        bsPopover(
          id = "help_mod2",
          title = NULL,
          content = "Characters or a numeric sequence specifying the moderator values used in post-hoc probing. Possible characters include mean, median, and JN.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "Hetero_i_cdda",
          label = "Separate Homoscedasticity Tests",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        selectInput(
          inputId = "DiffTest_i_cdda",
          label = "Differences Tests",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        selectInput(
          inputId = "Paral_i_cdda",
          label = "Multiple Cores",
          choices = c("TRUE", "FALSE"),
          selected = "FALSE"
        ),
        numericInput(
          inputId = "CoresN_i_cdda",
          label = "Number of Cores",
          value = 1
        ),
        numericInput(
          inputId = "NlFun_i_cdda",
          label = "Values Used for Power Transformation",
          value = NA
        ),
        selectInput(
          inputId = "HsicMethod_i_cdda",
          label = "HSCI Inference Method",
          choices = c("gamma", "eigenvalue", "boot", "permutation"),
          selected = "gamma"
        ),
        tags$div(
          id = "boot_section_i_cdda",
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          tags$label(
            "Number of Bootstrap Samples",
            style = "margin: 0; font-weight: 600; font-size: 15px;"
          ),
          tags$span(
            "?",
            id = "help_boot5",
            style = "
      display: inline-block;
      margin-left: 6px;
      width: 18px;
      height: 18px;
      line-height: 18px;
      text-align: center;
      border-radius: 50%;
      background-color: #007BFF;
      color: white;
      font-size: 12px;
      cursor: pointer;
    "
          )
        ),
        numericInput(
          "BootN_i_cdda",
          label = NULL,
          value = 200
        ),
        bsPopover(
          id = "help_boot5",
          title = NULL,
          content = "The bootstrap sample size should be greater than 1.",
          placement = "right",
          trigger = "hover"
        ),
        selectInput(
          inputId = "BootType_i_cdda",
          label = "Bootstrap Type",
          choices = c("perc", "bca"),
          selected = "perc"
        ),
        numericInput(
          inputId = "CI_i_cdda",
          label = "Bootstrap Confidence Intervals",
          value = .95,
          min = 0,
          max = 1
        ),
        br(),
        downloadButton("downloadReport2", "Export Report", icon = shiny::icon("file-pdf")) # Download CDDA Report
      ),
      mainPanel( # Main Panel for fourth tab
        tabsetPanel(
          tabPanel(
            title = "Variables",
            verbatimTextOutput("cdda_var"),
            uiOutput("plotControls1"),
            plotOutput("cdda_var_plot"),
            uiOutput("summaryControls1"),
            verbatimTextOutput("cdda_var_summary")
          ),
          tabPanel(
            title = "Independence",
            verbatimTextOutput("cdda_ind"),
            uiOutput("plotControls2"),
            conditionalPanel( # Only diff test = true, the plots show
              condition = "input.DiffTest_i_cdda == 'TRUE'",
              plotOutput("cdda_ind_plot")
            ),
            uiOutput("summaryControls2"),
            verbatimTextOutput("cdda_ind_summary")
          )
        )
      )
    )
  )
)

# 3. --- Define Server Logic ----
server <- function(input, output, session) {
  # --- First Tab - Upload and Read Data  ----

  # --- Reactive: Uploaded and Read Data
  data <- reactive({
    req(input$document)
    read.csv(input$document$datapath)
  })

  # --- Reactive: Conditional Output
  output$showMain1 <- reactive({
    !is.null(input$document)
  })
  outputOptions(output, "showMain1", suspendWhenHidden = FALSE)

  # --- renderDT: Show Data
  output$tablePreview <- renderDT({
    req(data())
    datatable(head(data(), 10))
  })

  # --- Second Tab - Descriptive Statistics and Visualization ----

  # --- renderDataTable: Summary Statistics Table
  output$descTable_mainvars <- renderDataTable({
    req(data(), input$cause, input$effect)
    selected_data <- data()[, c(input$cause, input$effect), drop = FALSE]
    desc <- psych::describe(selected_data)
    desc_out <- as.data.frame(desc)
    desc_out <- desc_out[, c("n", "mean", "sd", "min", "max", "skew", "kurtosis")]
    num_cols <- c("mean", "sd", "min", "max", "skew", "kurtosis")
    desc_out[num_cols] <- lapply(desc_out[num_cols], function(x) {
      round(x, 2)
    })
    desc_out <- cbind(Variable = rownames(desc_out), desc_out)
    rownames(desc_out) <- NULL
    desc_out
  })

  # --- renderPlot: Visualization
  output$hist <- renderPlot({
    req(data(), input$cause, input$effect)

    # First histogram
    p1 <- ggplot(data(), aes(x = .data[[input$cause]])) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 30,
        fill = "blue",
        color = "black"
      ) +
      geom_density(color = "red", linewidth = 1.2) +
      labs(
        title = paste("Histogram of", input$cause),
        x = input$cause,
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 18,
          face = "bold"
        ),
        panel.grid = element_blank()
      )

    # Second histogram
    p2 <- ggplot(data(), aes(x = .data[[input$effect]])) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 30,
        fill = "red",
        color = "black"
      ) +
      geom_density(color = "blue", linewidth = 1.2) +
      labs(
        title = paste("Histogram of", input$effect),
        x = input$effect,
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 18,
          face = "bold"
        ),
        panel.grid = element_blank()
      )

    # Display side by side
    grid.arrange(p1, p2, ncol = 2)
  })

  # --- Third Tab - Run DDA ----

  # --- Observe: Select Variables
  observe({
    req(data())

    # 1. Collect all chosen variables
    chosen_vars <- c(
      input$cause,
      input$effect,
      input$Covar_n_dda,
      input$Covar_c_dda
    )
    chosen_vars <- chosen_vars[chosen_vars != ""] # remove blanks
    chosen_vars <- chosen_vars[!is.null(chosen_vars)]

    # 2. Get all variable names from dataset
    all_vars <- names(data())

    # 3. Update one selectInput with remaining variables
    updateChoices <- function(id, current_selection) {
      remaining <- setdiff(all_vars, setdiff(chosen_vars, current_selection))
      updateSelectInput(
        session,
        inputId = id,
        choices = remaining,
        selected = current_selection
      )
    }

    # 4. Update all dropdowns
    updateChoices("cause", input$cause)
    updateChoices("effect", input$effect)
    updateChoices("Covar_n_dda", input$Covar_n_dda)
    updateChoices("Covar_c_dda", input$Covar_c_dda)
  })

  # --- ObserveEvent: Show/Hide Conditions
  observeEvent(list(
    input$HsicMethod_i_dda,
    input$DiffTest_i_dda,
    input$BootType_i_dda,
    input$CI_i_dda
  ), {
    if (input$HsicMethod_i_dda %in% c("boot", "permutation") || input$DiffTest_i_dda == "TRUE") {
      shinyjs::show("BootN_i_dda") # Make selectable
      updateNumericInput(session, "BootN_i_dda", value = 100)

      shinyjs::show("BootType_i_dda") # Make selectable
      updateSelectInput(session, "BootType_i_dda", selected = "perc")

      shinyjs::show("CI_i_dda") # Make selectable
      updateNumericInput(session, "CI_i_dda", value = "0.95")

      shinyjs::show("help_boot3")
      shinyjs::show("boot_section_i_dda")
      
    } else {
      shinyjs::hide("BootN_i_dda") # Make unclickable
      updateNumericInput(session, "BootN_i_dda", value = NA)

      shinyjs::hide("BootType_i_dda") # Make unclickable
      updateSelectInput(session, "BootType_i_dda", selected = NULL)

      shinyjs::hide("CI_i_dda") # Make unclickable
      updateNumericInput(session, "CI_i_dda", value = NA)

      shinyjs::hide("help_boot3")
      shinyjs::hide("boot_section_i_dda")
    }
  })

  observeEvent(input$DiffTest_i_dda, {
    if (input$DiffTest_i_dda == "TRUE") {
      shinyjs::show("Paral_i_dda") # Make selectable
    } else {
      shinyjs::hide("Paral_i_dda") # Make unclickable
      updateSelectInput(session, "Paral_i_dda", selected = "FALSE")

      shinyjs::hide("CoresN_i_dda") # Make unclickable
      updateNumericInput(session, "CoresN_i_dda", value = 1)
    }
  })

  observeEvent(input$Paral_i_dda, {
    if (input$Paral_i_dda == "TRUE") {
      shinyjs::show("CoresN_i_dda") # Make selectable
    } else {
      shinyjs::hide("CoresN_i_dda") # Make unclickable
      updateNumericInput(session, "CoresN_i_dda", value = 1)
    }
  })

  # --- Reactive: Run DDA - Variables
  rundda_var <- reactive({
    req(data(), input$cause, input$effect, input$BootN_v_dda, input$BootType_v_dda, input$CI_v_dda)

    # --- Prepare Formula ----
    # 1. Combine numeric covariates if not NULL
    if (!is.null(input$Covar_n_dda) && length(input$Covar_n_dda) > 0) {
      covar_num <- paste(input$Covar_n_dda, collapse = " + ")
    } else {
      covar_num <- NULL
    }

    # 2. Combine categorical covariates if not NULL
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      covar_cat <- paste(input$Covar_c_dda, collapse = " + ")
      formatted_vars <- paste('"', input$Covar_c_dda, '"', collapse = ", ")
    } else {
      covar_cat <- NULL
      formatted_vars <- NULL
    }

    # 3. Combine numeric and categorical covariates if not NULL
    covars <- paste(c(covar_num, covar_cat), collapse = " + ")

    covars <- if (is.null(covars)) {
      NULL
    } else {
      covars
    }

    # 4. Form fomula
    if (is.null(input$Covar_c_dda) && is.null(input$Covar_n_dda)) {
      formula_str <- paste(input$effect, "~", input$cause)
    } else {
      formula_str <- paste(input$effect, "~", paste(input$cause, "+", covars))
    }

    formula <- as.formula(formula_str)

    # --- Change Categorical Covariates to Factors ----
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      df <- data()

      for (v in input$Covar_c_dda) {
        if (v %in% names(df)) {
          df[[v]] <- as.factor(df[[v]])
        }
      }
    } else {
      df <- data()
    }

    # --- DDA Variables Condition ----
    if (is.na(input$BootN_v_dda) || input$BootN_v_dda == 1 || is.na(input$CI_v_dda) || is.null(input$BootType_v_dda)) {
      return(NULL)
    }

    out_var <- dda.vardist(
      formula,
      pred = input$cause,
      B = input$BootN_v_dda,
      boot.type = input$BootType_v_dda,
      conf.level = input$CI_v_dda,
      data = df
    )
    return(out_var)
  })

  # --- renderPrint: Show DDA - Variables Output
  output$dda_var <- renderPrint({
    req(rundda_var())
    print(rundda_var())
  })

  # --- Reactive: Run DDA - Residuals
  rundda_res <- reactive({
    req(data(), input$cause, input$effect, input$BootN_r_dda, input$BootType_r_dda, input$CI_r_dda)

    # --- Prepare Formula ----
    # 1. Combine numeric covariates if not NULL
    if (!is.null(input$Covar_n_dda) && length(input$Covar_n_dda) > 0) {
      covar_num <- paste(input$Covar_n_dda, collapse = " + ")
    } else {
      covar_num <- NULL
    }

    # 2. Combine categorical covariates if not NULL
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      covar_cat <- paste(input$Covar_c_dda, collapse = " + ")
    } else {
      covar_cat <- NULL
    }

    # 3. Combine numeric and categorical covariates if not NULL
    covars <- paste(c(covar_num, covar_cat), collapse = " + ")

    covars <- if (is.null(covars)) {
      NULL
    } else {
      covars
    }

    # 4. Form fomula
    if (is.null(input$Covar_c_dda) && is.null(input$Covar_n_dda)) {
      formula_str <- paste(input$effect, "~", input$cause)
    } else {
      formula_str <- paste(input$effect, "~", paste(input$cause, "+", covars))
    }

    formula <- as.formula(formula_str)

    # --- Change Categorical Covariates to Factors ----
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      df <- data()

      for (v in input$Covar_c_dda) {
        if (v %in% names(df)) {
          df[[v]] <- as.factor(df[[v]])
        }
      }
    } else {
      df <- data()
    }

    # --- DDA Residuals Condition ----
    if (is.na(input$BootN_r_dda) || input$BootN_r_dda == 1 || is.na(input$CI_r_dda) || is.null(input$BootType_r_dda)) {
      return(NULL)
    }

    out_res <- dda.resdist(
      formula,
      pred = input$cause,
      B = input$BootN_r_dda,
      boot.type = input$BootType_r_dda,
      conf.level = input$CI_r_dda,
      prob.trans = paste0('"', input$ProbTrans_r_dda, "'"),
      data = df
    )

    return(out_res)
  })

  # --- renderPrint: Show DDA - Residuals Output
  output$dda_res <- renderPrint({
    req(rundda_res())
    print(rundda_res())
  })

  # --- Reactive: Run DDA - Independence
  rundda_ind <- reactive({
    validate(
      need(data(), "Please upload or select a dataset."),
      need(input$cause, "Please select a cause variable."),
      need(input$effect, "Please select an effect variable."),
      need(input$NlFun_i_dda, "Please choose a value used for power transformation.")
    )

    req(
      data(), input$cause, input$effect, input$NlFun_i_dda
    )

    # --- Prepare Formula ----
    # 1. Combine numeric covariates if not NULL
    if (!is.null(input$Covar_n_dda) && length(input$Covar_n_dda) > 0) {
      covar_num <- paste(input$Covar_n_dda, collapse = " + ")
    } else {
      covar_num <- NULL
    }

    # 2. Combine categorical covariates if not NULL
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      covar_cat <- paste(input$Covar_c_dda, collapse = " + ")
    } else {
      covar_cat <- NULL
    }

    # 3. Combine numeric and categorical covariates if not NULL
    covars <- paste(c(covar_num, covar_cat), collapse = " + ")

    covars <- if (is.null(covars)) {
      NULL
    } else {
      covars
    }

    # 4. Form fomula
    if (is.null(input$Covar_c_dda) && is.null(input$Covar_n_dda)) {
      formula_str <- paste(input$effect, "~", input$cause)
    } else {
      formula_str <- paste(input$effect, "~", paste(input$cause, "+", covars))
    }

    formula <- as.formula(formula_str)

    # --- Change Categorical Covariates to Factors ----
    if (!is.null(input$Covar_c_dda) && length(input$Covar_c_dda) > 0) {
      df <- data()

      for (v in input$Covar_c_dda) {
        if (v %in% names(df)) {
          df[[v]] <- as.factor(df[[v]])
        }
      }
    } else {
      df <- data()
    }

    # --- DDA Independence Condition ----
    args <- list(
      formula = formula,
      pred = input$cause,
      hsic.method = input$HsicMethod_i_dda,
      nlfun = input$NlFun_i_dda,
      hetero = input$Hetero_i_dda,
      diff = as.logical(input$DiffTest_i_dda),
      data = df
    )

    args$parallelize <- ifelse(input$DiffTest_i_dda == "TRUE", FALSE, TRUE)
    args$cores <- ifelse(input$Paral_i_cdda == "FALSE", 0, input$CoresN_i_cdda)

    if ((input$HsicMethod_i_dda %in% c("boot", "permutation") || input$DiffTest_i_dda == "TRUE") &&
      (is.na(input$BootN_i_dda) || input$BootN_i_dda == 1 || is.na(input$CI_i_dda) || is.null(input$BootType_i_dda))) {
      return(NULL)
    }

    if (input$HsicMethod_i_dda %in% c("boot", "permutation") || input$DiffTest_i_dda == "TRUE") {
      if (!is.na(input$BootN_i_dda)) args$B <- input$BootN_i_dda
      if (!is.null(input$BootType_i_dda)) args$boot.type <- input$BootType_i_dda
      if (!is.na(input$CI_i_dda)) args$conf.level <- input$CI_i_dda
    }

    out_ind <- do.call(dda.indep, args)

    return(out_ind)
  })

  # --- renderPrint: Show DDA - Independent Output
  output$dda_ind <- renderPrint({
    req(rundda_ind())
    print(rundda_ind())
  })

  # --- downloadHandler: Download report
  output$downloadReport1 <- downloadHandler(
    filename = "DDA Report.pdf",
    content = function(file) {
      src <- normalizePath("DDA_Report.Rmd")

      owd <- setwd(tempdir()) # temporarily dir
      on.exit(setwd(owd))
      file.copy(src, "DDA_Report.Rmd", overwrite = TRUE)

      out <- render("DDA_Report.Rmd",
        output_format = pdf_document(),
        params = list(
          rundda_var = rundda_var(),
          rundda_res = rundda_res(),
          rundda_ind = rundda_ind()
        ),
        envir = new.env(parent = globalenv())
      )

      file.rename(out, file)
    }
  )

  # --- Fourth Tab - Run CDDA ----

  # --- Observe: Select Variables
  observe({
    req(data())

    # 1. Collect all chosen variables
    chosen_vars <- c(
      input$cause,
      input$effect,
      input$mod_cdda,
      input$Covar_n_cdda,
      input$Covar_c_cdda
    )
    chosen_vars <- chosen_vars[chosen_vars != ""] # remove blanks
    chosen_vars <- chosen_vars[!is.null(chosen_vars)]

    # 2. Get all variable names from dataset
    all_vars <- names(data())

    # 3. Update one selectInput with remaining variables
    updateChoices <- function(id, current_selection) {
      remaining <- setdiff(all_vars, setdiff(chosen_vars, current_selection))
      updateSelectInput(
        session,
        inputId = id,
        choices = remaining,
        selected = current_selection
      )
    }

    # 4. Update all dropdowns
    updateChoices("cause", input$cause)
    updateChoices("effect", input$effect)
    updateChoices("mod_cdda", input$mod_cdda)
    updateChoices("Covar_n_cdda", input$Covar_n_cdda)
    updateChoices("Covar_c_cdda", input$Covar_c_cdda)
  })

  # --- ObserveEvent: Show/Hide Conditions
  observeEvent(list(
    input$HsicMethod_i_cdda,
    input$DiffTest_i_cdda,
    input$BootType_i_cdda,
    input$CI_i_cdda
  ), {
    if (input$HsicMethod_i_cdda %in% c("boot", "permutation") || input$DiffTest_i_cdda == "TRUE") {
      shinyjs::show("BootN_i_cdda") # Make selectable
      updateNumericInput(session, "BootN_i_cdda", value = 2)

      shinyjs::show("BootType_i_cdda") # Make selectable
      updateSelectInput(session, "BootType_i_cdda", selected = "perc")

      shinyjs::show("CI_i_cdda") # Make selectable
      updateNumericInput(session, "CI_i_cdda", value = "0.95")

      shinyjs::show("BootHelpText_cdda")
      shinyjs::show("boot_section_i_cdda")
    } else {
      shinyjs::hide("BootN_i_cdda") # Make unclickable
      updateNumericInput(session, "BootN_i_cdda", value = NA)

      shinyjs::hide("BootType_i_cdda") # Make unclickable
      updateSelectInput(session, "BootType_i_cdda", selected = NULL)

      shinyjs::hide("CI_i_cdda") # Make unclickable
      updateNumericInput(session, "CI_i_cdda", value = NA)

      shinyjs::hide("BootHelpText_cdda")
      shinyjs::hide("boot_section_i_cdda")
    }
  })

  observeEvent(input$DiffTest_i_cdda, {
    if (input$DiffTest_i_cdda == "TRUE") {
      shinyjs::show("Paral_i_cdda") # Make selectable
    } else {
      shinyjs::hide("Paral_i_cdda") # Make unclickable
      updateSelectInput(session, "Paral_i_cdda", selected = "FALSE")

      shinyjs::hide("CoresN_i_cdda") # Make unclickable
      updateNumericInput(session, "CoresN_i_cdda", value = 1)
    }
  })

  observeEvent(input$Paral_i_cdda, {
    if (input$Paral_i_cdda == "TRUE") {
      shinyjs::show("CoresN_i_cdda") # Make selectable
    } else {
      shinyjs::hide("CoresN_i_cdda") # Make unclickable
      updateNumericInput(session, "CoresN_i_cdda", value = 1)
    }
  })

  # --- ObserveEvent: Show/Hide Plots
  observeEvent(input$DiffTest_i_cdda, {
    if (input$DiffTest_i_cdda == "TRUE") {
      shinyjs::show("plotControls2") # Show plots
    } else {
      shinyjs::hide("plotControls2") # Hide plots
    }
  })

  # --- renderUI: Show Plots and Plot Conditions for CDDA - Variables
  output$plotControls1 <- renderUI({
    req(runcdda_var())
    fluidRow(
      br(),
      h2("Plot Options"),
      column(
        3,
        selectInput(
          inputId = "stat_var",
          label = "Statistics",
          choices = c("coskew", "cokurt", "rhs", "rcc", "rtanh"),
          selected = "rhs"
        )
      ),
      column(
        6,
        splitLayout(
          cellWidths = c("50%", "50%"),
          numericInput("ylim_min_var", "Y-axis min", value = NA, step = 1),
          numericInput("ylim_max_var", "Y-axis max", value = NA, step = 1)
        )
      )
    )
  })

  # --- renderUI: Show Summary and Summary Conditions for CDDA - Variables
  output$summaryControls1 <- renderUI({
    req(runcdda_var())
    fluidRow(
      br(),
      h2("Summary Options"),
      column(
        3,
        selectInput("skew", "Skewness", choices = c("TRUE", "FALSE"), selected = "TRUE")
      ),
      column(
        3,
        selectInput("coskew", "Co-Skewness", choices = c("TRUE", "FALSE"), selected = "FALSE")
      ),
      column(
        3,
        selectInput("kurt", "Kurtosis", choices = c("TRUE", "FALSE"), selected = "TRUE")
      ),
      column(
        3,
        selectInput("cokurt", "Co-Kurtosis", choices = c("TRUE", "FALSE"), selected = "FALSE")
      )
    )
  })

  # --- renderUI: Show Plots and Plot Conditions for CDDA - Independence
  output$plotControls2 <- renderUI({
    req(runcdda_ind())
    fluidRow(
      br(),
      h2("Plot Options"),
      column(
        3,
        selectInput(
          inputId = "stat_ind",
          label = "Statistics",
          choices = c("hsic.diff", "dcor.diff", "mi.diff"),
          selected = "hsic.diff"
        )
      ),
      column(
        6,
        splitLayout(
          cellWidths = c("50%", "50%"),
          numericInput("ylim_min_ind", "Y-axis min", value = NA, step = 1),
          numericInput("ylim_max_ind", "Y-axis max", value = NA, step = 1)
        )
      )
    )
  })

  # --- renderUI: Show Summary and Summary Conditions for CDDA - Independence
  output$summaryControls2 <- renderUI({
    req(runcdda_ind())
    fluidRow(
      br(),
      h2("Summary Options"),
      column(
        6,
        selectInput("hsic", "HSIC", choices = c("TRUE", "FALSE"), selected = "TRUE"),
        selectInput("dcor", "Distance Correlation", choices = c("TRUE", "FALSE"), selected = "TRUE"),
        selectInput("nlfun", "Non-linear Correlation", choices = c("TRUE", "FALSE"), selected = "FALSE")
      ),
      column(
        6,
        selectInput("hsic.diff", "HSIC Difference", choices = c("TRUE", "FALSE"), selected = "FALSE"),
        selectInput("dcor.diff", "Distance Correlation Difference", choices = c("TRUE", "FALSE"), selected = "FALSE"),
        selectInput("mi.diff", "Mutual Information Difference", choices = c("TRUE", "FALSE"), selected = "FALSE")
      )
    )
  })

  # --- Reactive: Run CDDA - Variables
  runcdda_var <- reactive({
    req(data(), input$cause, input$effect, input$mod_cdda, input$BootN_v_cdda, input$CI_v_cdda, input$Modval_v_cdda)

    # --- Prepare Formula ----
    # 1. Combine numeric covariates if not NULL
    if (!is.null(input$Covar_n_cdda) && length(input$Covar_n_cdda) > 0) {
      covar_num <- paste(input$Covar_n_cdda, collapse = " + ")
    } else {
      covar_num <- NULL
    }

    # 2. Combine categorical covariates if not NULL
    if (!is.null(input$Covar_c_cdda) && length(input$Covar_c_cdda) > 0) {
      covar_cat <- paste(input$Covar_c_cdda, collapse = " + ")
    } else {
      covar_cat <- NULL
    }

    # 3. Combine numeric and categorical covariates if not NULL
    covars <- paste(c(covar_num, covar_cat), collapse = " + ")

    covars <- if (is.null(covars)) {
      NULL
    } else {
      covars
    }

    # 4. Form fomula
    if (is.null(input$Covar_c_cdda) && is.null(input$Covar_n_cdda)) {
      formula_str <- paste(input$effect, "~", paste(input$cause, "*", input$mod_cdda))
    } else {
      formula_str <- paste(input$effect, "~", paste(paste(input$cause, "*", input$mod_cdda), "+", covars))
    }

    formula <- as.formula(formula_str)

    # --- Change Categorical Covariates to Factors ----
    if (!is.null(input$Covar_c_cdda) && length(input$Covar_c_cdda) > 0) {
      df <- data()

      for (v in input$Covar_c_cdda) {
        if (v %in% names(df)) {
          df[[v]] <- as.factor(df[[v]])
        }
      }
    } else {
      df <- data()
    }

    # --- CDDA Variables Condition ----
    m <- lm(formula, data = df)
    pred <- input$cause
    mod <- input$mod_cdda

    if (is.na(input$BootN_v_cdda) || input$BootN_v_cdda == 1 || is.na(input$CI_v_cdda) || is.null(input$BootType_v_cdda)) {
      return(NULL)
    }

    out_var <- cdda.vardist(
      formula = m,
      pred = pred,
      mod = mod,
      modval = input$Modval_v_cdda,
      B = input$BootN_v_cdda,
      boot.type = input$BootType_v_cdda,
      conf.level = input$CI_v_cdda,
      data = df
    )
    return(out_var)
  })

  # --- renderPrint: Show CDDA - Variables Output
  output$cdda_var <- renderPrint({
    req(runcdda_var())
    print(runcdda_var())
  })

  # --- rendePlot: Show CDDA - Variables Plots
  output$cdda_var_plot <- renderPlot({
    req(runcdda_var())
    ymin <- input$ylim_min_var
    ymax <- input$ylim_max_var
    ylim_vals <- if (!is.na(ymin) && !is.na(ymax)) c(ymin, ymax) else NULL

    args <- list(
      x = runcdda_var(),
      stat = input$stat_var
    )
    if (!is.null(ylim_vals)) args$ylim <- ylim_vals
    do.call(plot, args)
  })

  # --- rendePrint: Show CDDA - Variables Summary
  output$cdda_var_summary <- renderPrint({
    req(runcdda_var())
    res <- summary(
      runcdda_var(),
      skew = input$skew,
      coskew = input$coskew,
      kurt = input$kurt,
      cokurt = input$cokurt
    )
    print(res)
  })

  # --- Reactive: Run CDDA - Independence
  runcdda_ind <- reactive({
    req(
      data(), input$cause, input$effect, input$mod_cdda, input$Modval_i_cdda, input$NlFun_i_cdda
    )

    # --- Prepare Formula ----
    # 1. Combine numeric covariates if not NULL
    if (!is.null(input$Covar_n_cdda) && length(input$Covar_n_cdda) > 0) {
      covar_num <- paste(input$Covar_n_cdda, collapse = " + ")
    } else {
      covar_num <- NULL
    }

    # 2. Combine categorical covariates if not NULL
    if (!is.null(input$Covar_c_cdda) && length(input$Covar_c_cdda) > 0) {
      covar_cat <- paste(input$Covar_c_cdda, collapse = " + ")
    } else {
      covar_cat <- NULL
    }

    # 3. Combine numeric and categorical covariates if not NULL
    covars <- paste(c(covar_num, covar_cat), collapse = " + ")

    covars <- if (is.null(covars)) {
      NULL
    } else {
      covars
    }

    # 4. Form fomula
    if (is.null(input$Covar_c_cdda) && is.null(input$Covar_n_cdda)) {
      formula_str <- paste(input$effect, "~", paste(input$cause, "*", input$mod_cdda))
    } else {
      formula_str <- paste(input$effect, "~", paste(paste(input$cause, "*", input$mod_cdda), "+", covars))
    }

    formula <- as.formula(formula_str)

    # --- Change Categorical Covariates to Factors ----
    if (!is.null(input$Covar_c_cdda) && length(input$Covar_c_cdda) > 0) {
      df <- data()

      for (v in input$Covar_c_cdda) {
        if (v %in% names(df)) {
          df[[v]] <- as.factor(df[[v]])
        }
      }
    } else {
      df <- data()
    }

    # --- CDDA Independence Condition ----
    m <- lm(formula, data = df)
    pred <- input$cause
    mod <- input$mod_cdda

    args <- list(
      formula = m,
      pred = pred,
      mod = mod,
      modval = input$Modval_i_cdda,
      hsic.method = input$HsicMethod_i_cdda,
      nlfun = input$NlFun_i_cdda,
      hetero = input$Hetero_i_cdda,
      diff = input$DiffTest_i_cdda,
      data = df
    )

    args$parallelize <- ifelse(input$DiffTest_i_cdda == "TRUE", FALSE, TRUE)
    args$cores <- ifelse(input$Paral_i_cdda == "FALSE", 0, input$CoresN_i_cdda)

    if (input$HsicMethod_i_cdda %in% c("boot", "permutation") || input$DiffTest_i_cdda == "TRUE") {
      if (!is.na(input$BootN_i_cdda)) args$B <- input$BootN_i_cdda
      if (!is.null(input$BootType_i_cdda)) args$boot.type <- input$BootType_i_cdda
      if (!is.na(input$CI_i_cdda)) args$conf.level <- input$CI_i_cdda
    }

    if ((input$HsicMethod_i_cdda %in% c("boot", "permutation") || input$DiffTest_i_cdda == "TRUE") &&
      (is.na(input$BootN_i_cdda) || input$BootN_i_cdda == 1 || is.na(input$CI_i_cdda) || is.null(input$BootType_i_cdda))) {
      return(NULL)
    }

    out_ind <- do.call(cdda.indep, args)

    return(out_ind)
  })

  # --- renderPrint: Show CDDA - Independence Output
  output$cdda_ind <- renderPrint({
    req(runcdda_ind())
    print(runcdda_ind())
  })

  # --- rendePlot: Show CDDA - Independence Plots
  output$cdda_ind_plot <- renderPlot({
    req(runcdda_ind())

    ymin <- input$ylim_min_ind
    ymax <- input$ylim_max_var
    ylim_vals <- if (!is.na(ymin) && !is.na(ymax)) c(ymin, ymax) else NULL

    args <- list(
      x = runcdda_ind(),
      stat = input$stat_ind
    )
    if (!is.null(ylim_vals)) args$ylim <- ylim_vals
    do.call(plot, args)
  })

  # --- rendePrint: Show CDDA - Independence Summary
  output$cdda_ind_summary <- renderPrint({
    req(runcdda_ind())
    cdda_ind_res <- runcdda_ind()
    res <- summary(
      cdda_ind_res,
      hsic = input$hsic,
      dcor = input$dcor,
      nlfun = input$nlfun,
      hsic.diff = input$hsic.diff,
      dcor.diff = input$dcor.diff,
      mi.diff = input$mi.diff
    )
    print(res)
  })

  # --- downloadHandler: Download report
  output$downloadReport2 <- downloadHandler(
    filename = "CDDA Report.pdf",
    content = function(file) {
      src <- normalizePath("CDDA_Report.Rmd")

      owd <- setwd(tempdir()) # temporarily dir
      on.exit(setwd(owd))
      file.copy(src, "CDDA_Report.Rmd", overwrite = TRUE)

      req(runcdda_var())
      ymin_var <- input$ylim_min_var
      ymax_var <- input$ylim_max_var
      ylim_vals_var <- if (!is.na(ymin_var) && !is.na(ymax_var)) c(ymin_var, ymax_var) else NULL

      args1 <- list(
        x = runcdda_var(),
        stat = input$stat_var
      )
      if (!is.null(ylim_vals_var)) args1$ylim_var <- ylim_vals_var

      cdda_var_plot <- function() {
        do.call(plot, args1)
      }

      req(runcdda_ind())
      ymin_ind <- input$ylim_min_ind
      ymax_ind <- input$ylim_max_ind
      ylim_vals_ind <- if (!is.na(ymin_ind) && !is.na(ymax_ind)) c(ymin_ind, ymax_ind) else NULL

      args2 <- list(
        x = runcdda_ind(),
        stat = input$stat_ind
      )
      if (!is.null(ylim_vals_ind)) args2$ylim_ind <- ylim_vals_ind

      cdda_ind_plot <- function() {
        do.call(plot, args2)
      }

      cdda_var_args <- list(
        object = runcdda_var(),
        skew = input$skew,
        coskew = input$coskew,
        kurt = input$kurt,
        cokurt = input$cokurt
      )

      cdda_var_summary <- function() {
        do.call(summary, cdda_var_args)
      }

      cdda_ind_args <- list(
        object = runcdda_ind(),
        hsic = input$hsic,
        dcor = input$dcor,
        nlfun = input$nlfun,
        hsic.diff = input$hsic.diff,
        dcor.diff = input$dcor.diff,
        mi.diff = input$mi.diff
      )

      cdda_ind_summary <- function() {
        do.call(summary, cdda_ind_args)
      }

      out <- render("CDDA_Report.Rmd",
        output_format = pdf_document(),
        params = list(
          runcdda_var = runcdda_var(),
          cdda_var_plot = cdda_var_plot,
          cdda_var_summary = cdda_var_summary,
          runcdda_ind = runcdda_ind(),
          cdda_ind_plot = cdda_ind_plot,
          cdda_ind_summary = cdda_ind_summary
        ),
        envir = new.env(parent = globalenv())
      )

      file.rename(out, file)
    }
  )
}

# 4. --- A Call to The ShinyApp ----
shinyApp(ui = ui, server = server)


# Places could potentially improve:
# 1. Show more descriptive statistics
# 2. Diagnosis of DDA
# 3. Why is there a Null in the CDDA output?
# 4. Add more help texts
# 4. Make the report.Rmd neat!

# Problem 1
# Use bca as boottype:
# DDA Variables Error Message: estimated adjustment 'a' is NA
# DDA RESIDUAL Error message : argument "prob.trans" is missing, with no default
# DDA Independence Error Message: estimated adjustment 'a' is NA
# Do you think a help text or a warning message should put there?

# Problem 2
# When I select boot for HSIC Inference method
# Error: object 'critical_value' not found
