

# Source the R file containing the function definitions
source("POD-TPI_prob.R")
source("TITE-setup.R")


# Define UI
ui <- fluidPage(
  
  titlePanel(
    tags$div(
      style = "text-align: center;",
      tags$h1(style = "color: purple;", "POD_TPI Decision"),
      tags$p(style = "font-size: 16px;", "Last Updated: 07/26/2024"),
      tags$p(style = "font-size: 16px; font-weight: bold;", "Author1, Author2, Author3"),
      tags$p(style = "font-size: 16px;", "Department of Population Health Sciences, Weill Cornell Medicine")
    )
  ),
  
  tabsetPanel(
    
    tabPanel("Description", includeHTML('./Description.html')),
    
  tabPanel("Main",
  sidebarLayout(
    sidebarPanel(
      numericInput("current_dose", "Current Dose:", value = 2),
      numericInput("W", "DLT assessment window:", value = 28),
      numericInput("MaxDose", "Max Dose:", value = 4),
      textInput("dose_vec", "Dose Vector (comma-separated):", value = "2, 2, 2, 2, 2, 2"),
      textInput("event_time_vec", "Event Time Vector (comma-separated):", value = "28, 28, 9, 26, 15, 8"),
      textInput("event_vec", "Event Vector (comma-separated):", value = "0, 0, 1, 1, NA, NA"),
      numericInput("n_d", "DLT Patients:", value = 2),
      numericInput("m_d", "non-DLT Patients:", value = 2),
      numericInput("r_d", "Pending Patients:", value = 2),
      numericInput("p_T", "Target DLT Probability:", value = 0.3),
      numericInput("epsilon1", "Left Equivalence Interval Bound of Target DLT Probability:", value = 0.05),
      numericInput("epsilon2", "Right Equivalence Interval Bound of Target DLT Probability:", value = 0.05),
      numericInput("niter", "Number of Iterations for MCMC:", value = 1000),
      selectInput("type_dose_decision", "Type Dose Decision:", choices = c("mTPI", "i3+3"), selected = "mTPI"),
      selectInput("type_p_prior", "Type P Prior:", choices = c("uniform", "semiparametric"), selected = "uniform"),
      selectInput("type_t_model", "Type T Model:", choices = c("uniform", "pwuniform", "dhazard", "pwhazard"), selected = "pwuniform"),
      numericInput("type_suspension", "Type of Suspension:", value = 2),
      numericInput("q1", "pi_E for Trial Suspension:", value = 1),
      numericInput("q2", "pi_D for Trial Suspension:", value = 0.15),
      actionButton("compute", "Compute")
    ),
    mainPanel(
      verbatimTextOutput("result")
    )
  )
  )
  )
)

# Define server logic
server <- function(input, output) {
  observeEvent(input$compute, {
    # Parse input vectors
    dose_vec <- as.numeric(unlist(strsplit(input$dose_vec, ",")))
    event_time_vec <- as.numeric(unlist(strsplit(input$event_time_vec, ",")))
    event_vec <- as.numeric(unlist(strsplit(input$event_vec, ",")))
    
    # Call the POD_TPI_prob function
    result <- POD_TPI_prob(
      current_dose = input$current_dose,
      W = input$W,
      MaxDose = input$MaxDose,
      dose_vec = dose_vec,
      event_time_vec = event_time_vec,
      event_vec = event_vec,
      n_d = input$n_d,
      m_d = input$m_d,
      r_d = input$r_d,
      p_T = input$p_T,
      epsilon = c(input$epsilon1, input$epsilon2),
      niter = input$niter,
      type_dose_decision = input$type_dose_decision,
      type_p_prior = input$type_p_prior,
      type_t_model = input$type_t_model,
      type_suspension = input$type_suspension,
      q1 = input$q1,
      q2 = input$q2
    )
    
    # Display the result
    output$result <- renderText(result)
  })
}

shinyApp(ui = ui, server = server)