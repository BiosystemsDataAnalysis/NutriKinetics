# NutriKinetics program 
library(ggplot2)
library(shinydashboard)
library(minpack.lm)
library(DT)

ui <- fluidPage(
  titlePanel("Nutrikinetics"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(4, 
          selectInput("subjects", 
                      label = "Number of subjects", width = "80px",
                      choices = c("1", "2", "3", "4", "5",
                                  "6", "7", "8", "9"),
                      selected = 9
          )
        ),
        column(8,
          selectInput("times", multiple = TRUE,
                      label = HTML("Measurement times"),
                      choices  = sort(c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 
                                 14, 16, 20, 24, 28, 32)),
                      selected = sort(c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 
                                  14, 16, 20, 24, 28, 32)))
        )
      ),
      fluidRow(
        box(width = 6,
          radioButtons("type", label = "Simulation type",
                    choices = c("Urine" = "urine",
                                "Plasma" = "plasma"),
                     selected = "urine"
                     )
        ),
        box(width = 6,
          radioButtons("error.type", label = "Error type",
                     choices = c("Constant" = "constant",
                                 "Linear" = "linear"),
                     selected = "constant")
          )
        ),
      
      conditionalPanel(condition = paste0("input['", "type", "'] == 'urine' "),
        div(style = "font-size: 12pt;",
          sliderInput("xmax", 
                      label = HTML("x<sub>max</sub> (&micro;M)"), 
                      value = 10, min = 0.1, max = 20, step = 0.1), 
          sliderInput("ke",
                      label = HTML("k<sub>e</sub> (h<sup>-1</sup>)"),
                      value = 0.1, min = 0.01, max = 0.2, step = 0.001),
          sliderInput("tau.u",
                      label = HTML("&tau; (h)"),
                      value = 6, min = 1, max = 10, step = 0.1)
        )
      ),
      conditionalPanel(condition = paste0("input['", "type", "'] == 'plasma' "),
                       div(style = "font-size: 12pt;",
                           sliderInput("Cd", 
                                       label = HTML("Cd (nM)"), 
                                       value = 40, min = 30, max = 50, step = 0.1), 
                           sliderInput("kel",
                                       label = HTML("k<sub>el</sub> (h<sup>-1</sup>)"),
                                       value = 0.1, min = 0.1, max = 2, step = 0.01),
                           sliderInput("ka",
                                       label = HTML("k<sub>a</sub> (h<sup>-1</sup>)"),
                                       value = 1, min = 0.1, max = 5, step = 0.1),
                           sliderInput("tau.p",
                                       label = HTML("&tau; (h)"),
                                       value = 0.4, min = .1, max = 1, step = 0.01)
                       )
      ),
      conditionalPanel(condition = paste0("input['", "error.type", "'] == 'linear' "),
        sliderInput("errorl",
                    label = HTML("&epsilon; % of calculated value"),
                    value = 0, min = 0, max = 50, step = 0.1)
        ),
      conditionalPanel(condition = paste0("input['", "error.type", "'] == 'constant' "),
                       sliderInput("errorc",
                                   label = HTML("&epsilon;"),
                                   value = 0, min = 0, max = 5, step = 0.1)
      )
  ),
    
    mainPanel(
      fluidRow(
        column(6, 
               div(style = "font-weight: bold;",
          textOutput("text1")),
          plotOutput("plot1")),
        column(6, 
          dataTableOutput("table"))
      ),
      div(style="padding: 20px;"),
      fluidRow(
        column(2, actionButton("update.all", "Update all subjects")
        ),
        column(1),
        column(2, actionButton("update.single", "Update this subject")
#        ),
#        column(1),
#        column(2, actionButton("update.error", "Update errors")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
## Functions

###########################
  
  urine <- function(times, params) {
    
    xmax           <- params[[1]]
    ke             <- params[[2]]
    tau            <- params[[3]]
    
    y              <- xmax*(1 - exp(-ke*(times - tau))) 
    y[times < tau] <- 0
    
    return(y)

  }
  
############################
  
  plasma <- function(times, params) {
    
    Cd             <- params[[1]]
    ka             <- params[[2]]
    kel            <- params[[3]]
    tau            <- params[[4]]
  
  # Make sure ka and kel are not equal
    
    if(abs(ka - kel) < 0.01) {
      ka <- kel + 0.01
    }
    
    y              <- Cd*(ka/(ka - kel))*(exp(-kel*(times - tau)) - 
                      exp(-ka*(times - tau)))
    y[times < tau] <- 0

    return(y)
  }

##############################
  # Difference between calculated and "measured" values
  
  urine.res <- function(parms, times, y) {
    res <- urine(times, parms) - y
    return(res)
  } 
  
############################## 
  # Difference between calculated and "measured" values
  
  plasma.res <- function(parms, times, y) {
    res <- plasma(times, parms) - y
    return(res)
  } 
  
  
##############################
  
  doFit <- function(pars, fn, times, y) {
    fn <- eval(parse(text = fn))
    fit.lm <-  tryCatch({
      fit   <- nls.lm(par = pars, lower = NULL, 
                       upper = NULL,
                       fn = fn, times = times,
                       control = nls.lm.control(maxiter = 100),
                       y = y)
    },
    error = function(e) {
      print("Fit error")
    }, 
    warning = function(e) {
      print("Fit does not converge")
    }
    )
    return(fit.lm)
  }
    
##############################
  
  plotData <- function(df) {
    p <- ggplot(data = df)
    p <- p + geom_point(aes(x = times, 
                            y = y, 
                            col = type), 
                        size = 2)
    p <- p + geom_line(data = df[df$type == "fit", ], aes(x = times, y = y))
    if(input$type == "plasma") {
      p <- p + labs(x = "Time (h)", y = "Concentration (nM)")
    }
    else if(input$type == "urine") {
      p <- p + labs(x = "Time (h)", y = "Excretion (\u03bc M)")
    }
    p <- p + theme(legend.position = c(0.9, 0.2))
    p <- p + theme(legend.title = element_blank(),
                   legend.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   axis.text = element_text(size = 14))
    return(p)
  }
  
#############################  
  
  doSubject <- function() {
    
    error.type   <- input$error.type
    times        <- as.numeric(input$times)
    
    if(input$type == "urine") {
      xmax       <- as.numeric(input$xmax)
      ke         <- as.numeric(input$ke)
      tau        <- as.numeric(input$tau.u)
      
      pars       <- list(xmax = xmax, ke = ke, tau = tau)
      exp.data   <- urine(times, pars)
      
    # Add noise
      if(error.type == "linear") {
        error      <- (as.numeric(input$errorl) * exp.data / 100.) * 
                       rnorm(length(exp.data), 0, 1)
      } else if(error.type == "constant" ) {
        error      <- as.numeric(input$errorc) * rnorm(length(exp.data), 0, 1)
      }
    # Make sure all data points are positive
      exp.data   <- abs(exp.data + error)
      
      fit.out    <- doFit(pars, "urine.res", times, exp.data) 
      fit.data   <- urine(times, fit.out$par)
      
    } else if (input$type == "plasma") {
        Cd       <- as.numeric(input$Cd)
        ka       <- as.numeric(input$ka)
        kel      <- as.numeric(input$kel)
        tau      <- as.numeric(input$tau.p) 
        
        pars     <- list(Cd = Cd, ka = ka, kel = kel, tau = tau)
        exp.data <- plasma(times, pars)
        
      # Add noise
        if(error.type == "linear") {
          error      <- (as.numeric(input$errorl) * exp.data / 100.) * 
            rnorm(length(exp.data), 0, 1)
        } else if(error.type == "constant" ) {
          error      <- as.numeric(input$errorc) * rnorm(length(exp.data), 0, 1)
        }
      # Make sure all data points are positive
        exp.data   <- abs(exp.data + error)
        
        fit.out  <- doFit(pars, "plasma.res", times, exp.data)
        fit.data <- plasma(times, fit.out$par)
    }
    plot.data    <- c(exp.data, fit.data)

    params.out   <- unlist(fit.out$par)

    return(list(plot.data, params.out))
  }
   
  updateFunction <- function()  {   
    n.subjects                <- as.numeric(input$subjects)
    times                     <- as.numeric(input$times)
    df.data$plots             <- as.data.frame(matrix(nrow = 2*length(times),
                                                      ncol = n.subjects + 2))
    colnames(df.data$plots)   <- c("times", "type", paste0(seq(n.subjects)))
    df.data$plot$times        <- c(times, times)
    df.data$plot$type         <- rep(c("data", "fit"), each = length(times)) 
    
    plots                     <- lapply(seq(n.subjects), 
                                        function(x) {doSubject()[[1]]
                                        })
    df.data$plots             <- do.call(cbind, plots)
    
    df.data$plots             <- cbind.data.frame(c(times, times), 
                                                  rep(c("data", "fit"), each = length(times)),
                                                  df.data$plots)
    colnames(df.data$plots)   <- c("times", "type", paste0(seq(n.subjects)))
    
    # Prepare data frame with parameters, df.data$table
    
    table.out                 <- lapply(seq(n.subjects), 
                                        function(x) {doSubject()[[2]]
                                        })
    df.data$table           <- do.call(rbind.data.frame, table.out) 
    
    # add means and std
    df.data$table             <- rbind(df.data$table, colMeans(df.data$table),
                                       sapply(df.data$table, sd))
    tmp                       <- data.frame(seq(n.subjects))
    tmp                       <- rbind(tmp, "Means", "Std. dev")
    df.data$table             <- cbind.data.frame(tmp, df.data$table)
    if(input$type == "urine") {
      colnames(df.data$table) <- c("Subject", "xmax", "ke", "tau")
    } else if(input$type == "plasma") {
      colnames(df.data$table) <- c("Subject", "Cd", "ka", "kel", "tau")
    }
    dd$select                 <- 1
  }
  
#################################  
  

## Reactives  
  
    df.data                         <- reactiveValues()
    dd                              <- reactiveValues(select = NULL)

    observe({
      if(input$update.all == 0) {   # run at start-up
         isolate({
           updateFunction()
         })
      }
    })  
    
    observeEvent(input$table_rows_selected ,{
      dd$select   <- isolate(input$table_rows_selected)
    })
  
    observeEvent(input$type, {
      updateFunction()
    })
    
    
    observeEvent(input$subjects, {
      updateFunction()
    })
    
    observeEvent(input$times, {
      updateFunction()
      updateSelectInput(session, "times", 
                        selected = sort(as.numeric(input$times)))
    })
    
    observeEvent(input$update.all, {
      updateFunction()
    })
 
    
  ## Update single subject
    
    observeEvent(input$update.single, {
 
      n.subjects              <- as.numeric(input$subjects)
      
      isolate({
      selected.subject        <- dd$select
      if(is.null(selected.subject)) {
        selected.subject      <- 1
      }
      if(selected.subject > n.subjects) {
        selected.subject      <- 1
      }
  
      subject.list                          <- doSubject()
    
      df.data$plots[, 2 + selected.subject] <- subject.list[[1]]
      df.data$table[selected.subject, -1]   <- subject.list[[2]]
      
   # update subjects, means and std
      df.data$table[n.subjects + 1, -1] <- colMeans(df.data$table[1:n.subjects, -1 ], 
                                                    na.rm = TRUE)
      df.data$table[n.subjects + 2, -1] <- sapply(df.data$table[1:n.subjects, - 1],
                                                  function(x) sd(x, na.rm = TRUE))
      })
    
 })    
    output$text1       <- renderText({
      n.subjects       <- as.numeric(input$subjects)
      displayed.plot   <- dd$select
  
      if(is.null(displayed.plot) ) {
        displayed.plot <- 1
      }
      if(displayed.plot > n.subjects) {
        displayed.plot <- 1
        dd$select      <- 1
      }
      plot.text <- paste0("Subject ", displayed.plot)
      plot.text
    })
    
    output$plot1       <- renderPlot({
      n.subjects       <- as.numeric(input$subjects)
      displayed.plot    <- dd$select
  
      if(is.null(displayed.plot)) {
        displayed.plot <- 1
      }
      if(displayed.plot > n.subjects) {
        displayed.plot <- 1
      }
      req((2 + n.subjects) <= ncol(df.data$plots))
      df               <- df.data$plots[, c(1, 2, 2 + displayed.plot)]
      
      colnames(df)     <- c("times", "type", "y")
      p                <- plotData(df)
      p
    })
    
    
    output$table <- DT::renderDataTable({
      par.table           <- df.data$table

      colnames(par.table) <- colnames(df.data$table)
      par.table[, 2] = round(x = par.table[,2],digits = 2)
      par.table[, 3] = round(x = par.table[,3],digits = 3)
      par.table[, 4] = round(x = par.table[,4],digits = 2)
      if(input$type == "plasma") {
        par.table[, 5] = round(x = par.table[, 5],digits = 3)
      }
      DT::datatable(par.table,
        selection = list(mode = "single", target = "row", 
                         selected = as.numeric(dd$select)), 
        rownames = FALSE,
        options = list(
          paging = FALSE,
          searching = FALSE,
          dom = 't',
          ordering = FALSE
      )
    ) 
  })
}

shinyApp(ui, server)