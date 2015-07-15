viewAnnotation <- function(){
  options(shiny.maxRequestSize = 20*1024^2)
  options(digits=10)
  shinyApp(
    ui = fluidPage(
      titlePanel("viewAnnotation"),
      sidebarLayout(
        sidebarPanel(
          fileInput('file1', 'Choose file to upload',
                    accept = c(
                      '.RData'
                    )
          ),
          tags$hr(),
          radioButtons('mode', 'Mode',
                       c(Positive='Positive_Mode',
                         Negative='Negative_Mode')),
          uiOutput("mz")
        ),
        mainPanel( 
          tabsetPanel(
            tabPanel("Accurate m/z", dataTableOutput('accurate')),
            tabPanel("Correlation Analysis", dataTableOutput('correlation')),
            tabPanel("Molecular Formulas",dataTableOutput('MF')),
            tabPanel("Theoretical Isotope Distributions",
                     fluidRow(uiOutput("selectMF")),
                     fluidRow(tableOutput('IsoDist')),
                     fluidRow(plotOutput('isoplot'))
            ),
            tabPanel("Putative Ionisation Products",dataTableOutput('PIP'))
          )
        )
      )
    ),
    
    server = function(input, output) {
      loadData <- reactive({
        if(is.null(input$file1)) return(NULL)
        load(input$file1$datapath)
        annot.res
      })
      availmz <- reactive({
        if(is.null(loadData())) return(NULL)
        annot.res <- loadData()
        annot.res[[input$mode]][["Accurate m/z"]][,1]
      })
      availMF <- reactive({
        if(is.null(loadData())) return(NULL)
        annot.res <- loadData()
        names(annot.res[[input$mode]][["Theoretical Isotope Distributions"]][[input$selectedmz]])
      })
      output$accurate <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Accurate m/z"]]
          rownames(res) <- NULL
          return(res)
        }
      })
      output$mz <- renderUI({
        selectInput('selectedmz', 'm/z', availmz())
      })
      output$selectMF <- renderUI({
        selectInput('selectedMF', 'MF', availMF())
      })
      output$correlation <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Correlation Analysis"]][[input$selectedmz]]
          names(res)[c(1,4,5)] <- c("Bin","Relation To","Relation From")
          if(nrow(res)>0){
            res[,4:5] <- apply(res[,4:5],2,function(x){x[which(x=="NA")] <- "";return(x)})
          }
          return(res)
        }
      })
      output$MF <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Molecular Formulas"]][[input$selectedmz]]
          return(res)
        }
      })
      output$IsoDist <- renderTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Theoretical Isotope Distributions"]][[input$selectedmz]][[input$selectedMF]]
          if(!is.null(res)){
            colnames(res)[3] <- "Relative Intensity"
          }
          return(res)
        }
      })
      output$isoplot <- renderPlot({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Theoretical Isotope Distributions"]][[input$selectedmz]][[input$selectedMF]]
          if(!is.null(res)){
            plot(res[,1],res[,3],type="h",col="Blue",xlab="m/z",ylab="Relative Intensity")
          }
        }
      })
      output$PIP <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          annot.res <- loadData()
          res <- annot.res[[input$mode]][["Putative Ionisation Products"]][[input$selectedmz]]
          return(res)
        }
      })
    }
  )
}