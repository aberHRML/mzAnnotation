#' viewAnnotationLC
#' @export

viewAnnotationLC <- function(){
  options(shiny.maxRequestSize = 20*1024^2)
  options(digits=10)
  shinyApp(
    ui = fluidPage(
      titlePanel("viewAnnotationLC"),
      sidebarLayout(
        sidebarPanel(
          fileInput('file1', 'Upload Annotation:',
                    accept = c(
                      '.RData'
                    )
          ),
          tags$hr(),
          uiOutput("mz"),
          tags$hr(),
          textInput("boxplots","Box Plot Folder Path:")
        ),
        mainPanel( 
          tabsetPanel(
            tabPanel("Peak Table", dataTableOutput('peakTab')),
            tabPanel("Box Plots", imageOutput('boxplot')),
            tabPanel('Correlation Analysis',dataTableOutput('corAnalysis')),
            tabPanel("Putative Ionisation Products",
                     fluidRow(dataTableOutput('PIP'))
            ),
            tabPanel("PIP Structures",
                     fluidRow(uiOutput("selectPIP")),
                     plotOutput('smile')
            )
          )
        )
      )
    ),
    
    server = function(input, output) {
      loadData <- reactive({
        if(is.null(input$file1)) return(NULL)
        load(input$file1$datapath)
        featAnnot
      })
      availmz <- reactive({
        if(is.null(loadData())) return(NULL)
        featAnnot <- loadData()
        names(featAnnot$featPIP)
      })
      availPIP <- reactive({
        if(is.null(loadData())) return(NULL)
        featAnnot <- loadData()
        featAnnot <- featAnnot$featPIP[[input$selectedmz]][,2]
        featAnnot <- gsub('"','',featAnnot)
        featAnnot
      })
      output$mz <- renderUI({
        selectInput('selectedmz', 'm/z', availmz())
      })
      output$selectPIP <- renderUI({
        selectInput('selectedPIP', 'Putative Ionisation Product', availPIP())
      })
      output$peakTab <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          featAnnot <- loadData()
          res <- featAnnot[["peakTable"]]
          rownames(res) <- NULL
          res <- data.frame(res)
          res[,c('rt','rtmin','rtmax')] <- round(res[,c('rt','rtmin','rtmax')]/60,3)
          return(res)
        }
      })
      output$boxplot <- renderImage({
        files <- list.files(input$boxplots)
        file <- grep(input$selectedmz,files)
        if(length(file)>0){
          path <- paste(input$boxplots,files[file],sep="/")
        } else {
          path <- ""
        }
        list(src=path)
      },deleteFile = FALSE)
      output$corAnalysis <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          featAnnot <- loadData()
          res <- featAnnot$corAnalysis[[input$selectedmz]]
          colnames(res)[1] <- 'ID'
          colnames(res)[6:7] <- c('Relation To','Relation From')
          res$r <- round(res$r,3)
          if(nrow(res)>0){
            res[,6:7] <- apply(res[,6:7],2,function(x){x[which(x=="NA")] <- "";return(x)})
          }
          return(res)
        }
      })
      output$PIP <- renderDataTable({
        if (is.null(loadData())){
          return(NULL)
        } else {
          featAnnot <- loadData()
          res <- featAnnot$featPIP[[input$selectedmz]]
          names(res)[c(4,7,8)] <- c("Accurate Mass", "m/z","PPM Error")
          res[,4] <- round(res[,4],5)
          res[,7] <- round(res[,7],5)
          res[,8] <- round(res[,8],3)
          return(res)
        }
      })
      output$smile <- renderPlot({
        if (is.null(loadData())){
          return(NULL)
        } else {
          featAnnot <- loadData()
          res <- featAnnot$featPIP[[input$selectedmz]]
          if(nrow(res)>0){
            res[,2] <- gsub('"','',res[,2])
            smile <- res[which(res[,2]==input$selectedPIP),5]
            par(mar=c(0,0,0,0))
            sm <- parse.smiles(smile)[[1]]
            temp1 <- view.image.2d(sm,500,500)
            plot(NA,NA,xlim=c(1,100),ylim=c(1,100),xaxt='n',yaxt='n',xlab='',ylab='') 
            rasterImage(temp1,1,1,100,100) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
          }
        }},width=400,height=400)
    }
  )
}