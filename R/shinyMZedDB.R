#' shinyMZedDB
#' @description A shiny application for the putative annotation of high resolution m/z.
#' @export
#' @author Jasen Finch
#' @importFrom shiny shinyApp navbarPage tabPanel fluidRow column numericInput checkboxInput reactive
#' @importFrom shiny textInput plotOutput sidebarPanel mainPanel tableOutput renderPlot renderTable
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom rcdk parse.smiles view.image.2d
#' @importFrom graphics plot par rasterImage
#' @examples
#' \dontrun{
#' shinyMZedDB()
#' }

shinyMZedDB <- function(){
  shinyApp(
    ui = navbarPage("Shiny-MZedDB",
                    tabPanel("Putative Ionisation Product",
                             # Create a new Row in the UI for selectInputs
                             fluidRow(
                               column(3, numericInput("acc_mz","Accurate Mass:",117.078979)),
                               column(2, numericInput("ppm","PPM:",1)),
                               column(2, checkboxInput("mode_p", "Positive Mode",value = F)),
                               column(2, checkboxInput("mode_ne", "Neutral",value = T)),
                               column(2, checkboxInput("mode_n", "Negative Mode",value = F)),
                               column(2, checkboxInput("pip_iso", "Isotopes",value = F))
                               #column(4, selectInput('add_sel_i', 'Adducts:', adducts, multiple=TRUE, selectize=FALSE))                                                    
                             ),
                             fluidRow(
                               dataTableOutput(outputId = "pipTable")
                             ),
                             fluidRow(
                               plotOutput('Structure')
                             )
                    ),
                    tabPanel("Molecular Formula Generator",
                             fluidRow(
                               column(3, numericInput("acc_mass","Accurate Mass:",117.078979)),
                               column(2, numericInput("accuracy","PPM:",1)),
                               column(2, numericInput("charge", "Charge:",0)),
                               column(2, checkboxInput("applygr", "Apply 7GR",value = T))
                             ),
                             # new row for elemental composition
                             fluidRow(
                               column(2, numericInput("Cmax","C:",10)),
                               column(2, numericInput("Hmax","H:",20)),
                               column(2, numericInput("Nmax","N:",5)),
                               column(2, numericInput("Omax","O:",5)),
                               column(2, numericInput("Pmax","P:",0)),
                               column(2, numericInput("Smax","S:",0)),
                               column(2, numericInput("Namax","Na:",0)),
                               column(2, numericInput("Kmax","K:",0)),
                               column(2, numericInput("Clmax","Cl:",0)),
                               column(2, numericInput("iCmax","C13:",0)),
                               column(2, numericInput("iOmax","O18:",0)),
                               column(2, numericInput("iKmax","K41:",0)),
                               column(2, numericInput("iClmax","Cl37:",0))
                             ),
                             # Create a new row for the table.
                             fluidRow(
                               dataTableOutput(outputId = "mf_table")
                             ),
                             fluidRow(
                               plotOutput('isoDistsPlot')
                             ),
                             fluidRow(
                               dataTableOutput('isoDistsTable')
                             )
                    ),
                    tabPanel("Adduct Calculator",
                             fluidRow(
                               dataTableOutput(outputId = "MZDB")
                             ),
                             fluidRow(
                               tableOutput("adduct")
                             ) 
                    )
    ),
    server = function(input, output) {
      
      getPIPs <- reactive({
        if (input$mode_p) {
          mode <- "p"
        }
        if (input$mode_n) {
          mode <- "n"
        }
        if (input$mode_ne) {
          mode <- "ne"
        }
        if (input$pip_iso) {
          iso <- c('C13','O18','S34')
        } else {
          iso <- NULL
        }
        pip_tab <- PIPsearch(input$acc_mz,mode,input$ppm,iso = iso)
        pip_tab
      })
      
      output$pipTable <- renderDataTable({
        pip_tab <- getPIPs()
        pip_tab
      },server = T,selection = 'single',rownames = F) 
      
      output$Structure <- renderPlot({
        r <- input$pipTable_rows_selected
        if (!is.null(r)) {
          res <- getPIPs()
          res <- res[r,]
          smile <- gsub('"','',res$Smiles)
          par(mar = c(0,0,0,0))
          sm <- parse.smiles(smile)[[1]]
          temp1 <- view.image.2d(sm,500,500)
          plot(NA,NA,xlim = c(1,100),ylim = c(1,100),xaxt = 'n',yaxt = 'n',xlab = '',ylab = '')
          rasterImage(temp1,1,1,100,100)
        } else {
          NULL
        }
      },width = 400,height = 400)
      
      getMF <- reactive({
        mz <- as.numeric(input$acc_mass)
        maxi <- c(C = input$Cmax,iC = input$iCmax,H = input$Hmax,iH = 0,N = input$Nmax,iN = 0,O = input$Omax,iO = input$iOmax,F = 0 ,Na = input$Namax,Si = 0,P = input$Pmax,S = input$Smax,Cl = input$Clmax,iCl = input$iClmax,Br = 0,iBr = 0,K = input$Kmax,iK = input$iKmax)
        
        res <- generateMF(mz,ppm = input$accuracy,charge = input$charge,applygr = input$applygr,composition = maxi)
        
        res
      })
      
      output$mf_table <- renderDataTable({
        res <- getMF()
        res
      },server = T,selection = 'single',rownames = F) 
      
      output$isoDistsPlot <- renderPlot({
        r <- input$mf_table_rows_selected
        if (!is.null(r)) {
          mf <- getMF()
          mf <- mf[r,]
          if (grepl('i',mf)) {
            NULL
          } else {
            res <- isoDistr(as.character(mf$MF),chrg = input$charge)
            plot(res[,1],res[,3],type = "h",xlab = "m/z",ylab = "Relative Intensity",col = "Blue")
          }
        } else {
          NULL
        }
      })
      
      output$isoDistsTable <- renderDataTable({
        r <- input$mf_table_rows_selected
        if (!is.null(r)) {
          mf <- getMF()
          mf <- mf[r,]
          if (grepl('i',mf)) {
            NULL
          } else {
            res <- isoDistr(as.character(mf$MF),chrg = input$charge)
            res <- res[,-2]
            colnames(res)[2] <- "Relative Intensity"
            res
          }
        } else {
          NULL
        }
      })
      
      output$MZDB <- renderDataTable({
        res <- MZedDB$MZedDB_ALL
        res <- res[,-c(3,6,7,8,10)]
        names(res)[4:5] <- c('Accurate Mass', 'SMILE')
        res
      },server = T,selection = 'single',rownames = F)
      
      output$adduct <- renderTable({
        r <- input$MZDB_rows_selected
        if (!is.null(r)) {
        mz <- MZedDB$MZedDB_ALL$Accurate.Mass[r]
        metrules <- MZedDB$MZedDB_METRULES[r,]
        print(class(metrules))
        print(metrules)
        rules <- MZedDB$ADDUCT_FORMATION_RULES
        res <- apply(rules,1,function(rule,id,metrules){
          Nch <- metrules$Nch
          Nacc <- metrules$Nacc
          Ndon <- metrules$Ndon
          Nnhh <- metrules$Nnhh
          Noh <-  metrules$Noh
          Ncooh <- metrules$Ncooh
          Ncoo <-  metrules$Ncoo
          if (eval(parse(text = as.character(rule[6])))) {
            if (as.numeric(rule[2]) != 0) {
              mz <- mz/as.numeric(rule[2])
            }
            mz <- mz + as.numeric(rule[4])
            mz <- round(mz,5)
          } else {
            mz <- 'Not Possible'
          }
          return(mz)
        },id = mz,metrules = metrules)
        res <- cbind(rules$Name,res)
        colnames(res) <- c('Adduct','m/z')
        res
        } else {
          NULL
        }
      })
    }
  )
}