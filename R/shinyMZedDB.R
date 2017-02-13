#' shinyMZedDB
#' @description A shiny application for the putative annotation of high resolution m/z.
#' @export
#' @author Jasen Finch
#' @importFrom shiny shinyApp navbarPage tabPanel fluidRow column numericInput checkboxInput reactive
#' @importFrom shiny textInput plotOutput sidebarPanel mainPanel tableOutput renderPlot renderTable
#' @importFrom shiny selectInput uiOutput renderUI
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
                               column(3, numericInput("acc_mz","Accurate Mass:",118.08626)),
                               column(2, numericInput("ppm","PPM:",1)),
                               column(2, selectInput("modeI", "Mode:",c('positive','neutral','negative'))),
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
                               column(2, selectInput("Mode", "Mode:",c('positive','neutral','negative'))),
                               column(2, selectInput('Isotope', 'Isotope:',c('none','C13','2C13','O18','S34'))),
                               column(2, uiOutput('adducts')),
                               column(2, checkboxInput("applygr", "Apply 7GR",value = T))
                             ),
                             # new row for elemental composition
                             fluidRow(
                               column(2, numericInput("Cmax","C:",10)),
                               column(2, numericInput("Hmax","H:",20)),
                               column(2, numericInput("Nmax","N:",5)),
                               column(2, numericInput("Omax","O:",5)),
                               column(2, numericInput("Pmax","P:",0)),
                               column(2, numericInput("Smax","S:",0))
                             ),
                             # Create a new row for the table.
                             fluidRow(
                               dataTableOutput(outputId = "mf_table")
                             ),
                             # fluidRow(
                             #   plotOutput('isoDistsPlot')
                             # ),
                             # fluidRow(
                             #   dataTableOutput('isoDistsTable')
                             # ),
                             fluidRow(
                               dataTableOutput(outputId = 'mfHitsTable')
                             ),
                             fluidRow(
                               plotOutput('MFstructure')
                             )
                    ),
                    tabPanel("Adduct Calculator",
                             fluidRow(
                               dataTableOutput(outputId = "MZDB")
                             ),
                             fluidRow(
                               tableOutput("adduct")
                             ) 
                    ),
                    tabPanel('Relationship Prediction',
                             fluidRow(
                               column(3,numericInput('mz1','m/z 1:',135.0288)),
                               column(3,numericInput('mz2','m/z 2:',136.03216)),
                               column(2, selectInput("mode", "Mode:",c('positive','negative')))
                             ),
                             fluidRow(
                               dataTableOutput(outputId = 'predictedRelTable')
                             )
                             )
    ),
    server = function(input, output) {
      
      ###################### PIPs ##################################
      getPIPs <- reactive({
        if (input$modeI == 'positive') {
          mode <- "p"
        }
        if (input$modeI == 'negative') {
          mode <- "n"
        }
        if (input$modeI == 'neutral') {
          mode <- "ne"
        }
        if (input$pip_iso) {
          iso <- c('C13','2C13','O18','S34')
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
      
      ############################## MFs ############################################
      
      getAdducts <- reactive({
        mode <- input$Mode
        adducts <- MZedDB$ADDUCT_FORMATION_RULES
        if (mode == 'positive') {
          adducts <- adducts$Name[adducts$Nelec < 0] 
        }
        if (mode == 'neutral') {
          adducts <- adducts$Name[adducts$Nelec == 0] 
        }
        if (mode == 'negative') {
          adducts <- adducts$Name[adducts$Nelec > 0] 
        }
        adducts <- rev(adducts)
        return(adducts)
      })
      
      output$adducts <- renderUI({
        selectInput('Adduct','Adduct:',getAdducts())
      })
      
      getMZ <- reactive({
        mz <- input$acc_mass
        addrule <- MZedDB$ADDUCT_FORMATION_RULES[which(MZedDB$ADDUCT_FORMATION_RULES$Name == input$Adduct),]
        mz <- ((mz - addrule$Add[1]) * addrule$Charge[1]) / addrule$xM 
        if (input$Isotope != 'none') {
          isorule <- MZedDB$ISOTOPE_RULES[which(MZedDB$ISOTOPE_RULES$Isotope == input$Isotope),]
          mz <- mz - isorule$Mass.Difference
        }
        return(mz)
      })
      
      getMF <- reactive({
        mz <- getMZ()
        maxi <- c(C = input$Cmax,iC = 0,H = input$Hmax,iH = 0,N = input$Nmax,iN = 0,O = input$Omax,iO = 0,F = 0 ,Na = 0,Si = 0,P = input$Pmax,S = input$Smax,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0)
        res <- generateMF(mz,ppm = input$accuracy,charge = 0,applygr = input$applygr,composition = maxi)
        if (nrow(res) > 0){
          mz <- sapply(as.numeric(res$`m/z`),function(x){
            addrule <- MZedDB$ADDUCT_FORMATION_RULES[which(MZedDB$ADDUCT_FORMATION_RULES$Name == input$Adduct),]
            x <- (x * addrule$xM[1])/addrule$Charge[1] + addrule$Add[1]
            if (input$Isotope != 'none') {
              isorule <- MZedDB$ISOTOPE_RULES[which(MZedDB$ISOTOPE_RULES$Isotope == input$Isotope),]
              x <- x + isorule$Mass.Difference
            }
            return(x)
          })
          ppmErr <- sapply(mz,function(x,y){
            x <- (y - x)/x * 10^6
            return(x)
          },y = input$acc_mass)
          res <- data.frame(MF = res$MF, AccurateMass = res$`m/z`, mz = mz, PPMerror = ppmErr)
          colnames(res)[3:4] <- c('m/z','PPM Error')
          res
        } else {
          res <- data.frame(MF = res$MF, AccurateMass = res$`m/z`, mz = double(), PPMerror = double())
          colnames(res)[3:4] <- c('m/z','PPM Error')
          res
        }
      })
      
      output$mf_table <- renderDataTable({
        res <- getMF()
        return(res)
      },server = T,selection = 'single',rownames = F) 
      
      
      # output$isoDistsPlot <- renderPlot({
      #   r <- input$mf_table_rows_selected
      #   if (!is.null(r)) {
      #     mf <- getMF()
      #     mf <- mf[r,]
      #     if (grepl('i',mf)) {
      #       NULL
      #     } else {
      #       res <- isoDistr(as.character(mf$MF),chrg = input$charge)
      #       plot(res[,1],res[,3],type = "h",xlab = "m/z",ylab = "Relative Intensity",col = "Blue")
      #     }
      #   } else {
      #     NULL
      #   }
      # })
      # 
      # output$isoDistsTable <- renderDataTable({
      #   r <- input$mf_table_rows_selected
      #   if (!is.null(r)) {
      #     mf <- getMF()
      #     mf <- mf[r,]
      #     if (grepl('i',mf)) {
      #       NULL
      #     } else {
      #       res <- isoDistr(as.character(mf$MF),chrg = input$charge)
      #       res <- res[,-2]
      #       colnames(res)[2] <- "Relative Intensity"
      #       res
      #     }
      #   } else {
      #     NULL
      #   }
      # })
      
      getMFmatches <- reactive({
        r <- input$mf_table_rows_selected
        if (!is.null(r)) {
          mf <- getMF()
          mf <- mf$MF[r]
          DB <- MZedDB$MZedDB_ALL
          DB <- DB[which(DB$MF == mf),]
          if (nrow(DB) > 0) {
            DB <- DB[,c(1,2,4,5,9)]
            colnames(DB)[4:5] <- c('Accurate Mass', 'Smile')
            DB
          } else {
            NULL
          }
        } else {
          NULL
        }
      })
      
      output$mfHitsTable <- renderDataTable({
        r <- input$mf_table_rows_selected
          if (!is.null(r)) {
            getMFmatches()
          } else {
            NULL
          }
      },server = T,selection = 'single',rownames = F)
      
      output$MFstructure <- renderPlot({
        r <- input$mfHitsTable_rows_selected
        if (!is.null(r)) {
          res <- getMFmatches()
          res <- res[r,]
          smile <- gsub('"','',res$Smile)
          par(mar = c(0,0,0,0))
          sm <- parse.smiles(smile)[[1]]
          temp1 <- view.image.2d(sm,500,500)
          plot(NA,NA,xlim = c(1,100),ylim = c(1,100),xaxt = 'n',yaxt = 'n',xlab = '',ylab = '')
          rasterImage(temp1,1,1,100,100)
        } else {
          NULL
        }
      },width = 400,height = 400)
      
      ###################### Adduct Calculator ###########################
      
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
      
      ################### Relationship Prediction ############################
      
      output$predictedRelTable <- renderDataTable({
        if (input$mode == 'positive') {
          mode <- 'p'
        }
        if (input$mode == 'negative') {
          mode <- 'n'
        }
        res <- relationshipPredictor(c(input$mz1,input$mz2),mode)
      })
    }
  )
}