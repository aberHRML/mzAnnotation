#' shinyMZedDB
#' @export

shinyMZedDB <- function(){
  shinyApp(
    ui = navbarPage("Shiny-MZedDB",
                    tabPanel("Putative Ionisation Product",
                             # Create a new Row in the UI for selectInputs
                             fluidRow(
                               column(3, numericInput("acc_mz","Accurate Mass:",117.078979)),
                               column(2, numericInput("ppm","PPM:",1)),
                               column(2, checkboxInput("mode_p", "Positive Mode",value=F)),
                               column(2, checkboxInput("mode_ne", "Neutral",value=T)),
                               column(2, checkboxInput("mode_n", "Negative Mode",value=F)),
                               column(2, checkboxInput("pip_iso", "Isotopes",value=F))
                               #column(4, selectInput('add_sel_i', 'Adducts:', adducts, multiple=TRUE, selectize=FALSE))                                                    
                             ),
                             fluidRow(
                               dataTableOutput(outputId="pip_table")
                             )
                    ),
                    tabPanel("Molecular Formula Generator",
                             fluidRow(
                               column(3, numericInput("acc_mass","Accurate Mass:",117.078979)),
                               column(2, numericInput("accuracy","PPM:",1)),
                               column(2, numericInput("charge", "Charge:",0)),
                               column(2, checkboxInput("applygr", "Apply 7GR",value=T))
                             ),
                             # new row for elemental composition
                             fluidRow(
                               column(2, numericInput("C_min","C:",0),numericInput("C_max","",10)),
                               column(2, numericInput("H_min","H:",0),numericInput("H_max","",20)),
                               column(2, numericInput("N_min","N:",0),numericInput("N_max","",5)),
                               column(2, numericInput("O_min","O:",0),numericInput("O_max","",5)),
                               column(2, numericInput("P_min","P:",0),numericInput("P_max","",0)),
                               column(2, numericInput("S_min","S:",0),numericInput("S_max","",0)),
                               column(2, numericInput("Na_min","Na:",0),numericInput("Na_max","",0)),
                               column(2, numericInput("K_min","K:",0),numericInput("K_max","",0)),
                               column(2, numericInput("Cl_min","Cl:",0),numericInput("Cl_max","",0)),
                               column(2, numericInput("iC_min","iC:",0),numericInput("iC_max","",0)),
                               column(2, numericInput("iO_min","iO:",0),numericInput("iO_max","",0)),
                               column(2, numericInput("iK_min","iK:",0),numericInput("iK_max","",0)),
                               column(2, numericInput("iCl_min","iCl:",0),numericInput("iCl_max","",0))
                             ),
                             # Create a new row for the table.
                             fluidRow(
                               dataTableOutput(outputId="mf_table")
                             )
                    ),
                    tabPanel("Isotope Distributions",
                             # Create a new Row in the UI for selectInputs
                             fluidRow(
                               column(3, textInput("mf","Molecular Formula:","C4H5O5")),
                               column(2, numericInput("chrge","Charge:",-1))
                             ),
                             # Create a new row for the plot.
                             fluidRow(
                               plotOutput(outputId = "plot_rel")
                             ),
                             # Create a new row for the table.
                             fluidRow(
                               dataTableOutput(outputId="iso_table")
                             ) 
                             
                    ),
                    tabPanel("Draw Smiles",
                             sidebarPanel(
                               textInput("smile","Smiles:","O[C@@H](CC([O-])=O)C([O-])=O")
                             ),
                             mainPanel(
                               plotOutput(outputId = "structure"),width=5
                             )
                    ),
                    tabPanel("MZedDB",
                             fluidRow(
                               dataTableOutput(outputId="MZDB")
                             ) 
                    ),
                    tabPanel("Adduct Calculator",
                              fluidRow(
                               column(3, textInput("id","Database ID:","D29160"))
                             ),
                             fluidRow(
                               dataTableOutput(outputId="entry")
                             ),
                             fluidRow(
                               tableOutput("adduct")
                             ) 
                    )
    ),
    server = function(input, output) {
      
      # Filter data based on selections
      output$pip_table <- renderDataTable({
        if(input$mode_p){
          mode <- "p"
        }
        if(input$mode_n){
          mode <- "n"
        }
        if(input$mode_ne){
          mode <- "ne"
        }
        if(input$pip_iso){
          iso <- c('C13','O18','S34')
        } else {
          iso <- NULL
        }
        pip_tab <- getPIP(input$acc_mz,mode,input$ppm,iso=iso)
        #pip_tab <- pip_tab[,-c(3,6,7,8,10)]
        #colnames(pip_tab)[3:8] <- c("MF","Accurate Mass","Smiles","Adduct","Adduct m/z","PPM Error")
        pip_tab
      }) 
      
      output$mf_table <- renderDataTable({
        maxi <- c(C=input$C_max,iC=input$iC_max,H=input$H_max,iH=0,N=input$N_max,iN=0,O=input$O_max,iO=input$iO_max,F=0,Na=input$Na_max,Si=0,P=input$P_max,S=input$S_max,Cl=input$Cl_max,iCl=input$iCl_max,Br=0,iBr=0,K=input$K_max,iK=input$iK_max)
        mini <- c(C=input$C_min,iC=input$iC_min,H=input$H_min,iH=0,N=input$N_min,iN=0,O=input$O_min,iO=input$iO_min,F=0,Na=input$Na_min,Si=0,P=input$P_min,S=input$S_min,Cl=input$Cl_min,iCl=input$iCl_min,Br=0,iBr=0,K=input$K_min,iK=input$iK_min)
        prec <- input$accuracy/10^6*input$acc_mass*1000
        res <- mfGen(input$acc_mass,maxi,mini,prec,input$charge,applygr=input$applygr)
        if(length(res)==0){
          res <- matrix(nrow=0,ncol=5)
          colnames(res) <- c("Clean MF", "MF", "RDB", "m/z", "PPM Error")
        } else{
          res <- data.frame(matrix(unlist(res), nrow = length(res),byrow=T))
          res <- res[,-c(5,6)]
          colnames(res) <- c("Clean MF", "MF", "RDB", "m/z","PPM Error")
          res[,5] <- sapply(as.numeric(as.character(res[,4])),function(x,mass){x <- as.numeric(x);x <- (x-mass)/mass*10^6; return(round(x,5))},mass=input$acc_mass)
        }   
        res
      }) 
      
      output$iso_table <- renderDataTable({
        res <- isoDistr(input$mf,chrg = input$chrge)
        res <- res[,-2]
        colnames(res)[2] <- "Relative Intensity"
        res
      })
      
      output$plot_rel <- renderPlot({
        res <- isoDistr(input$mf,chrg = input$chrge)
        plot(res[,1],res[,3],type="h",xlab="m/z",ylab="Relative Intensity",col="Blue")
      })
      
      output$structure <- renderPlot({
        par(mar=c(0,0,0,0))
        sm <- parse.smiles(input$smile)[[1]]
        temp1 <- view.image.2d(sm,500,500)
        plot(NA,NA,xlim=c(1,100),ylim=c(1,100),xaxt='n',yaxt='n',xlab='',ylab='') 
        rasterImage(temp1,1,1,100,100) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
      })
      
      output$MZDB <- renderDataTable({
        data('MZedDB') 
        res <- MZedDB$MZedDB_ALL
        res <- res[,-c(3,6,7,8,10)]
        names(res)[4:5] <- c('Accurate Mass', 'SMILE')
        res
      })
      
      output$entry <- renderDataTable({
        data('MZedDB') 
        res <- MZedDB$MZedDB_ALL
        res <- res[,-c(3,6,7,8,10)]
        names(res)[4:5] <- c('Accurate Mass', 'SMILE')
        res <- res[which(res$ID==input$id),]
        res
      })
      
      output$adduct <- renderTable({
        data('MZedDB') 
        id <- MZedDB$MZedDB_ALL
        id <- id[,-c(3,6,7,8,10)]
        names(id)[4:5] <- c('Accurate Mass', 'SMILE')
        id <- id[which(id$ID==input$id),]
        
        metrules <- MZedDB$MZedDB_METRULES
        metrules <- metrules[which(metrules$ID==id$ID),]

        rules <- MZedDB$ADDUCT_FORMATION_RULES
        res <- apply(rules,1,function(rule,id,metrules){
          Nch <- metrules$Nch
          Nacc <- metrules$Nacc
          Ndon <- metrules$Ndon
          Nnhh <- metrules$Nnhh
          Noh <-  metrules$Noh
          Ncooh <- metrules$Ncooh
          Ncoo <-  metrules$Ncoo
          if(eval(parse(text=as.character(rule[6])))){
            mz <- id$`Accurate Mass`*as.numeric(rule[3])
            if (as.numeric(rule[2]) != 0){
              mz <- mz/as.numeric(rule[2])
            }
            mz <- mz + as.numeric(rule[4])
            mz <- round(mz,5)
          } else {
            mz <- 'Not Possible'
          }
          return(mz)
        },id=id,metrules=metrules)
        res <- cbind(rules$Name,res)
        colnames(res) <- c('Adduct','m/z')
        #res <- data.frame(Adduct=rules$Name, `m/z`=res)
        res
      })
    }
    )
}