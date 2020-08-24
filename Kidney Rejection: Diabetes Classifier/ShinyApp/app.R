#import libraries

library(shiny)
library(shinydashboard)
library(BiocManager)
library(affy)
library(GEOquery)
library(tidyverse)
library(limma)
library(GiANT)
library("annotate")
library("u133x3p.db")
library("hgu133plus2.db")
library(glmnet)
library(ggplot2)
library(ggpubr)

library(BiocManager)
options(repos = BiocManager::repositories())



load("classifiers.RData")
diabetes_top10_genes<-read_rds("diabetes_top10_genes.rds")
rejection_encode_filtered<-readRDS("rejection_encode_filtered.rds")
transplant_top10_genes<-readRDS("transplant_10_genes.rds")
disease<-readRDS("disease.rds")



##pie chart
#number of patients (predicted with) diabetes that have rejection / no rejections 

data=as.data.frame(as.matrix(table(rejection_encode_filtered)))
data$rejection_status=ifelse(rownames(data)=="0","no rejection","rejection")
piedata<-data
bp<-ggplot(piedata, aes(fill=rejection_status,y=V1,x="")) + geom_bar(width =1, stat="identity") +ggtitle("Rejection status for patients predicted with Type 2 Diabetes")+ labs(y = "Percentage",x="Diabetes patients(predicted)")+geom_text(aes(label = paste(round(V1 / sum(V1) * 100, 1), "%")),position = position_stack(vjust = 0.5))+theme_classic() +theme(plot.title = element_text(hjust=0.5),axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 

blank_theme <- theme_minimal()+theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
)

pie<-bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Dark2")+ blank_theme +theme(axis.text.x=element_blank()) 

##scatter plot

genemeans = t(apply(transplant_top10_genes[,-1], 1, mean))


colors <- c("Other patients(mean)" = "black", "You"="red")

### Confusion matrix 


X = as.matrix(t(diabetes_top10_genes))
y = disease
X2 = as.matrix(t(transplant_top10_genes))
y2 = rejection_encode_filtered

draw_confusion_matrix <- function(cmtrx,model_name) {
  
  total <- sum(cmtrx$table)
  res <- as.numeric(cmtrx$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    
    if (amount == 0)
      return("#FFFFFF")
    
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(paste('CONFUSION MATRIX -', model_name), cex.main=2)
  
  # create the matrix
  
  classes = colnames(cmtrx$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1])) # CM
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3])) # CM
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2])) # CM
  rect(250, 305, 340, 365, col=getColor("green", res[4])) # CM
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cmtrx results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics
  
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cmtrx$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cmtrx$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cmtrx$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cmtrx$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cmtrx$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cmtrx$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cmtrx$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cmtrx$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cmtrx$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cmtrx$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information
  text(30, 35, names(cmtrx$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cmtrx$overall[1]), 3), cex=1.4)
  text(70, 35, names(cmtrx$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cmtrx$overall[2]), 3), cex=1.4)
}


set.seed(1)
folds <- sample(rep(1:3, length.out = nrow(X)), size = nrow(X), replace = F)
table(folds)


#rf 3fold-diabetes
CV_rf <- lapply(1:3, function(x){ 
  model <- randomForest::randomForest(x = X[folds != x,], y = as.factor(y[folds != x]))
  preds <- predict(model,  X[folds == x,])
  return(data.frame(preds, real = y[folds == x]))
})

CV_rf <- do.call(rbind, CV_rf)

rf_cm <- caret::confusionMatrix(CV_rf$preds, CV_rf$real)

CM1<-draw_confusion_matrix(rf_cm, "Random forest - Diabetes")

#transplant
CV_svm <- lapply(1:3, function(x){ 
  model <- e1071::svm(x = X2[folds != x,], y = as.factor(y2[folds != x]))
  preds <- predict(model,  X2[folds == x,])
  return(data.frame(preds, real = as.factor(y2[folds == x])))
})

CV_svm <- do.call(rbind, CV_svm)

svm_cm <- caret::confusionMatrix(CV_svm$preds, CV_svm$real)

CM2<-draw_confusion_matrix(svm_cm, "SVM (Transplant)")


##Product



###

ui = dashboardPage(
  dashboardHeader(title = "Kidney Rejection Predictor",titleWidth = 400),
  dashboardSidebar(
    

    sidebarMenu(
      
      menuItem("Information", tabName = "Information", icon = icon("info")),
      menuItem("Data", tabName = "Data", icon = icon("database")),
      menuItem("Classifier", tabName = "Classifier", icon = icon("clipboard")),
      menuItem("Validation", tabName = "Validation", icon = icon("folder-open") )
      
    )
      
    ),
    
    dashboardBody(
      
      tabItems(
        
        tabItem(tabName = "Information",
                tabsetPanel(
                  tabPanel("About this Tool",
                           tags$head(tags$style('h2 {color:cornflowerblue;}')), h2(strong("Hello and welcome to our Kidney Rejection Risk Calculator.")),
                           h4("In this calculator, you will enter your specific Gene Expression values for the Top 10 Genes we have identified as playing a role in Kidney Rejection, you will also do this for Type II Diabetes. "),
                           h4("This calculator will then predict whether or not you will experience Acute Kidney Rejection or acception."),
                           h4("The calculator will also make a prediction on whether you have Type II Diabetes based on the Gene Expression values you have input for the Top 10 Genes provided."),
                           br(),
                           h2(strong("How to navigate")),
                           h4("Head to the Data Tab to enter the required top 10 genes and click on the 'Enter' button to start your prediction."),
                           h4("The validation tab displays the predictors key performance accuracies."),
                           h4("Once you've entered the genes for yourself or for your patients, head to the Classifier Tab and switch between the Transplant Predictor and Diabetes Predictor for your prediction results!  "),
                           h4("There you'll also find some informative charts and plots that will help you understand your prediction results.")
                  ),
                  tabPanel("Methods",
                           
                             h4(code("GSE20966"),"-This dataset was used to train the diabetes classifier. It included more than 61 000 gene outputs for 20 patients. 10 patients in the dataset had Type 2 Diabetes and 10 did not"),
                             h4(code("GSE14346"),"-This dataset was used to train the kidney rejection classifier.It included more than 54 000 genes for 75 patients. 38 patients had acute kidney reject and 37 had no rejection."),
                             br(),
                             h4(strong("Step1:"),"Pre processed data sets - normalising, filtering out NA’s. "),
                             h4(strong("Step2:"),"Filtered common genes present in both data sets."),
                             h4(strong("Step3:"),"Filtered multiple probes by finding the maximum"),
                             h4(strong("Step4:"),"Filtered both datasets to top 10 genes with greatest significance."),
                             h4(strong("Step5:"),"Trained and tested Random Forest (RF) Type 2 Diabetes Classifier."),
                             h4(strong("Step6:"),"Filtered kidney dataset based on the Diabetes classifier."),
                             h4(strong("Step7:"),"Trained and tested Support Vector Machine (SVM) Kidney Transplant classifier specifically for Type 2 Diabetes Patients.")
                           
                             
                  )
          
        )
        
        
          
          
        ),
        
        tabItem(tabName = "Data",
            tabsetPanel(
                tabPanel("What to Upload",
                    h3(strong("This tool requires the input of the values of 10 different genes explained below:")),
                    
                    #h2(strong("Input Genes")), #Kidney Transplant Rejection and Type II Diabetes Predictor
                    h3("Kidney Transplant Rejection Predictor and Type II Diabetes Predictor draws and predicts on the same top 10 genes."),
                    h4(""),
                    h4(strong("CCDC69:")," A Protein Coding Gene"),
                    h4(strong("RASGRF1:"),"The protein encoded by this gene is a guanine nucleotide exchange factor (GEF) similar to the Saccharomyces cerevisiae CDC25 gene product. Functional analysis has demonstrated that this protein stimulates the dissociation of GDP from RAS protein. "),
                    h4(strong("MDFIC:"),"This gene product is a member of a family of proteins characterized by a specific cysteine-rich C-terminal domain, which is involved in transcriptional regulation of viral genome expression."),
                    h4(strong("FXYD3:"),"This gene belongs to a small family of FXYD-domain containing regulators of Na+/K+ ATPases which share a 35-amino acid signature sequence domain, beginning with the sequence PFXYD, and containing 7 invariant and 6 highly conserved amino acids. This gene encodes a cell membrane protein that may regulate the function of ion-pumps and ion-channels."),
                
                    h4(strong("SLC2A2:"),"This gene encodes an integral plasma membrane glycoprotein of the liver, islet beta cells, intestine, and kidney epithelium."),                  
                    h4(strong("TSHR:"),  "Thyroid Stimulating Hormone Receptor is a Protein Coding gene."),
                    h4(strong("IGFBP3:"), "This gene is a member of the insulin-like growth factor binding protein (IGFBP) family and encodes a protein with an IGFBP domain and a thyroglobulin type-I domain"),
                    h4(strong("COA3:"), "This gene encodes a member of the cytochrome c oxidase assembly factor family."),
                    h4(strong("PCOLCE2:"), "Procollagen C-Endopeptidase Enhancer 2 is a Protein Coding gene. "),
                    h4(strong("NEFL:"), "Neurofilaments are type IV intermediate filament heteropolymers composed of light, medium, and heavy chains. Neurofilaments comprise the axoskeleton and they functionally maintain the neuronal caliber. They may also play a role in intracellular transport to axons and dendrites. This gene encodes the light chain neurofilament protein. "),
                    
                    
                    
            
                    
                                               
                ), tabPanel("Upload Information",
                            fluidRow(
                              
                              box( title = "Kidney Transplant Prediction Gene Values",  status = "warning", solidHeader = T,
                                   
                                #h4(strong("Kidney Transplant Prediction Gene Values")),      
                                numericInput(inputId = "gene1", label = "CCDC69 gene value", value = 0),  
                                numericInput(inputId = "gene2", label = "RASGRF1 gene value", value = 0),
                                numericInput(inputId = "gene3", label = "MDFIC gene value", value = 0),
                                numericInput(inputId = "gene4", label = "FXYD3 gene value", value = 0),
                                numericInput(inputId = "gene5", label = "SLC2A2 gene value", value = 0),
                                numericInput(inputId = "gene6", label = "TSHR gene value", value = 0),
                                numericInput(inputId = "gene7", label = "IGFBP3 gene value", value = 0),
                                numericInput(inputId = "gene8", label = "COA3 gene value", value = 0),
                                numericInput(inputId = "gene9", label = "PCOLCE2 gene value", value = 0),
                                numericInput(inputId = "gene10", label = "NEFL gene value", value = 0),
                                
                                actionButton(inputId = "enter", label = "Enter genes and predict for rejection"),
                                
                                width  = 6
                              ),
                              
                              box(title= "Diabetes Prediction",  status = "warning", solidHeader = T,
                                #h4(strong("Diabetes Prediction")), 
                              
                                actionButton(inputId = "denter", label = "Predict for Diabetes"),
                                br(),
                                tags$em("Note: The same top 10 genes are used for the prediction of both patient's Diabetes (Type II) and Kidney Rejection status. "),
                                
                                width  = 6
                              )
                                
                                
                                
                              )
                              
                              
                            
                        
                        
                       
                        
                        
               )
        
          )
          
        ), 
        tabItem(tabName = "Classifier",
                   
        
        
            
            tabsetPanel(
       
              tabPanel("Transplant Predictor",
                       #box( title = "Warning",width = 10,background = "maroon",
                       #     h5('Please enter your Gene Data to output the Gene comparison plot!')
                       #),
                       textOutput(outputId = "warning"),
                       fluidRow(
                         valueBoxOutput(outputId = "class"),
                       
                       box(title = "Pie Chart",  status = "primary", solidHeader = T,
                           width = 10, collapsible = T,
                           plotOutput("pie")),
                       
                       
                       
                       box(
                         title = "Gene Comparison Plot",  status = "primary",solidHeader = T,
                         width = 10, collapsible = T,
                         plotOutput("scatter_transplant")
                       )
                      
                  
                       )
                       ),
              tabPanel( "Diabetes Predictor", 
                        textOutput(outputId = "warning2"),
                        #box( title = "Warning",width = 10,background = "maroon",
                         #    h5('Please enter your Gene Data to output the Gene comparison plot!')
                       # ),
                        fluidRow(
                            valueBoxOutput(outputId = "diab")
                        ),
                        fluidRow( 
                           box( background = "maroon",
                             h5('This Classifier is Not for Diagnosis Purposes!')
                                ),
                        ),
                        
                        fluidRow( box(
                            title = "Gene Comparison Plot", status = "primary",
                            solidHeader = TRUE,
                            width = 10, collapsible = T,
                            fluidRow(
                            plotOutput("scatter_diabetes")
                            )
                          )),
                        
                       
                          
                          
                          #FR
                                
                                )
              
          
          ), height = 30
       
        
        
        
                   
          ),
        tabItem(tabName = "Validation",
                
                fluidRow(
                  
                  box(title = "Confusion Matrix-Random Forest-Transplant Predictor",
                      status = "primary", solidHeader = TRUE,
                      width = 9,collapsible = T,
                      fluidRow(
                        plotOutput(outputId = "CM2")
                        
                      )),
                 # box(
                  #  title = "Confusion Matrix-Random Forest-Transplant Predictor", solidHeader = T,
                  #  width = 8, collapsible = T,
                 #   plotOutput(outputId = "CM2")
                 # ),
                 
                 
                   box(title = "Confusion Matrix-Random Forest-Diabetes Predictor",
                       status = "primary", solidHeader = TRUE,
                       width = 9,collapsible = T,
                       fluidRow(
                         plotOutput(outputId = "CM")
                         
                       )),
                  #box(
                  #  title = "Confusion Matrix-Random Forest-Diabetes Predictor", solidHeader = T,
                  #  width = 8, collapsible = T,
                 
                  #),
                  
                
                                     
                  
                  
                  
                ),
                
                
                
                
                
                
        )
       
    
    
    
    )
  
  )
)
  
  

server = function(input, output) {
 
  
  t_results = observeEvent(input$enter, {
    
    
    values = matrix(c(input$gene1, input$gene2, input$gene3, input$gene4, input$gene5, input$gene6, input$gene7, input$gene8, input$gene9, input$gene10), nrow = 1)
    
    print(1)
    results = predict(t_svm_res, values)
    if(results == 1){
      
      output$class = renderValueBox({
        valueBox("Rejection Predicted", "Kidney Transplant Prediction",icon = icon("exclamation-triangle"), color = "yellow")
        } )
      
    }
    
    else if(results == 0){
      
      output$class = renderValueBox({
        valueBox("No Rejection Predicted","Kidney Transplant Prediction", icon = icon("exclamation-triangle"), color = "yellow")
        } )
      
    }
  }
    
  )
  
  
    di_results = observeEvent(input$denter, {
      #output$warning = renderText("hhhihihi")
     
      values = matrix(c(input$gene10, input$gene4, input$gene7, input$gene3, input$gene2, input$gene5, input$gene9, input$gene8, input$gene6, input$gene1), nrow = 1)
      
      
      dresults = predict(d_rf_res, values)
      print(dresults)
      if(dresults == "non-diabetic control"){

        output$diab = renderValueBox({
          valueBox("Negative","Type II Diabetes Prediction", icon = icon("exclamation-triangle"), color = "yellow")
          } )
        
      }
      
      else if(dresults == "type 2 diabetes"){
        
        output$diab = renderValueBox({
          valueBox("Positive","Type II Diabetes Prediction", icon = icon("exclamation-triangle"), color = "yellow")
          
      
          } )
        
      }

    
    
  }
  )                       
                           
    observeEvent(input$denter,{
      
       
      showNotification(
        #title = "Notification",
        "You've entered your gene values!"
      )
    })  
    
    observeEvent(input$enter,{
      
      
      showNotification(
        #title = "Notification",
        "You've entered your gene values!"
      )
    })  
   
    

  output$class = renderText({"Please Input Your Gene Data In The Data Tab"})

  output$pie <- renderPlot({
    
      bp<-ggplot(piedata, aes(fill=rejection_status,y=V1,x="")) + geom_bar(width =1, stat="identity") +ggtitle("Rejection status for patients predicted with Type 2 Diabetes")+ labs(y = "Percentage",x="Diabetes patients(predicted)")+geom_text(aes(label = paste(round(V1 / sum(V1) * 100, 1), "%")),position = position_stack(vjust = 0.5))+theme_classic() +theme(plot.title = element_text(hjust=0.5),axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
    
      piechart<- bp+coord_polar("y", start=0) + scale_fill_brewer(palette="Dark2")+ blank_theme +theme(axis.text.x=element_blank(),panel.background = element_rect(fill = "lightyellow")) +labs(fill = "Rejection Status")
      piechart
    
    
    
    
  })
 
  
  t_scatter = observeEvent(input$enter, {
    
    
    df <- data.frame(gene_names = rownames(transplant_top10_genes),
                     values = as.numeric(genemeans),
                     inputvalues=as.numeric(c(input$gene1, input$gene2, input$gene3, input$gene4, input$gene5, input$gene6, input$gene7, input$gene8, input$gene9, input$gene10)))
    
  output$scatter_transplant <- renderPlot({
  
  ggplot(df, aes(x = gene_names))+ 
    geom_point(aes(x = gene_names,y = values,color = "Other patients(mean)"),size = 1.5) +
    geom_point(aes(x = gene_names,y = inputvalues,color = "You" , shape =17),size = 1.5)+
    scale_shape_identity()+
    ggtitle("Gene comparisons")+
    theme(plot.title = element_text(hjust = 0.5))+  
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.title = element_text( size = 12),panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                                                              size = 2, linetype = "solid")
    )+ labs(color='Legend')  +xlab("Genes") + ylab("Values")
  
})
  
})
  
  
  d_scatter = observeEvent(input$denter, {
    
    
    df <- data.frame(gene_names = rownames(diabetes_top10_genes),
                     values = as.numeric(genemeans),
                     inputvalues=as.numeric(c(input$gene10, input$gene4, input$gene7, input$gene3, input$gene2, input$gene5, input$gene9, input$gene8, input$gene6, input$gene1)))
    
    output$scatter_diabetes <- renderPlot({
      
      ggplot(df, aes(x = gene_names))+ 
        geom_point(aes(x = gene_names,y = values,color = "Other patients(mean)"),size = 1.5) +
        geom_point(aes(x = gene_names,y = inputvalues,color = "You" , shape =17),size = 1.5)+
        scale_shape_identity()+
        ggtitle("Gene comparisons")+
        theme(plot.title = element_text(hjust = 0.5))+  
        theme(
          legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top"),
          legend.title = element_text( size = 12),panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                                                                  size = 2, linetype = "solid")
        )+ labs(color='Legend')  +xlab("Genes") + ylab("Values")
      
    })
    
  })  
  
  
  output$userpanel = renderUI({
    sidebarUserPanel(
      span("test")
    )
    
    
  })
  
  output$CM <-renderPlot({
    
    draw_confusion_matrix(rf_cm, "Random forest - Diabetes")
  })
  
  output$CM2 <-renderPlot({
    
    draw_confusion_matrix(svm_cm, "SVM -Transplant")
  })
  
}

shinyApp(ui = ui, server = server)