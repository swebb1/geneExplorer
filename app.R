#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#

library(shiny)
library(shinydashboard)
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(DT)

a<-readRDS("annotation.rds") ###Sort out factors
r<-readRDS("rna_data.rds") ### sort out missing data and merged data
p<-readRDS("prot.rds") 

cols<-c(brewer.pal(n = 9,name = "Set1")[3],"grey39")
cols2<-c(brewer.pal(n = 9,name = "Set1")[2],"grey39")

ui<-shinyUI(dashboardPage(
  title="Gene explorer",
  skin="blue",
  dashboardHeader(title = "Gene Explorer", titleWidth = 400),
  dashboardSidebar(width=220,
                   sidebarMenu(
                     #sidebarSearchForm(textId = "searchText", buttonId = "searchButton",label = "Search..."),
                     menuItem("Home",tabName="home",icon=shiny::icon("home")),
                     menuItem("Multiple Gene Analysis",tabName="multi",icon=shiny::icon("chart-bar")),
                     menuItem("Data Table",tabName="data",icon=shiny::icon("database"))
                   )
  ),
  dashboardBody(
    includeCSS("www/custom.css"),
    tabItems(
      tabItem(tabName="home",
              fluidRow(
                box(
                  title="Gene Search",width = 9,status="primary",solidHeader=F,
                  div(style = 'overflow-x: scroll', dataTableOutput('anno'))
                ),
                box(
                  title="Selected Gene",width = 3,status="primary",solidHeader=F,
                  div(style = 'overflow-x: scroll', dataTableOutput('linksTable'))
                ),
                box(
                  title="RNA-seq Data",width = 6,status="primary",solidHeader=F,
                  numericInput("sig",label = "Significance threshold",value = 0.05,min = 0,max = 1,width=100),
                  plotOutput("rnaPlot")
                ),
                box(
                  title="Proteomics Data",width = 9,status="primary",solidHeader=F,
                  numericInput("sig2",label = "Significance threshold",value = 0.05,min = 0,max = 1,width=100),
                  div(style = 'overflow-x: scroll', dataTableOutput('protTable'))
                ),
                box(
                  title="Selected Protein",width = 3,status="primary",solidHeader=F,
                  div(style = 'overflow-x: scroll', dataTableOutput('plinksTable'))
                ),
                box(
                  title="Proteomics Plot",width = 6,status="primary",solidHeader=F,
                  #uiOutput("protSelect"),
                  #uiOutput("pepSelect"),
                  plotOutput("protPlot")
                )#,
                #box(
                #  title="Download table",width = 12,status="primary",solidHeader=F
                #  #uiOutput("downloadFiles")
                #)
              )
      ),
      tabItem(tabName="multi",
              fluidRow(
                box(
                  title = "Multiple Gene Analysis", width = 12, status = "primary",solidHeader=F,
                  textAreaInput("multi",label = "Paste gene list",resize = "both",height = 250),
                  plotOutput("multiPlot")
                  #div(style = 'overflow-x: scroll', dataTableOutput('table'))
                )
              )
      )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output){
  
  output$anno = DT::renderDT({
    return(DT::datatable(a,selection = 'single',filter="bottom",
                         options = list(bSortClasses = TRUE,
                                        aLengthMenu = c(1,5,10,20,50),
                                        pageLength = 5
                                        )))
  })

  genes<-reactive({
    if(length(input$anno_rows_selected)==0){
      genes<-a[input$anno_rows_all[1],]
    }
    else{
      genes<-a[input$anno_rows_selected[1],]
    }
    return(genes)
  })
  
  peptides<-reactive({
    gene<-genes()
    ps<-subset(p,p$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    return(ps)
  })
  
  pep<-reactive({
    peptides<-peptides()
    if(length(input$protTable_rows_selected)==0){
      peps<-peptides[input$anno_rows_all[1],]
    }
    else{
      peps<-peptides[input$protTable_rows_selected[1],]
    }
    return(peps)
  })
  
  output$linksTable = DT::renderDT({
    gene<-genes()
    ensembl<-paste0("<a href='https://www.ensembl.org/Homo_sapiens/Gene/Idhistory?g=",gene$Gene.ID,"' target='_blank'>Ensembl</a>")
    geneC<-paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene$Gene.ID,"' target='_blank'>Gene Cards</a>")
    hpa<-paste0("<a href='https://www.proteinatlas.org/",gene$Gene.ID,"' target='_blank'>Human Protein Atlas</a>")
    
    df<-data.frame(links=c(ensembl,geneC,hpa))
    names(df)<-gene$Gene.Symbol
    return(DT::datatable(df,escape = F,options=list(paging=F,searching=F)))
  })
  
  output$plinksTable = DT::renderDT({
    pep<-pep()
    uniprot<-paste0("<a href='https://www.uniprot.org/uniprot/",unique(pep$Uniprot.ID),"' target='_blank'>Uniprot</a>")
    pp2a<-paste0("<a href='http://proviz.ucd.ie/proviz.php?uniprot_acc=",unique(pep$Uniprot.ID),"&tools=pp2a' target='_blank'>PP2A</a>")
    grnai<-paste0("<a href='http://www.genomernai.org/v17/geneSearch/%5B%22",unique(pep$Uniprot.ID),"%22%5D'  target='_blank'>Genome RNAi</a>" )
    
    df<-data.frame(links=c(uniprot,pp2a,grnai))
    names(df)<-unique(pep$Uniprot.ID)
    return(DT::datatable(df,escape = F,options=list(paging=F,searching=F)))
  })

  output$rnaPlot <- renderPlot({
    sig=input$sig
    values=cols
    names(values)<-c(paste0(" Significant <= ",sig),paste0(" > ",sig))
    gene<-genes()
    rs<-subset(r,r$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    if(dim(rs)[1]==0){
      return(NULL)
    }
    g1<-ggplot(rs %>% mutate(fill = ifelse(padj<=sig,paste0(" Significant <= ",sig),paste0(" > ",sig))),
        aes(text=paste("padj:",padj)))+
      geom_col(aes(x=Experiment,y=log2FoldChange,fill=fill))+
      scale_fill_manual(values = values)+
      theme_bw()+
      theme(text = element_text(size=20),legend.title = element_blank())+
      ggtitle(as.character(gene$Gene.Symbol))+
      facet_grid(Lab~.,scales="free",space="free")+
      coord_flip()
    g1
  })
  
  # output$protSelect<-renderUI({
  #   ps<-peptides()
  #   choices<-unique(ps$Uniprot.ID)
  #   selectInput("protSelect",label="Select Protein:", choices=choices,multiple = F,width = 200)
  # })
  # 
  # output$pepSelect<-renderUI({
  #   ps<-peptides()
  #   choices<-unique(ps$ID)
  #   selectInput("pepSelect",label="Select Peptide:", choices=choices,multiple = F,width = 200)
  # })
  
  output$protTable <- DT::renderDT({
    ps<-peptides()
    if(dim(ps)[1]==0){
      return(NULL)
    }
    dt<-DT::datatable(ps,selection="single",filter="bottom",
                         options = list(bSortClasses = TRUE,
                                        aLengthMenu = c(1,5,10,20,50),
                                        pageLength = 50
                         ))
    dt<-dt %>% formatStyle(
      'pval',
      target = 'row',
      backgroundColor = styleInterval(c(0.05), c(cols[1], cols[2]))
    )
    dt
  })
    
      
  output$protPlot <- renderPlot({
    sig2=input$sig2
    values=cols
    names(values)<-c(paste0(" Significant <= ",sig2),paste0(" > ",sig2))
    gene<-genes()
    ps<-peptides()
    if(dim(ps)[1]==0){
      return(NULL)
    }
    pep<-pep()
    pss<-subset(ps,ps$Uniprot.ID==unique(pep$Uniprot.ID) & ps$ID==unique(pep$ID))
    seq<-unique(pss$Sequence)
    g1<-ggplot(pss %>% mutate(fill = ifelse(pval<=sig2,paste0(" Significant <= ",sig2),paste0(" > ",sig2))),
               aes(text=paste("pval:",pval)))+
      geom_col(aes(x=Experiment,y=FoldChange,fill=fill))+
      scale_fill_manual(values = values)+
      theme_bw()+
      theme(text = element_text(size=20),legend.title = element_blank())+
      ggtitle(paste(gene$Gene.Symbol,seq))+
      facet_grid(ID~Uniprot.ID)+#,scales="free",space="free")+
      coord_flip()
    ##NEED SEQUENCE SOMEWHERE!!
    g1
  })
  
  output$multiPlot <- renderPlot({
    sig=input$sig
    genes<-strsplit(input$multi,"\n")[[1]]
    if(length(genes>=1)){
      rs<-subset(r,r$Gene.ID %in% genes) %>% spread(measurement,value)
      
      g1<-ggplot(rs,aes(x=Experiment,y=log2FoldChange))+
        geom_boxplot(fill=cols2[1])+geom_jitter(colour=cols2[2])+
        theme_bw()+
        theme(text = element_text(size=20),legend.title = element_blank())+
        ggtitle("Gene Set")+
        facet_grid(Lab~.,scales="free",space="free")+
        coord_flip()
      return(g1)
    }
    else{
      return(NULL)
    }
  })
  
  
}

shinyApp(ui = ui, server = server)




