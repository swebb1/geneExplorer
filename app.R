#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(DT)

a<-readRDS("../data/annotation.rds") ###Sort out factors
r<-readRDS("../data/rna_data.rds") ### sort out missing data

cols<-c(brewer.pal(n = 9,name = "Set1")[3],"grey39")

ui<-shinyUI(dashboardPage(
  title="Gene explorer",
  skin="red",
  dashboardHeader(title = "Gene Explorer", titleWidth = 400),
  dashboardSidebar(width=220,
                   sidebarMenu(
                     #sidebarSearchForm(textId = "searchText", buttonId = "searchButton",label = "Search..."),
                     menuItem("Home",tabName="home",icon=shiny::icon("home")),
                     menuItem("Multiple Gene Analysis",tabName="multi",icon=shiny::icon("home")),
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
                  title="RNA-seq Data",width = 12,status="primary",solidHeader=F,
                  numericInput("sig",label = "Significance threshold",value = 0.05,min = 0,max = 1,width=100),
                  plotOutput("rnaPlot")
                ),
                box(
                  title="Spencer Data",width = 12,status="primary",solidHeader=F,
                  plotlyOutput("spencerPlot")
                ),
                box(
                  title="Ly/Barr Data",width = 12,status="primary",solidHeader=F
                  #plotOutput("lyPlot")
                ),
                box(
                  title="Download table",width = 12,status="primary",solidHeader=F
                  #uiOutput("downloadFiles")
                )
              )
      ),
      tabItem(tabName="multi",
              fluidRow(
                box(
                  title = "Multiple Gene Analysis", width = 12, status = "primary",solidHeader=TRUE
                  #div(style = 'overflow-x: scroll', dataTableOutput('table'))
                )
              )
      ),
      tabItem(tabName="data",
              fluidRow(
                box(
                  title = "Data Table", width = 12, status = "primary",solidHeader=TRUE
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
                                        aLengthMenu = c(5,10,20,50),
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
  
  output$linksTable = DT::renderDT({
    gene<-genes()
    geneC<-paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene$Gene.ID,"' target='_blank'>Gene Cards</a>")
    hpa<-paste0("<a href='https://www.proteinatlas.org/",gene$Gene.ID,"' target='_blank'>Human Protein Atlas</a>")
    df<-data.frame(links=c(geneC,hpa))
    names(df)<-gene$Gene.Symbol
    return(DT::datatable(df,escape = F,options=list(paging=F,searching=F)))
  })

  output$rnaPlot <- renderPlot({
    sig=input$sig
    values=cols
    names(values)<-c(paste0(" Significant <= ",sig),paste0(" > ",sig))
    gene<-genes()
    rs<-subset(r,r$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    
    g1<-ggplot(rs %>% mutate(fill = ifelse(padj<=sig,paste0(" Significant <= ",sig),paste0(" > ",sig))),
        aes(text=paste("padj:",padj)))+
      geom_col(aes(x=Experiment,y=log2FoldChange,fill=fill))+
      scale_fill_manual(values = values)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.title = element_blank())+
      ggtitle(as.character(gene$Gene.Symbol))+
      facet_grid(Lab~.,scales="free",space="free")+
      coord_flip()
    g1
  })
  
  output$gilPlot <- renderPlotly({
    gene<-genes()
    gs<-subset(g,g$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    gs$adj.p.value<-factor(gs$adj.p.value)
    gs$adj.p.value<-fct_recode(gs$adj.p.value,`Significant <=0.05`="1",`>0.05`="0")
    
    g1<-ggplot(gs,aes(x=Experiment,y=log2FoldChange,fill=adj.p.value,text=paste("padj:",padj,"baseMean:",baseMean)))+
      geom_col(position = position_dodge(width=1))+
      theme_bw()+scale_fill_manual(values = c(`Significant <=0.05`=cols[1],`>0.05`=cols[2]))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))+
      ggtitle(as.character(gene$Gene.Symbol))
    ggplotly(g1)
  })

  output$spencerPlot <- renderPlotly({
    gene<-genes()
    ss<-subset(s,s$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    ss$adj.p.value<-factor(ss$adj.p.value)
    ss$adj.p.value<-fct_recode(ss$adj.p.value,`Significant <=0.05`="1",`>0.05`="0")
    
    g2<-ggplot(ss,aes(x=Experiment,y=log2FoldChange,fill=adj.p.value,text=paste("padj:",padj)))+
      geom_col(position = position_dodge2(width=1,preserve = "single"))+
      theme_bw()+scale_fill_manual(values = c(`Significant <=0.05`=cols[1],`>0.05`=cols[2]))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))+
      ggtitle(as.character(gene$Gene.Symbol))
    ggplotly(g2)
  })
}

shinyApp(ui = ui, server = server)




