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

a<-readRDS("../data/annotation.rds")
g<-readRDS("../data/gil_data.rds")
s<-readRDS("../data/spencer_data.rds")

cols<-brewer.pal(n = 9,name = "Set1")[c(3,1)]

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
                  title="Gene Information",width = 12,status="primary",solidHeader=F,
                  textInput("gene",label = "Gene",value = a$Gene.ID[1]),
                  div(style = 'overflow-x: scroll', dataTableOutput('linksTable'))
                ),
                box(
                  title="Gil Data",width = 12,status="primary",solidHeader=F,
                  plotlyOutput("gilPlot")
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
server <- function(input, output) {
  
  output$linksTable = renderDataTable({
    mtcars
    gene<-input$gene
    df<-subset(a,a$Gene.ID %in% gene)
    df$Gene_Card<-paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene,"' target='_blank'>Gene Cards</a>")
    return(datatable(df,escape = F))
  },options = list(bSortClasses = TRUE,iDisplayLength = 5)
  )
  
  output$gilPlot <- renderPlotly({
    #
    gene<-input$gene
    gs<-subset(g,g$Gene.ID %in% gene) %>% spread(measurement,value)
    gs$adj.p.value<-factor(gs$adj.p.value)
    gs$adj.p.value<-fct_recode(gs$adj.p.value,`Significant <=0.05`="1",`>0.05`="0")
    
    g1<-ggplot(gs,aes(x=Experiment,y=log2FoldChange,fill=adj.p.value,text=paste("padj:",padj,"baseMean:",baseMean)))+geom_bar(stat="identity")+
      theme_bw()+scale_fill_manual(values = c(`Significant <=0.05`=cols[1],`>0.05`=cols[2]))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
    ggplotly(g1)
  })

  output$spencerPlot <- renderPlotly({
    #
    gene<-input$gene
    ss<-subset(s,s$Gene.ID %in% gene) %>% spread(measurement,value)
    ss$adj.p.value<-factor(ss$adj.p.value)
    ss$adj.p.value<-fct_recode(ss$adj.p.value,`Significant <=0.05`="1",`>0.05`="0")
    
    g2<-ggplot(ss,aes(x=Experiment,y=log2FoldChange,fill=adj.p.value,text=paste("padj:",padj)))+geom_bar(stat="identity")+
      theme_bw()+scale_fill_manual(values = c(`Significant <=0.05`=cols[1],`>0.05`=cols[2]))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
    ggplotly(g2)
  })
}

shinyApp(ui = ui, server = server)




