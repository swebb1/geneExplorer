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
#library(Biostrings)
library(stringr)

a<-readRDS("annotation.rds") ###Sort out factors
r<-readRDS("rna_data.rds") ### sort out missing data and merged data
pr<-readRDS("prot.rds")
p<-readRDS("pep.rds")
##Move to compile and create RDS
fa = readAAStringSet("../data/Human_Ref_Proteome_REVIEWED_2018-06-01.fasta")
names(fa)<-str_split(names(fa),"\\|",n = 3,simplify = T)[,2]

cols<-c(brewer.pal(n = 9,name = "Set1")[3],"grey39")
cols2<-c(brewer.pal(n = 9,name = "Set1")[2],"grey39")

ui<-shinyUI(dashboardPage(
  title="Data explorer",
  skin="blue",
  dashboardHeader(title = "Data Explorer", titleWidth = 400),
  dashboardSidebar(width=220,
                   sidebarMenu(
                     #sidebarSearchForm(textId = "searchText", buttonId = "searchButton",label = "Search..."),
                     menuItem("Home",tabName="home",icon=shiny::icon("home")),
                     menuItem("Multiple Gene Analysis",tabName="multi",icon=shiny::icon("chart-bar"))#,
                     #menuItem("Data Table",tabName="data",icon=shiny::icon("database"))
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
                  div(style = 'overflow-x: scroll', dataTableOutput('pepTable'))
                ),
                box(
                  title="Selected Protein",width = 3,status="primary",solidHeader=F,
                  div(style = 'overflow-x: scroll', dataTableOutput('plinksTable'))
                ),
                box(
                  title="Proteomics: Peptide Plot",width = 6,status="primary",solidHeader=F,
                  #uiOutput("protSelect"),
                  #uiOutput("pepSelect"),
                  textOutput("pepText"),
                  plotOutput("pepPlot")
                ),
                box(
                  title="Proteomics: Protein Groups Plot",width = 6,status="primary",solidHeader=F,
                  #uiOutput("protSelect"),
                  #uiOutput("pepSelect"),
                  #textOutput("pepText"),
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
  
  proteins<-reactive({
    gene<-genes()
    ps<-subset(pr,pr$Gene.ID %in% gene$Gene.ID) %>% spread(measurement,value)
    return(ps)
  })
  
  pep<-reactive({
    peptides<-peptides()
    if(dim(peptides)[1]==0){
      return(NULL)
    }
    if(length(input$pepTable_rows_selected)==0){
      peps<-peptides[input$anno_rows_all[1],]
    }
    else{
      peps<-peptides[input$pepTable_rows_selected[1],]
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
  
   seqloc<-reactive({
     #Get peptide sequence and ID
     pep<-pep()
     pepseq<-pep$Sequence
     id<-pep$Seq.ID

     ##Cutout alterations
     seq<-str_replace_all(pepseq,c("_"="","\\([a-o,q-z]+\\)"=""))
     seq<-str_replace_all(seq,"\\([a-z]+\\)","p")
     seqp<-str_replace_all(seq,"p","")
     ##Get protein sequence ##!! Loop for all isoforms?!!!
     fas<-fa[id]
     ##Search for sequence within protein and get starting location
     vm<-matchPattern(seqp,fas[[1]],max.mismatch = 3,with.indels = T)
     start<-start(vm@ranges)
     ##Get locations of phosphorylation within sequence and protein sequence
     loc<-as.data.frame(str_locate_all(seq,"p")[1])$start
     rloc<-loc+(start-1)
     rloc
   })
  
  output$pepTable <- DT::renderDT({
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
  
  output$pepText<-renderText({
    gene<-genes()
    pep<-pep()
    if(is.null(pep)){
      return(as.character(gene$Gene.Symbol))
    }
    ploc<-seqloc()    
    text<-paste(c(as.character(gene$Gene.Symbol),pep$Uniprot.ID,pep$Sequence,ploc),collapse=" ")
    return(text)
  })
      
  output$pepPlot <- renderPlot({
    sig2=input$sig2
    values=cols
    names(values)<-c(paste0(" Significant <= ",sig2),paste0(" > ",sig2))
    
    ps<-peptides()
    if(dim(ps)[1]==0){
      return(NULL)
    }
    pep<-pep()
    pss<-subset(ps,ps$Uniprot.ID==unique(pep$Uniprot.ID) & ps$ID==unique(pep$ID))
    
    g1<-ggplot(pss %>% mutate(fill = ifelse(pval<=sig2,paste0(" Significant <= ",sig2),paste0(" > ",sig2))),
               aes(text=paste("pval:",pval)))+
      geom_col(aes(x=Experiment,y=FoldChange,fill=fill))+
      scale_fill_manual(values = values)+
      theme_bw()+
      theme(text = element_text(size=20),legend.title = element_blank())+
      facet_grid(ID~Uniprot.ID)+#,scales="free",space="free")+
      coord_flip()
    g1
  })
  
  output$protPlot <- renderPlot({
    sig2=input$sig2
    values=cols
    names(values)<-c(paste0(" Significant <= ",sig2),paste0(" > ",sig2))
    
    ps<-proteins()
    if(dim(ps)[1]==0){
      return(NULL)
    }
    prot.name=unique(ps$Uniprot.ID)
    g1<-ggplot(ps %>% mutate(fill = ifelse(pval<=sig2,paste0(" Significant <= ",sig2),paste0(" > ",sig2))),
               aes(text=paste("pval:",pval)))+
      geom_col(aes(x=Experiment,y=FoldChange,fill=fill))+
      scale_fill_manual(values = values)+
      theme_bw()+
      ggtitle(prot.name)+
      theme(text = element_text(size=20),legend.title = element_blank())+
      facet_grid(~Uniprot.ID)+#,scales="free",space="free")+
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




