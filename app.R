# J. de Ligt
# Based on https://plot.ly/r/shinyapp-explore-diamonds/#code

library(shiny)
library(plotly)

# Human data 
library("org.Hs.eg.db")

multmerge <- function(mypath){
  # Gather all RPKM files available
  filenames = list.files(path=mypath, full.names=TRUE, pattern="*_RPKM.txt")
  # Read files
  datalist = lapply(filenames, function(x){read.table(file=x, header=T, sep='\t', row.names=1)})
  Reduce(function(x,y) {merge(x,y)}, datalist)
}
  
rpkmdata <- data.frame(multmerge("counts/"))
metadata <- data.frame(read.table("RNAMetaData.txt", sep='\t', header=T, row.names=1))

genesymbols <- select(org.Hs.eg.db, rownames(rpkmdata), "SYMBOL", "ENSEMBL")
genes <- genesymbols[!duplicated(genesymbols$ENSEMBL),]

nms <- colnames(metadata)


ui <- fluidPage(

    headerPanel("Organoid RNAseq explorer"),
    sidebarPanel(
        selectInput('x', 'X', choices = nms, selected = "Tissue"),
        selectInput('gene', 'Gene', choices = genes$ENSEMBL, selected = "ENSG00000139292"),
        selectInput('color', 'Color', choices = nms, selected = "Tissue"),

        selectInput('facet_row', 'Facet Row', c(None = '.', nms), selected = "Tissue"),
        selectInput('facet_col', 'Facet Column', c(None = '.', nms)),
    ),
    mainPanel(
      plotlyOutput('boxPlot', height = "900px")
    )
)

server <- function(input, output) {

  #add reactive data information. Dataset = built in diamonds data
  dataset <- reactive({
    rpkmdata[input$gene,]
  })

  output$boxPlot <- renderPlotly({

    # build graph with ggplot syntax
    p <- ggplot(dataset(), aes_string(x = input$x, y=RPKM, color=input$color)) + 
       geom_boxplot() + geom_point()

    # if at least one facet column/row is specified, add it
    facets <- paste(input$facet_row, '~', input$facet_col)
    if (facets != '. ~ .') p <- p + facet_grid(facets)

    ggplotly(p) %>% 
      layout(height=800, autosize=TRUE)

  })

}

shinyApp(ui, server)
