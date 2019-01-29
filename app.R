# J. de Ligt
# Based on https://plot.ly/r/shinyapp-explore-diamonds/#code

library(shiny)
library(plotly)

# Human data
library(org.Hs.eg.db)

multmerge <- function(mypath){
  # Gather all RPKM files available
  filenames = list.files(path=mypath, full.names=TRUE, pattern="*_RPKM.txt")
  # Read files
  #, nrows=10
  datalist = lapply(filenames, function(x){read.table(file=x, header=T, sep='\t', row.names=1)})
  Reduce(function(x,y) {transform(merge(x,y,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)}, datalist)
}

rpkmdata <- data.frame(multmerge("counts/"))
metadata <- data.frame(read.table("RNAMetaData.txt", sep='\t', header=T, row.names=1))

genesymbols <- select(org.Hs.eg.db, rownames(rpkmdata), "SYMBOL", "ENSEMBL")
genes <- genesymbols[!duplicated(genesymbols$ENSEMBL),]
genesymbols <- NULL
#colnames(genes) <- c("value","label")

nms <- colnames(metadata)
metadata$Time <- as.factor(metadata$Time)
print(head(rpkmdata))
print(head(metadata))

#ENSG00000000003
#ENSG00000139292
ui <- fluidPage(

    headerPanel("Organoid RNAseq explorer"),
    sidebarPanel(
        selectInput('x', 'X', choices=nms, selected="Time"),
        selectizeInput('gene', 'Gene', choices=genes, selected="TSPAN6", multiple=F),
        selectInput('color', 'Color', choices=nms, selected="Tissue"),

        selectInput('facet_row', 'Facet Row', c(None='.', nms)),
        selectInput('facet_col', 'Facet Column', c(None='.', nms), selected="Tissue"),
        width=3
    ),
    mainPanel(
      plotlyOutput('boxPlot', height="900px")
    )
)

server <- function(input, output) {

  #add reactive data information. Dataset = built in diamonds data
  dataset <- reactive({

    id <- NULL
    if ( ! input$gene %in% rownames(rpkmdata)) {
      id <- which(genes$SYMBOL==input$gene)
    } else {
      id <- input$gene
    }

    dat <- data.frame(merge(t(rpkmdata[id,]), metadata, by=0))
    rownames(dat) <- dat$Row.names
    colnames(dat)[1] <- "Sample"
    colnames(dat)[2] <- "RPKM"

    return(dat)

  })

  output$boxPlot <- renderPlotly({

    # build graph with ggplot syntax
    p <- ggplot(dataset(), aes_string(x=input$x, y="RPKM", color=input$color)) +
       geom_boxplot(size=.75) + geom_point() + ggtitle(input$gene)

    # if at least one facet column/row is specified, add it
    facets <- paste(input$facet_row, '~', input$facet_col)
    if (facets != '. ~ .') p <- p + facet_grid(facets)

    # display with plotly
    ggplotly(p) %>%
      layout(height=600, autosize=T, width=800)

  })

}

shinyApp(ui, server)
