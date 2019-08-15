library(shiny)
library(shinythemes)
library(DT)

fluidPage(
  titlePanel("VCF Annotations"),
  mainPanel(
    DT::dataTableOutput("VCF_Table")
  )
)