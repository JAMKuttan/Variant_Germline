# Load Packages
library(shiny)
library(DT)

function(input,output) {
  output$VCF_Table <- DT::renderDataTable(
    DT::datatable(Final)
  )
}