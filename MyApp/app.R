#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#install.packages("shiny")
#install.packages("shinythemes")
#library(shiny)
#library(shinythemes)

######### UI #########

ui <- fluidPage(theme = shinytheme("cosmo"), 
                titlePanel("Julia's Shiny App"),
                sidebarLayout(
                  sidebarPanel(
                    
                    numericInput(inputId = "b1", label = "Leeftijd", value = 0),
                    br(),
                    radioButtons("b2", label = "Hypertensie?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    radioButtons("ct1", label = "Aortaboog elongatie aanwezig?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    radioButtons("ct2", label = "Aortaboog type C aanwezig?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    radioButtons("ct3", label = ">90 hoek aanwezig in Inn of Cca?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    radioButtons("ct4", label = ">90 hoek aanwezig na bifurcatie?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    radioButtons("ct5", label = ">70% stenose aanwezig in de ICA?",
                                 choices = list("Ja" = 1, "Nee" = 0)),
                    br(),
                    br(),
                    actionButton("cal", label = "Bereken kans"),
                  ),
                  
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Voorspelling",
                               wellPanel(
                                 br(),
                                 h3(strong("Kans op falen:")),
                                    br(),
                                    br(),
                                 tableOutput("values")
                               )),
                    tabPanel("Aortaboog classificaties", fluidRow(
                      column(width = 6, wellPanel(p(
                        br(),
                        br(),
                        div(img(src = "variant.png", height = 400, width = 250), style = "text-align: center;"),
                        br(),
                        br(),
                        "Aortaboog elongatie: de truncus brachiocelphalicus ontstaat", strong("inferieur"), "ten opzichte van de aortaboog."))),
                      column(width = 6, wellPanel(p(
                        br(),
                        br(),
                        div(img(src = "variant.png", height = 400, width = 250), style = "text-align: center;"),
                        br(),
                        br(),
                        "Aortaboog variant: de", strong("linker carotide"), "ontstaat op de truncus brachiocephalicus.",
                      ))))),
                    tabPanel("Achtergrond informatie model", fluidRow( 
                      column(width = 12, wellPanel(
                        br(),
                        div(img(src = "mrclean.jpg", height = 100, width = 260), style = "text-align: center;"),
                        br(),
                        br(),
                        h4(strong("Model ontwikkelaars")),
                        "Het model dat gebruikt wordt voor het berekenen van de kans op het succesvol 
                        volbrengen van de transfemorale route is ontwikkeld door M.A.A. van Walderveen, G.Holswilder en J. van Hees.
                        ",
                        br(),
                        br(),
                        h4(strong("Stappen in het ontwikkelen van het model")),
                        "Ten grondslag aan het model ligt de", strong(a("MR CLEAN Registry."), href = "https://www.mrclean-trial.org/"), "Dit register bestaat uit meer dan
                        3000 herseninfarct patienten behandeld in 16 verschillende ziekenhuizen verspreid door heel Nederland.",
                        br(),
                        br(),
                        "Voor het ontwikkelen van het model is een deel van dit register gebruikt. Een ander
                        deel van dit register is gebruikt om the model te valideren (testen of het model betrouwbare voorspellingen kon doen).
                        "
                      ))))
                  )
                )
)
)


####### SERVER #######

server <- function(input, output, session) { 
 
  formule <- reactive({
    
    data.frame(Name = c("b1",
                        "b2", 
                        "ct1", 
                        "ct2", 
                        "ct3", 
                        "ct4", 
                        "ct5",
                        "metfactor"),
    Value <- as.numeric(c(input$b1, input$b2, input$ct1, input$ct2, input$ct3, input$ct4, input$ct5)))
    beta <- as.numeric(c(0.025781720, -0.6600788, 1.5103744, 1.0567201, 1.1750418, 0.9186653, 1.4339806))
    metfactor <- value %*% beta
    
  })
  
  output$values <- renderPrint({
    formule()
  })
}

###### RUN APP #######
shinyApp(ui = ui, server = server)


