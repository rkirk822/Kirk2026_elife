

shinyUI(fluidPage(

    titlePanel("Mini Display & Analysis"),

    wellPanel(
      fluidRow(column(width=2,
        numericInput("group", "Group", value=1, min=1, max=20, step=1)),
        column(width=8,
          checkboxGroupInput("checkGroup", label = h3("mini #"), 
          choices = list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"=8,"9"=9,
              "10"=10), selected = 1, inline = T))),
      fluidRow(column(width=2, 
        radioButtons("rb", "trace type?", c("raw","filtered"))),
        column(width=2,
          radioButtons("gOrb", "?", c("Good","Bad"))),
        column(width=2,
          checkboxInput("bFit", "Show Fit?")),
        column(width=3, 
          actionButton("summ", "summarize")),
        column(width=3, 
          shinyDirButton("reload", "change cell", title="",buttonType = "default", class = NULL, icon = NULL, style = NULL)),
      )),

    fluidRow(column(width=8,
      plotlyOutput("mPlot"))),
    fluidRow(column(width=12, 
      DTOutput("mTable")))
))

