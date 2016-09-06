#----------------------------------------------------------------------------#
# App parameters

ts.default.dir = ""
behav.default.dir = ""

f.ts.stem = "SUB"
f.behav.stem = "fMRI_Cup_sub"
f.type = ".mat"
ROI.default = 1
trim.default = 20
p.default = c(12, 20, 3, 1, 6, 0, 32)
T.default = 16
RT.default = 2
trial.time = 580
feedback.time = 0.5

#subj.list <- read.delim("210sublist.txt", header=F)[,1]
subj.list = c("2001", "2002")
n.subjs = 2
n.timepts = 290
n.ROIs = 375

col.width = 5
#----------------------------------------------------------------------------#

ui <- shinyUI(fixedPage(

    titlePanel("fMRI Data Analysis"),

    wellPanel("This app allows for the visualization and exploration of a sample of the
              study fMRI data.  The observed blood oxygend level dependent (BOLD) responses
              for all 375 brain regions of interest (ROIs) for two subjects
              can be examined.  Different smoothing windows (using a Gaussian kernel smoother) and
              leading trim lengths can be selected.  The stimulus series for the
              Cups Task can also be displayed.  The most recent version of the
              code for this app is available on GitHub at repo ",
              tags$a("ds-pluta/fMRI_analysis", href="https://github.com/ds-pluta/fMRI_analysis"), "."),
    fluidRow(
        column(12,
               plotOutput("main.plot")
        )
    ),
    fixedRow(
        column(col.width,
               wellPanel(
                   h3("Subject ID and ROI"),
                   numericInput("subj.index", "Subject Index", value=1, min=1, max=n.subjs, width='25%'),
                   tags$br(),
                   tags$b("Subject ID"),
                   tags$br(),
                   textOutput("subj.id", inline=T),
                   tags$br(),
                   tags$hr(),
                   numericInput("ROI.id", "ROI ID", ROI.default, value=1, min=1, max=n.ROIs, width='25%'),
                   tags$br()
               )
        ),
        column(col.width,
               wellPanel(
                   h3("Plot Window"),
                   numericInput("x.lwr", "x-axis lower limit:",
                                min = 0, max = 270, value = 0, step=1),
                   numericInput("x.upr", "x-axis upper limit:",
                                min = 0, max = 270, value = 280, step=1),
                   numericInput("y.lwr", "y-axis lower limit:",
                                min = -2, max = 0, value = -0.3, step=0.05),
                   numericInput("y.upr", "y-axis upper limit:",
                                min = 0, max = 2, value = 0.3, step=0.05)
               )
        ),
        column(col.width,
               wellPanel(
                   h3("Plot Parameters"),
                   checkboxInput("show.subj.BOLD", label="Show preprocessed BOLD signal", value=TRUE),
                   checkboxInput("show.detrended", label="Show detrended signal", value=FALSE),
                   checkboxInput("show.stimulus", label="Show stimulus", value=FALSE),
                   checkboxInput("show.expected.BOLD", label="Show expected BOLD signal", value=FALSE)
               )
        ),
        column(col.width,
               wellPanel(
                   h3("Preprocessing"),
                   numericInput("smoothing.window", "Smoothing window size", value=5, min=3, max=33, step=2),
                   numericInput("trim.length", "Time Series Trim Length", value=20, min=0, max=100, step=2)
               )
        )
    )
))
