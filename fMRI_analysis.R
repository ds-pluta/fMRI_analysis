library(shiny)
library(R.matlab)

#----------------------------------------------------------------------------#
# App parameters

ts.default.dir = "~/UCI/ImagingGenetics/Data/fingerprint_210subs/time_series/cup/"
behav.default.dir = "~/UCI/ImagingGenetics/Data/fingerprint_210subs/fMRIbehav/cup/"
f.ts.stem = "SUB"
f.behav.stem = "fMRI_Cup_sub"
f.type = ".mat"
id.default = 2001
ROI.default = 1
trim.default = 20
hrf.p.default = c(12, 20, 3, 1, 6, 0, 32)
T.default = 16
RT.default = 2
trial.time = 580
feedback.time = 0.5


n.timepts = 290
n.rois = 375

#----------------------------------------------------------------------------#
# Functions

spm_Gpdf <- function(x, h, l){
    return(dgamma(x, h, l))
}

spm_hrf <- function(RT=2, p, T){

    # % Return a hemodynamic response function
    # % FORMAT [hrf,p] = spm_hrf(RT,p,T)
    # % RT   - scan repeat time
    # % p    - parameters of the response function (two Gamma functions)
    #     %
    # %                                                           defaults
    # %                                                          (seconds)
    # %        p(1) - delay of response (relative to onset)          6
    # %        p(2) - delay of undershoot (relative to onset)       16
    # %        p(3) - dispersion of response                         1
    # %        p(4) - dispersion of undershoot                       1
    # %        p(5) - ratio of response to undershoot                6
    # %        p(6) - onset (seconds)                                0
    # %        p(7) - length of kernel (seconds)                    32
    # %
    # % T    - microtime resolution [Default: 16]
    # %
    # % hrf  - hemodynamic response function
    # % p    - parameters of the response function

    dt <- RT/T
    u <- 0:ceiling(p[7]/dt) - p[6]/dt
    hrf <- spm_Gpdf(u,p[1]/p[3],dt/p[3]) - spm_Gpdf(u,p[2]/p[4],dt/p[4])/p[5]
    indices <- 0:floor(p[7]/RT)*T + 1
    hrf <- hrf[indices]
    hrf <- hrf/sum(hrf)
    return(hrf)
}

preprocess <- function(ts.dir=ts.default.dir, f.ts.stem, f.type, subj.id, trim.length, ROI.id,
                       filter.weights, filter.method="moving average"){
    #######################################################################
    # Subject data stored as matrix, COLS are ROIS, ROWS are time points  #
    #   subj.raw: original data                                           #
    #   subj.trim: original series trimmed by trim.length                 #
    #   subj.resids: linear trend removed from trimmed data               #
    #   subj: high-frequency noise removed, end result of preprocessing   #
    #######################################################################

    # Read one subject's data into matrix, all ROIs, all timepoints
    f.name <- paste0(ts.dir, f.ts.stem, subj.id, f.type)
    subj.raw <- readMat(f.name)[[1]][, ROI.id]
    subj.raw <- subj.raw - mean(subj.raw)

    # Trim leading time points
    subj.trim <- subj.raw[-(1:trim.length)]

    # Remove linear trend
    time.pts <- seq(2*(trim.length + 1), 2*n.timepts, 2)
    fit <- lm(subj.trim ~ time.pts)
    subj.detrend <- fit$residuals

    # Remove high-freq noise with small-window smoothing
    if(filter.method == "moving average")
        subj <- filter(subj.detrend, filter.weights, method="conv")
    subj <- subj[c(-1, -2, -269, -270)]

    return(list(raw=subj.raw, trimmed=subj.trim, detrended=subj.detrend, final=subj))
}

get.stimulus <- function(behav.dir=behav.default.dir, f.behav.stem, subj.id,
                         trim.length, use.reaction.time=FALSE, time.scale=1,
                         p=hrf.p.default, T=T.default, RT, trial.time=580,
                         stim.length=2.5) {
    f.behav <- paste0(behav.dir, f.behav.stem, subj.id, f.type)
    behav.dat <- readMat(f.behav)[[1]]
    stim.onset <- behav.dat[,11]
    if(use.reaction.time) {
        reaction.time <- round(behav.dat[,10], 2)
    }
    ts.stimulus <- rep(0, trial.time/time.scale)
    for(i in 1:stim.length/time.scale)
        ts.stimulus[stim.onset + i] <- 1
    return(ts.stimulus)
}

get.expected.BOLD <- function(ts.stimulus, hrf) {
    BOLD.active <- convolve(ts.stimulus, rev(hrf), type="o")
    BOLD.active <- BOLD.active - mean(BOLD.active)
    BOLD.rest <- convolve(1 - ts.stimulus, rev(hrf), type="o")
    BOLD.rest <- BOLD.rest - mean(BOLD.rest)
    return(list(BOLD.active=BOLD.active, BOLD.rest=BOLD.rest))
}

get.gaussian.weights <- function(lwr, upr, n.pts, stdev=1, show.plot=F){
    delta <- (upr - lwr)/(n.pts - 1)
    weights <- dnorm(seq(lwr, upr, delta), sd=stdev)
    weights <- weights/sum(weights)
    if(show.plot) plot(seq(lwr, upr, delta), weights)
    return(weights)
}

#----------------------------------------------------------------------------#


ui <- shinyUI(fluidPage(

    titlePanel("fMRI Data Analysis"),

    sidebarLayout(
        sidebarPanel(
            numericInput("subj.id", "Subject ID", id.default),
            numericInput("ROI.id", "ROI ID", ROI.default)
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Data Exploration",
                    fluidRow(
                        column(12,
                            plotOutput("main.plot")
                        )
                    ),
                    fluidRow(
                        tabsetPanel(
                            tabPanel("Preprocessing",
                                     numericInput("trim.length", "Trim Length", trim.default),
                                     radioButtons("smoothing", label="Smoothing method", choices=list("Moving average"=1, "Moving median"=2), selected=1),
                                     #numericInput("smoothing.window", "Smoothing window size (in secs)", value=5)
                                     textInput("filter.weights", "Filter Weights", "0.05, 0.1, 0.20, 0.3, 0.20, 0.1, 0.05")
                                     ),
                            tabPanel("HRF Parameters",
                                     fluidRow(column(3,
                                            numericInput("subj.scale", label="Subject fMRI Scale", value=0.4)
                                          ),
                                         column(3,
                                            sliderInput("p1",
                                                         "Delay of response:",
                                                         min = 1, max = 20, value = 6, step=0.5),
                                            sliderInput("p2",
                                                         "Delay of undershoot:",
                                                         min = 1, max = 30, value = 16, step=0.5),
                                            sliderInput("p3",
                                                         "Dispersion of response:",
                                                         min = 0, max = 5, value = 1, step=0.5),
                                            sliderInput("p4",
                                                         "Dispersion of undershoot:",
                                                         min = 0, max = 5, value = 1, step=0.5)),
                                          column(3,
                                            sliderInput("p5",
                                                         "Ratio of response to undershoot:",
                                                         min = 0, max = 15, value = 6, step=0.5),
                                            sliderInput("p6",
                                                         "Onset (seconds):",
                                                         min = 0, max = 20, value = 0, step=0.5),
                                            sliderInput("p7",
                                                         "Length of kernel:",
                                                         min = 1, max = 50, value = 32, step=0.5),
                                            sliderInput("T",
                                                         "Microtime resolution:",
                                                         min = 1, max = 20, value = 16, step=0.5))
                                          )),
                            tabPanel("Plot Options",
                                 fluidRow(
                                     column(6,
                                        checkboxInput("show.subj.BOLD", label="Show preprocessed BOLD signal", value=TRUE),
                                        checkboxInput("show.expected.BOLD", label="Show expected BOLD signal", value=FALSE),
                                        checkboxInput("show.stimulus", label="Show stimulus", value=FALSE),
                                        checkboxInput("show.fit", label="Show model fit", value=FALSE),
                                        checkboxInput("show.raw", label="Show raw data", value=FALSE),
                                        checkboxInput("show.trimmed", label="Show trimmed series", value=FALSE),
                                        checkboxInput("show.detrended", label="Show detrended signal", value=FALSE)
                                     ),
                                    column(6,
                                        sliderInput("y.lwr",
                                                    "y-axis lower limit:",
                                                    min = -100.0, max = 0, value = -50, step=0.1),
                                        sliderInput("y.upr",
                                                    "y-axis upper limit:",
                                                    min = 0, max = 100, value = 50, step=0.1),
                                        sliderInput("x.lwr",
                                                    "x-axis lower limit:",
                                                    min = 0, max = 270, value = 0, step=1),
                                        sliderInput("x.upr",
                                                    "x-axis upper limit:",
                                                    min = 0, max = 270, value = 270, step=1)
                    )))))),
                tabPanel("HRF Plot",
                    fluidRow(
                        column(12, "Plot")
                        )
                ),
                tabPanel("Data Directory",
                    fluidRow(
                        column(12,
                           textInput("ts.dir", label=h3("Time Series Directory"), value=ts.default.dir, width="100%"),
                           textInput("behav.dir", label="Behavioral Data Directory", value=behav.default.dir, width="100%")
))))))))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
    output$main.plot <- renderPlot({
        subj.dat <- preprocess(input$ts.dir, f.ts.stem, f.type, input$subj.id, input$trim.length, input$ROI.id, as.numeric(strsplit(input$filter.weights, ",")[[1]]))
        ts.stimulus <- get.stimulus(f.behav.stem=f.behav.stem, subj.id=2001, trim.length=20, RT=1)
        p <- c(input$p1, input$p2, input$p3, input$p4, input$p5, input$p6, input$p7)
        hrf <- spm_hrf(RT=1, p, input$T)

        plot.window(xlim=c(input$x.lwr, input$x.upr), ylim=c(input$y.lwr, input$y.upr))
        axis(1)
        axis(2)
        if(input$show.subj.BOLD)
            lines(subj.dat$final)
        if(input$show.detrended)
            lines(subj.dat$detrended, col="red", lty=2)
        if(input$show.raw)
            lines(subj.dat$raw, col="green")
        if(input$show.trimmed)
            lines(subj.dat$trimmed, col="blue", lty=3)
        if(input$show.stimulus)
            points(1:length(ts.stimulus)/2, ts.stimulus)
        if(input$show.expected.BOLD) {
            BOLD <- get.expected.BOLD(ts.stimulus, hrf)
            lines(BOLD[[1]])
        }
    })
})

# Run the application
shinyApp(ui = ui, server = server)
