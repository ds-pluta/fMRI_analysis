library(shiny)
library(R.matlab)

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
                       filter.weights){
    #######################################################################
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
    subj <- filter(subj.detrend, filter.weights, method="conv")

    return(list(raw=subj.raw, trimmed=subj.trim, detrended=subj.detrend, final=subj,
                detrended.scaled=subj.detrend/var(subj, na.rm=T), final.scaled=subj/var(subj, na.rm=T)))
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

get.gaussian.weights <- function(n.pts, stdev=1, show.plot=F){
    delta <- 4/(n.pts - 1)
    lwr <- -2
    upr <- 2
    weights <- dnorm(seq(lwr, upr, delta), sd=stdev)
    weights <- weights/sum(weights)
    if(show.plot) plot(seq(lwr, upr, delta), weights)
    return(weights)
}

inputRow <- function (inputId, label, value="", type="text", ...)
{
    div(style="display:inline-block",
        tags$label(label, `for` = inputId),
        tags$input(id=inputId, type=type, value = value,...))
}

#----------------------------------------------------------------------------#

server <- shinyServer(function(input, output) {
    output$subj.id <- renderPrint({subj.list[input$subj.index]})
    output$main.plot <- renderPlot({
        subj.dat <- preprocess(ts.default.dir, f.ts.stem, f.type, subj.list[input$subj.index],
                               trim.length=20, input$ROI.id,
                               get.gaussian.weights(n.pts=input$smoothing.window, stdev=1))
        ts.stimulus <- get.stimulus(f.behav.stem=f.behav.stem, subj.id=2001, trim.length=20, RT=1)
        p <- c(input$p1, input$p2, input$p3, input$p4, input$p5, input$p6, input$p7)
        hrf <- spm_hrf(RT=1, p.default, T.default)

        plot.window(xlim=c(input$x.lwr, input$x.upr), ylim=c(input$y.lwr, input$y.upr), xlab="Time Pt (2s)", ylab="Signal Intensity (Scaled)", xaxs="r", yaxs="r")
        axis(1)
        axis(2)
        legend("topright", legend=c("Observed BOLD Signal"), col="black", lty=1, lwd=2)
        if(input$show.subj.BOLD)
            lines(subj.dat$final.scaled, xlab="Time Pt (2s)", ylab="Signal Intensity (Scaled)", lwd=2, col=rgb(0.1, 0.1, 0.1, 0.8))
        if(input$show.detrended) {
            lines(subj.dat$detrended.scaled, col=rgb(0.7, 0.2, 0.2, 0.4), lwd=1.5, lty=1, xlab="Time Pt (2s)", ylab="Signal Intensity (Scaled)")
            #legend("topright", legend=c("Observed BOLD Signal", "Expected BOLD Signal", "Detrended Observed Signal"), col=c("black", rgb(0, 0.4, 0.8, 0.8), rgb(0.7, 0.2, 0.2, 0.4)), lty=c(1,2, 1), lwd=2)
        }
        if(input$show.stimulus)
            points(1:(length(ts.stimulus)/2 - input$trim.length), 0.9*input$y.upr*ts.stimulus[seq(2*input$trim.length + 1, length(ts.stimulus), 2)], pch=73)
        if(input$show.expected.BOLD) {
            BOLD <- get.expected.BOLD(ts.stimulus, hrf)
            lines(1:(length(ts.stimulus)/2 - input$trim.length), BOLD[[1]][seq(2*input$trim.length + 1, length(ts.stimulus), 2)],
                  col=rgb(0, 0.4, 0.8, 0.8), lty=2, lwd=2)
            legend("topright", legend=c("Observed BOLD Signal", "Expected BOLD Signal"), col=c("black", rgb(0, 0.4, 0.8, 0.8)), lty=c(1,2), lwd=2)
        }
    })
})
