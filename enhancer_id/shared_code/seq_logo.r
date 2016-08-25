seqLogoGrid2 <- function (pwm, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE, xfontsize = 10, 
    yfontsize = 10, xmargin.scale = 1, ymargin.scale = 1, title = "", 
    titlefontsize = 15) 
{
    require(PWMEnrich)
    if (class(pwm) == "pwm") {
        pwm <- pwm@pwm
    }
    else if (class(pwm) == "data.frame") {
        pwm <- as.matrix(pwm)
    }
    else if (class(pwm) != "matrix") {
        stop("pwm must be of class matrix or data.frame")
    }
    if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) 
        stop("Columns of PWM must add up to 1.0")
    chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(pwm)
    if (ic.scale) {
        ylim <- 2
        ylab <- "Information content"
        facs <- PWMEnrich:::pwm2ic(pwm)
    }
    else {
        ylim <- 1
        ylab <- "Probability"
        facs <- rep(1, npos)
    }
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        column <- pwm[, j]
        hts <- 0.95 * column * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0) 
                letters <- PWMEnrich:::addLetter(letters, letter, x.pos, 
                  y.pos, ht, wt)
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }
    bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 2) * ymargin.scale
    leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 2) * xmargin.scale
    rightMargin = 2 * xmargin.scale
    topMargin = 2 * ymargin.scale
    pushViewport(plotViewport(c(bottomMargin, leftMargin, topMargin, 
        rightMargin)))
    pushViewport(dataViewport(0:ncol(pwm), 0:ylim, name = "vp1"))
    grid.polygon(x = unit(letters$x, "native"), y = unit(letters$y, 
        "native"), id = letters$id, gp = gpar(fill = letters$fill, 
        col = "transparent"))
    grid.text(title, y = 1.1, vjust = 1, gp = gpar(fontsize = titlefontsize))
    if (xaxis) {
        grid.xaxis(at = seq(0.5, ncol(pwm) - 0.5), label = 1:ncol(pwm), 
            gp = gpar(fontsize = xfontsize))
        grid.text("Position", y = unit(-3, "lines"), gp = gpar(fontsize = xfontsize))
    }
    if (yaxis) {
        grid.yaxis(gp = gpar(fontsize = yfontsize))
        grid.text(ylab, x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = yfontsize))
    }
    popViewport()
    popViewport()
}
