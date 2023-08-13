library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)

theme_set(theme_bw()+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size =16), plot.title = element_text(size = 20,hjust=0.5), strip.text = element_text(size=18), legend.text = element_text(size=10), legend.title = element_text(size=11,hjust=0.5)))

#### Main function: bivar_colorizer()
#### Arguments:
#### plotdata: a data.table, data.frame, or other object coercible into a data.table containing x and y axis variables and two continuous variables by which you want to co-colorize (keyvar1, keyvar2), below
#### xax: string specifying the x-axis column for plotting.
#### yax: string specifying the y-axis column for plotting.
#### keyvar1: string specifying the variable you would like to have be the legend's COLUMN-wise variable (left will be min value, right the max value)
#### keyvar2: string specifying the variable you would like to have be the legend's ROW-wise variable (bottom will be min value, top the max value)
#### keyvar.groups: string or vector of strings specifying grouping variables in plotdata for which to calculate the ranges of values when defining colors. (E.g., if you have measurements from 2 samples with different min-max ranges for the 2 variables being used to color code, you can get them on the same quantile-to-color scale by defining the column name with the sample id). If a vector of column names, the data will be subdivided in the order of those column names.
#### keyvar1.name.minmax: string specifying the legend label for the columnar legend colors, defining what the minimum value (left) and maximum value (right) signify. default for this could be "low-high".
#keyvar2.name.minmax: string specifying the legend label for the rowwise legend colors, defining what the minimum value (bottom) and maximum value (top) signify. default for this could be "low-high".
#### pseudocontinuous: TRUE: make an effectively continuous color gradient for the low to high range. FALSE: make a 4x4 color grid representing ~quartiles of the data (the points in the main plot will still be a fully continuous range of colors).
#### geomfun: which ggplot2 function to be used for the main plot: "geom_point","geom_jitter","geom_tile"
#### reds.keyvar: which coloring variable to color in the red spectrum. MUST BE one of keyvar1 or keyvar2. Defaults to be keyvar1 if not specified.
#### greens.keyvar which of the two coloring variables to place along the green spectrum, if any. MUST BE identical to keyvar1, keyvar2, OR set to NA if you want to use red-blue. defaults to keyvar2 if not specified.
#### blues.keyvar: which of the two coloring variables to place along the blue spectrum, if any. MUST BE identical to keyvar1, keyvar2, OR set to NA if you want to use red-blue. defaults to NA if not specified
#### custpltarea: default should work decently, but if needed, can input a custom list of 2 vectors of 4 numbers ranged 0:1 each, respectively for the main plot specifications and legend specifications in the shared plotting space. each vector specifies, in order: x (how far to thr right to begin the plot/legend, as a fraction of 1), y specifying how far up to begin plot legend as a fraction of 1, w specifying how much width of the plotting space should be dedicated to the plot/legend as a fraction of 1, and h, how much of the plot space height to use for the plot/legend.


### define helper functions modified from biscale (contbi_legbuild,palidate,palidate_names) to build the legend since biscale can't do dim>9 ; and some unmodified helper functions that are used along with.

bi_theme <- function(base_family = "sans", base_size = 24, bg_color = "#ffffff", font_color = "#000000", ...) {

    ggplot2::theme_minimal(base_family = base_family, base_size = base_size) +
        ggplot2::theme(

            # text defaults
            text = ggplot2::element_text(color = font_color),

            # remove all axes
            axis.line = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),

            # add a grid that blends into plot background
            panel.grid.major = ggplot2::element_line(color = bg_color),
            panel.grid.minor = ggplot2::element_blank(),

            # background colors
            plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
            panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
            legend.background = ggplot2::element_rect(fill = bg_color, color = NA),

            # borders and margins
            plot.margin = ggplot2::unit(c(.5, .5, .2, .5), "cm"),
            panel.border = ggplot2::element_blank(),
            panel.spacing = ggplot2::unit(c(-.1, 0.2, .2, 0.2), "cm"),

            # titles
            plot.title = ggplot2::element_text(size = ggplot2::rel(1.25), hjust = 0.5, color = font_color, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, color = font_color,
                                                  margin = ggplot2::margin(b = -0.1, t = -0.1, l = 2, unit = "cm"),
                                                  face = "bold", debug = FALSE),
            legend.title = ggplot2::element_text(color = font_color),
            legend.text = ggplot2::element_text(hjust = 0, color = font_color),

            # captions
            plot.caption = ggplot2::element_text(size = ggplot2::rel(.6), hjust = .5,
                                                 margin = ggplot2::margin(t = 0.2, b = 0, unit = "cm"),
                                                 color = font_color),
            ...
        )

}

bi_theme_legend <- function(base_family = "sans", base_size = 24, bg_color = "#ffffff", font_color = "#000000", ...) {

    ggplot2::theme_minimal(base_family = base_family, base_size = base_size) +
        ggplot2::theme(

            # text defaults
            text = ggplot2::element_text(color = font_color),

            # axes
            axis.text = ggplot2::element_text(size = ggplot2::rel(.8)),

            # add a grid that blends into plot background
            panel.grid.major = ggplot2::element_line(color = bg_color),
            panel.grid.minor = ggplot2::element_blank(),

            # background colors
            plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
            panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
            legend.background = ggplot2::element_rect(fill = bg_color, color = NA),

            # borders and margins
            plot.margin = ggplot2::unit(c(.5, .5, .2, .5), "cm"),
            panel.border = ggplot2::element_blank(),
            panel.spacing = ggplot2::unit(c(-.1, 0.2, .2, 0.2), "cm"),

            # titles
            plot.title = ggplot2::element_text(size = ggplot2::rel(1.25), hjust = 0.5, color = font_color, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, color = font_color,
                                                  margin = ggplot2::margin(b = -0.1, t = -0.1, l = 2, unit = "cm"),
                                                  face = "bold", debug = FALSE),
            legend.title = ggplot2::element_text(color = font_color),
            legend.text = ggplot2::element_text(hjust = 0, color = font_color),

            # captions
            plot.caption = ggplot2::element_text(size = ggplot2::rel(.6), hjust = .5,
                                                 margin = ggplot2::margin(t = 0.2, b = 0, unit = "cm"),
                                                 color = font_color),
            ...
        )

}

# validate palette input (unchanged from biscale)
palidate <- function(pal, dim, flip_axes, rotate_pal){

    # validate dim argument
    if (is.numeric(dim) == FALSE){
        stop("An integer scalar must be supplied for 'dim' that is greater than or equal to '2'.")
    }

    if (dim < 2 | (dim %% 1 == 0) == FALSE){
        stop("An integer scalar must be supplied for 'dim' that is greater than or equal to '2'.")
    }

    # validate logical arguments
    if (is.logical(flip_axes) == FALSE){
        stop("A logical scalar must be supplied for 'flip_axes'. Please provide either 'TRUE' or 'FALSE'.")
    }

    if (is.logical(rotate_pal) == FALSE){
        stop("A logical scalar must be supplied for 'rotate_pal'. Please provide either 'TRUE' or 'FALSE'.")
    }

    # validate palette
    if(length(pal) == 1){

        ## error if built-in palette supplied is not valid
        if (pal %in% c("DkViolet", "DkViolet2", "GrPink", "GrPink2", "DkBlue", "DkBlue2", "DkCyan", "DkCyan2", "Brown", "Brown2", "Bluegill", "BlueGold", "BlueOr", "BlueYl", "PinkGrn", "PurpleGrn", "PurpleOr") == FALSE){
            stop("The given palette is not one of the allowed options for bivariate mapping. Please see bi_pal's help file for a list of included palettes.")
        }

        ## error if legacy palette used with 4x4 plots
        if (dim == 4 & pal %in% c("DkViolet", "GrPink", "DkBlue", "DkCyan", "Brown")){
            stop(paste0("The legacy '", pal, "' palette does not support 4x4 bivarite mapping. Please use '", pal, "2' instead."))
        }

        ## error if dim is >4 with built-in palette
        if (dim > 4){
            stop("The palettes built-in to biscale only support bivariate maps where 'dim' is 2, 3, or 4.")
        }

    } else if (length(pal) > 1){

        ## ensure custom palette has correct number of values
        if (length(pal)/dim != dim){
            stop("The custom palette provided does not have the correct number of entries for the given dimensions.")
        }

        ## ensure formatting of hex values is correct
        bi_palidate_names(pal = pal, dim = dim)

        stopifnot(length(grep(pal,pattern="#"))==length(pal))
        stopifnot(sum(nchar(pal) %in% c(7,9))==length(pal))

        ## error for flipping and roating axes if dim > 4
        if (dim > 4 & flip_axes == TRUE){
            stop("Flipping axes for custom palettes is only available when 'dim' is 2, 3, or 4.")
        }

        if (dim > 4 & rotate_pal == TRUE){
            stop("Rotation for custom palettes is only available when 'dim' is 2, 3, or 4.")
        }

    }

}

bi_palidate_names <- function(pal, dim){

    x <- rep(x = 1:dim, times = dim)
    y <- sort(rep(x = 1:dim, times = dim))

    std <- paste(x, y, sep = "-")

    if (all(names(pal) == std) == FALSE){
        stop("Custom palette contains formatting errors - at least one entry name is incorrect.")
    }
    rm(x,y,std)

}

# flip palette axes
bi_pal_flip <- function(pal){
    if (length(pal) == 4){

        flipped <- pal

        flipped['1-2'] <- pal['2-1']
        flipped['2-1'] <- pal['1-2']

    } else if(length(pal) == 9){

        flipped <- pal

        flipped['1-2'] <- pal['2-1']
        flipped['1-3'] <- pal['3-1']
        flipped['2-1'] <- pal['1-2']
        flipped['2-3'] <- pal['3-2']
        flipped['3-1'] <- pal['1-3']
        flipped['3-2'] <- pal['2-3']

    } else if(length(pal) == 16){

        flipped <- pal

        flipped['1-2'] <- pal['2-1']
        flipped['1-3'] <- pal['3-1']
        flipped['1-4'] <- pal['4-1']

        flipped['2-4'] <- pal['4-2']
        flipped['3-4'] <- pal['4-3']

        flipped['4-2'] <- pal['2-4']
        flipped['4-3'] <- pal['3-4']

        flipped['3-2'] <- pal['2-3']
        flipped['2-3'] <- pal['3-2']

        flipped['2-1'] <- pal['1-2']
        flipped['3-1'] <- pal['1-3']
        flipped['4-1'] <- pal['1-4']

    }
    return(flipped)
}


# rotate axes
bi_pal_rotate <- function(pal){

    if (length(pal) == 4){

        rotated <- pal

        rotated['1-1'] <- pal['2-2']
        rotated['1-2'] <- pal['2-1']
        rotated['2-2'] <- pal['1-1']
        rotated['2-1'] <- pal['1-2']

    } else if (length(pal) == 9){

        rotated <- pal

        rotated['1-1'] <- pal['3-3']
        rotated['1-2'] <- pal['3-2']
        rotated['1-3'] <- pal['3-1']
        rotated['2-1'] <- pal['2-3']
        rotated['2-3'] <- pal['2-1']
        rotated['3-1'] <- pal['1-3']
        rotated['3-2'] <- pal['1-2']
        rotated['3-3'] <- pal['1-1']

    } else if (length(pal) == 16){

        rotated <- pal

        rotated['1-1'] <- pal['4-4']
        rotated['1-2'] <- pal['4-3']
        rotated['1-3'] <- pal['4-2']
        rotated['1-4'] <- pal['4-1']
        rotated['2-1'] <- pal['3-4']
        rotated['2-2'] <- pal['3-3']
        rotated['2-3'] <- pal['3-2']
        rotated['2-4'] <- pal['3-1']
        rotated['3-1'] <- pal['2-4']
        rotated['3-2'] <- pal['2-3']
        rotated['3-3'] <- pal['2-2']
        rotated['3-4'] <- pal['2-1']
        rotated['4-1'] <- pal['1-4']
        rotated['4-2'] <- pal['1-3']
        rotated['4-3'] <- pal['1-2']
        rotated['4-4'] <- pal['1-1']

    }
    return(rotated)

}

contbi_legbuild <- function(leg, dim, xlab, ylab, size, pad_width, pad_color, breaks, arrows, family){
    require(stringr)
    # global bindings
    # bi_fill = x = y = NULL

    # nse
    xQN <- as.name(xlab)
    yQN <- as.name(ylab)

    # create tibble for plotting
    leg <- as.data.table(data.frame(
        bi_class = names(leg),
        bi_fill = leg
    ))

    # get first classifier number
    leg$x <- as.integer(str_split_i(leg$bi_class, "-", 1))
    # " second " "
    leg$y <- as.integer(str_split_i(leg$bi_class, "-", 2))

    # create ggplot2 legend object
    ## initial build
    legend <- ggplot(data=leg,mapping = aes(x=x,y=y,fill=bi_fill))+geom_tile(lwd=pad_width,col=pad_color) +
        scale_fill_identity()

    ## optionally add breaks
    if (is.null(breaks) == FALSE){

        breaks_include <- TRUE

        if (length(breaks$bi_x) == dim){

            breaks_seq <- seq(from = 1, to = dim, by = 1)

        } else if (length(breaks$bi_x) == dim+1){

            breaks_seq <- seq(from = 0.5, to = dim+0.5, by = 1)

        }

        legend <- legend +
            ggplot2::scale_x_continuous(
                breaks = breaks_seq,
                labels = breaks$bi_x,
                expand = c(.015, .015)) +
            ggplot2::scale_y_continuous(
                breaks = breaks_seq,
                labels = breaks$bi_y,
                expand = c(.015, .015))

    }
    else {

        breaks_include <- FALSE

    }
    ## add arrows
    if (arrows == TRUE) {

        # add labels
        legend <- legend +
            ggplot2::labs(x = substitute(paste(xQN, ""%->%"")), y = substitute(paste(yQN, ""%->%"")))

    }
    else if (arrows == FALSE){

        # add labels
        legend <- legend +
            ggplot2::labs(x = xQN, y = yQN)

    }

    ## final legend elements
    legend <- legend +
        ggplot2::theme(axis.title = ggplot2::element_text(size = size)) +
        ggplot2::coord_fixed()

    ## add theme
    if (breaks_include == TRUE){
        legend <- legend + bi_theme_legend(base_size = size, base_family = family)
    } else if (breaks_include == FALSE){
        legend <- legend + bi_theme(base_size = size, base_family = family)
    }

    # return output
    return(legend)

}

### directly called function in plotting code belowq
contbi_legend <- function(pal, dim = 3, xlab, ylab, size = 10, flip_axes = FALSE, rotate_pal = FALSE, pad_width = NA, pad_color = '#ffffff',breaks = NULL, arrows = TRUE, base_family = "sans"){

    # global binding
    # bi_class = bi_fill = x = y = NULL

    # check parameters
    if (missing(pal) == TRUE){
        stop("A palette name or a custom palette vector must be specified for the 'pal' argument. Please see bi_pal's help file for a list of included palettes.")
    }

    if (is.logical(arrows) == FALSE){
        stop("A logical scalar must be supplied for 'arrows'. Please provide either 'TRUE' or 'FALSE'.")
    }

    if (missing(xlab) == TRUE){
        xlab <- "x var "
    }

    if (is.character(xlab) == FALSE){
        stop("The 'xlab' argument must be a character string.")
    }

    if (missing(ylab) == TRUE){
        ylab <- "y var "
    }

    if (is.character(ylab) == FALSE){
        stop("The 'ylab' argument must be a character string.")
    }

    if (is.numeric(size) == FALSE){
        stop("The 'size' argument must be a numeric value.")
    }

    # validate palette
    palidate(pal = pal, dim = dim, flip_axes = flip_axes, rotate_pal = rotate_pal)

    # create palette
    if (length(pal) == 1){
        legp <- bi_pal_pull(pal = pal, dim = dim, flip_axes = flip_axes, rotate_pal = rotate_pal)}
    else if (length(pal) > 1){
        legp <- pal
    }

    # build legend
    out <- contbi_legbuild(leg = legp, dim = dim, xlab = xlab, ylab = ylab, size = size, pad_width = pad_width, pad_color = pad_color, breaks = breaks, arrows = arrows, family = base_family)

    # return output
    return(out)

}




bivar_colorizer <- function(plotdata=d,xax="",yax="",keyvar1="",keyvar2="",keyvar.groups=NULL,keyvar1.name.minmax="",keyvar2.name.minmax="",pseudocontinuous=TRUE,geomfun=c("geom_point","geom_jitter","geom_tile"),reds.keyvar=keyvar1,greens.keyvar=keyvar2,blues.keyvar=NA,custpltarea=NA,return.color.appended.tab=T,return.ggcomponents=T,show.plot=T){
    library(data.table)
    library(ggplot2)
    library(cowplot)
    plotdata <- as.data.table(plotdata)

    # check for continuous variables to be colorized; check that only two RGB colors have been defined and that one of them is red (because blue + green is not going to be very useful)
    stopifnot(is.numeric(plotdata[,get(keyvar1)])&is.numeric(plotdata[,get(keyvar2)]))
    stopifnot(sum(is.na(c(blues.keyvar,reds.keyvar,greens.keyvar)))==1|is.na(reds.keyvar))

    # add columns converting color/fill key variables to fraction of their range, within one or more grouping variables if specified
    plotdata2 <- copy(plotdata)

    if (is.null(keyvar.groups)){
        plotdata2[,keyclr1:=((get(keyvar1))-min(get(keyvar1),na.rm=T))/(max(get(keyvar1),na.rm=T)-min(get(keyvar1),na.rm = T))]
        plotdata2[,keyclr2:=((get(keyvar2))-min(get(keyvar2),na.rm=T))/(max(get(keyvar2),na.rm=T)-min(get(keyvar2),na.rm = T))]
    }

    else{
        plotdata2[,keyclr1:=((get(keyvar1))-min(get(keyvar1),na.rm=T))/(max(get(keyvar1),na.rm=T)-min(get(keyvar1),na.rm = T)),by=keyvar.groups]
        plotdata2[,keyclr2:=((get(keyvar2))-min(get(keyvar2),na.rm=T))/(max(get(keyvar2),na.rm=T)-min(get(keyvar2),na.rm = T)),by=keyvar.groups]
    }

    ## generate a hexified RGB value using different channels from the proportions above in the user-specified colors. set the unused channel to zero.
    ## this is kinda messy. there's gotta be a better way to do it, but i'm not able to come up with it for now.
    if((reds.keyvar==keyvar1&greens.keyvar==keyvar2)==TRUE){
        plotdata2[,keyhex:=apply(.SD,MARGIN = 1,simplify = T,FUN=function(x){
            y <- as.data.table(x)
            ans <- rgb(red =y[1],blue=0,green=y[2],alpha=max(y[2],y[1]))
            ans}),.SDcols=c("keyclr1","keyclr2")]
    }
    else if ((reds.keyvar==keyvar2&greens.keyvar==keyvar1)==TRUE){
        plotdata2[,keyhex:=apply(.SD,MARGIN = 1,simplify = T,FUN=function(x){
            y <- as.data.table(x)
            ans <- rgb(red =y[2],blue=0,green=y[1],alpha=max(y[2],y[1]))
            ans}),.SDcols=c("keyclr1","keyclr2")]
    }
    else if (!is.na(blues.keyvar)&reds.keyvar==keyvar1&blues.keyvar==keyvar2){
        plotdata2[,keyhex:=apply(.SD,MARGIN = 1,simplify = T,FUN=function(x){
            y <- as.data.table(x)
            ans <- rgb(red =y[1],blue=y[2],green=0,alpha=max(y[2],y[1]))
            ans}),.SDcols=c("keyclr1","keyclr2")]
    }
    else if (!is.na(blues.keyvar)&reds.keyvar==keyvar2&blues.keyvar==keyvar1){
        plotdata2[,keyhex:=apply(.SD,MARGIN = 1,simplify = T,FUN=function(x){
            y <- as.data.table(x)
            ans <- rgb(red =y[2],blue=y[1],green=0,alpha=max(y[2],y[1]))
            ans}),.SDcols=c("keyclr1","keyclr2")]
    }
    else{
        stop("A variable to span the reds is not defined.")
    }

    # now the stuff for making the color legend, which ggplot2 has no functionality for when we have a two-variable two-gradient scheme like this.
    ## use biscale to make a 100x100 or 4x4 classification pairing, then unload it because we need to used modified internals of the package defined at the top of this source
    plotdata3 <- copy(plotdata2)
    coldim <- ifelse(pseudocontinuous,yes=100,no=4)
    plotdata3$keyvar1.bins <-
        cut(plotdata3[,keyclr1],breaks = seq(from=0,to=1,by=(1/coldim)),include.lowest = T)
    plotdata3$keyvar2.bins <-
        cut(plotdata3[,keyclr2],breaks=seq(from=0,to=1,by=(1/coldim)),include.lowest = T)

    ## make quantile-quantile classifications based on breaks and n(breaks)-1 (per biscale documentation) above. don't aactually load the package
    plotdata4 <- biscale::bi_class(plotdata3,x=keyvar1.bins,y=keyvar2.bins,dim = coldim)

    ### now the pain in the ass part--make a "discrete" classification scheme that spans like 200 breaks so we can get an effectively continuous legend from biscale (which only works in categorical terms/legend grids), but is the only package to do this in R.
    # the biscale classifications come out in terms of color legend space as the color legend's COLUMN, FOLLOWED BY the color legend's ROW (fucking confusingly as hell)
    i<-1 #i is the ROW index of the color legend grid, starting from BOTTOM row by default.
    custpal <- c()
    for(i in c(1:100)){
        j<-1 # j is the COLUMN index of the color legend grid, starting from left
        # therefore, bottom left="1-1", top left = "coldim-1",topright="coldim-coldim",bottom right = "1-coldim"
        # and therefore, we are assembling the legend in rows, left to right, from the bottom up
        for (j in c(1:100)){
            curbiclass <- paste0(j,"-",i)
            # make [min,min] white
            #if(curbiclass==""){stop("what")}
            if (i==1&j==1){
                custpal <- c("1-1"="#FFFFFF")
                next(j)
            }
            else if(nrow(plotdata4[bi_class==curbiclass])==0&j>1){
                # if the current class doesn't exist, just repeat the color of the last classification (since we are building left->right,bottom to top, this would be the column preceding (left of) column j), the value we last entered. if we're in the first column, we have to do this differently (next conditional).
                #(e.g., use class 1-88's color for 1-89 when no 1-89 classes exist), if such are present
                custpal <- c(custpal,
                             custpal[length(custpal)])
                names(custpal)[length(custpal)] <- curbiclass
                next(j)
            }
            else if (nrow(plotdata4[bi_class==curbiclass])==0&j==1&i>1){
                ## if we are in the first column, go one slot down since we can't go one slot left (because j=1).
                ## grab the last color class (e.g., if none of class 1-51, grab color of class 1-50; if no 1-50, 1-49-1,...)
                k <- i-1
                fetchbiclass <- paste0(j,"-",k)
                tmp <- custpal[fetchbiclass]
                custpal <- c(custpal,tmp)
                names(custpal)[length(custpal)] <- curbiclass
                rm(fetchbiclass,k)
                next(j)
            }
            else if (nrow(plotdata4[bi_class==curbiclass])==1){
                custpal <- c(custpal,
                             unique(plotdata4[bi_class==curbiclass,keyhex]))
                names(custpal)[length(custpal)] <- curbiclass
                next(j)
            }
            else{
                ## if multiple rows in a given class, get the most common color to put in the key (since our rgb-hex colors are going to be more continuous), or if a tie, get the first one for the max occurences
                newval <- plotdata4[bi_class==curbiclass,.N,by="keyhex"][N==max(N),keyhex][1]
                custpal <- c(custpal, newval)
                names(custpal)[length(custpal)] <- curbiclass
                next(j)
            }
        }
    }
    rm(i,j)

    # yuck, ok so now thats done and we have our corresponding palette to generate the 2D legend.

    # now we can finally make the plot. we need to define the coloring variable in aes() and then subsequently clarify that we want to use the identity of this variable as the color (scale_color_identity).

    # with the requested geom
    if (geomfun=="geom_point"){
        pltmain <- ggplot(plotdata4,aes(x=.data[[xax]],y=.data[[yax]],col=keyhex),show.legend=FALSE)+geom_point()+scale_color_identity(plotdata4$keyhex)
    }
    else if (geomfun=="geom_jitter"){
        pltmain <- ggplot(plotdata4,aes(x=.data[[xax]],y=.data[[yax]],col=keyhex),show.legend=FALSE)+geom_jitter()+scale_color_identity(plotdata4$keyhex)
    }
    else if (geomfun=="geom_tile"){
        pltmain <- ggplot(plotdata4,aes(x=.data[[xax]],y=.data[[yax]],col=keyhex),show.legend=FALSE)+geom_tile()+scale_color_identity(plotdata4$keyhex)
    }

    if(!is.null(keyvar.groups)){
        pltmain <- pltmain+facet_wrap(.~.data[[keyvar.groups]])
    }

    # build legend with biscale-derived code (the package's function doesn't work if dim > 9 because of a silly way of splitting the classification strings, so we write our own function for this at the top.)
    legn <- contbi_legend(pal=custpal,dim=coldim,xlab=keyvar1.name.minmax,ylab=keyvar2.name.minmax)

    # create final plot with cowplot and user specified dims if desired

    if(!is.list(custpltarea)){
        custpltarea <- list(c(0,0,0.6,1),c(0.625,0.5,0.365,0.5))
    }
    pltx <- c(custpltarea[[1]][1],custpltarea[[2]][1])
    plty <- c(custpltarea[[1]][2],custpltarea[[2]][2])
    pltw <- c(custpltarea[[1]][3],custpltarea[[2]][3])
    plth <- c(custpltarea[[1]][4],custpltarea[[2]][4])

    returnplot <- ggdraw()+draw_plot(pltmain,x=pltx[1],y=plty[1],width = pltw[1],height=plth[1])+draw_plot(legn,x=pltx[2],y=plty[2],width=pltw[2],plth[2])
    # returnplot <- pltmain

    if(show.plot==T){
        returnplot
    }
    if(return.color.appended.tab==T&return.ggcomponents==T){
        return(list(returnplot,pltmain,legn,plotdata4))
    }
    else if (return.color.appended.tab==T&return.ggcomponents==F){
        return(list(returnplot,plotdata4))
    }
    else if (return.color.appended.tab==F&return.ggcomponents==T){
        return(list(returnplot,pltmain,legn))
    }
    else{return(returnplot)}
}
