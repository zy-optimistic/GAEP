#! /usr/bin/Rscript
#library(optparse)

#Rscript snp_coverage_dotplot.R <input> <output_prefix> <png|pdf> options
#-w width
#-h height
#-c max coverage
#-s max SNV density

plot_scatter_with_density <- function( data, x_colname, y_colname, mean_cov, x_lim, y_lim, xbins, ybins ){
	library(ggplot2)
	library(ggExtra)
	commonTheme = list(labs(color="Density",fill="Density",
                            x="Coverage",
                            y="SNV density (%)"),
                            theme_bw(),
                            theme(legend.position=c(0,1),
                            legend.justification=c(0,1)
							)
						)
	p <- ggplot(data, aes(x=x_colname, y=y_colname)) +
	     xlim(0, x_lim) +
	     ylim(0, y_lim) +
		 geom_vline(aes(xintercept=mean_cov), colour="#BB0000", linetype="dashed") +
	     #geom_density_2d()
	     #stat_density_2d(aes(fill = ..level..), geom = "polygon")
	     geom_point(alpha = 0.1, na.rm=TRUE) +
	     stat_density_2d(aes(colour = ..level..), size=1.5, na.rm=TRUE) + 
	     scale_colour_gradient(low="green",high="red")  +
	     scale_linetype_manual(values = 5) +
	     commonTheme
	p1 <- ggMarginal(p, type="histogram",xparams = list(bins=xbins),yparams = list(bins=ybins))
	return(p1)
}

help_info <- function(){
	write("Usage: Rscript snp_coverage_dotplot.R <input> <output_prefix> <png|pdf> [options].\n", stderr())
	write("Options:", stderr())
	write("    -w width of plot", stderr())
	write("    -h height of plot", stderr())
	write("    -c max coverage", stderr())
	write("    -s max SNV density", stderr())
	return(1)
}


arg <- commandArgs(T)
pwidth  <- 0
pheight <- 0
cov_max <- 0
snv_max <- 0
input <- character()

i <- 1
while (i <= length(arg)) {
	if ( substr(arg[i], 1, 1) == "-" ) {
		if (substr(arg[i], 2, nchar(arg[i])) == "w") {
			if (length(grep("^[1-9]+$", arg[i+1], perl=TRUE))) {
				pwidth <- as.numeric(arg[i+1])
			}else {
				stop("Use non-numeric argument: ", arg[i+1])
			}
			i <- i + 2
		}else if (substr(arg[i], 2, nchar(arg[i])) == "h") {
			if (length(grep("^[1-9]+$", arg[i+1], perl=TRUE))) {
				pheight <- as.numeric(arg[i+1])
			}else {
				stop("Use non-numeric argument: ", arg[i+1])
			}
			i <- i + 2
		}else if (substr(arg[i], 2, nchar(arg[i])) == "c") {
			if (length(grep("^[1-9]+$", arg[i+1], perl=TRUE))) {
				cov_max <- as.numeric(arg[i+1])
			}else {
				stop("Use non-numeric argument: ", arg[i+1])
			}
			i <- i + 2
		}else if (substr(arg[i], 2, nchar(arg[i])) == "s") {
			if (length(grep("^[1-9]+$", arg[i+1], perl=TRUE))) {
				snv_max <- as.numeric(arg[i+1])
			}else {
				stop("Use non-numeric argument: ", arg[i+1])
			}
			i <- i + 2
		}else if (substr(arg[i], 2, nchar(arg[i])) == "help") {
			help_info()
			stop("Print help information.")
		}else {
			stop("Unknown option: ", substr(arg[i], 2, nchar(arg[i])))
		}
	}else {
		input<-append(input, arg[i])
		i <- i + 1
	}
}

if (length(input) < 3) {
	help_info()
	stop("Invalid arguments.")
}


data <- read.delim(input[1],header = FALSE)
if (cov_max == 0) {
	cov_max <- as.integer(mean(data$V4) * 2.5)
}
if (snv_max == 0) {
	snv_max <- mean(data$V5) * 4
}
if (snv_max < 1) {
    snv_max <- 1
}


#data <- data[data[,4]<cov_max,]
#data <- data[data[,5]<snv_max,]

if ( input[3] == "png" ) {
	out <- paste(input[2], ".png", sep = "")
	if (pwidth == 0) {
		pwidth <- 1000
	}
	if (pheight == 0) {
		pheight <- 1000
	}
	png(file = out , width=pwidth, height=pheight)
}else if ( input[3] == "pdf" ) {
	out <- paste(input[2], ".pdf", sep = "")
	if (pwidth == 0) {
		pwidth <- 20
	}
	if (pheight == 0) {
		pheight <- 20
	}
	pdf(file = out, width=pwidth, height=pheight)
}

print(cov_max)
print(snv_max)
cat(mean(data$V4), cov_max ,snv_max, "\n")
plot_scatter_with_density(data, data$V4 ,data$V5 , mean(data$V4), cov_max ,snv_max ,100 ,1000)

dev.off()
