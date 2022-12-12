library(data.table)
library(ggplot2)
library(ggridges)
library(scales)
library(gridExtra)

plot_count = function(dat = NULL, xlab = NULL) {
    ggplot(data = dat, aes(x = Coord, y = Max)) +
        geom_col() +
        theme_classic() +
        scale_x_continuous(limits = c(0, max(dat$Coord)), labels=comma, breaks = seq(from = 0, to = max(dat$Coord), by = 10000000)) +
        scale_y_continuous(labels = comma) +
        xlab(xlab) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              axis.line  = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0)))
}

dat1 = fread("input/grch38-queries-per-10k-window", header=T, nThread=4)
dat2 = fread("input/chm13v2-queries-per-10k-window", header=T, nThread=4)
dat3 = fread("input/chm13v2only-queries-per-10k-window", header=F, nThread=4)

colnames(dat1) = c("Coord", "avg", "num0", "Min", "Max")
colnames(dat2) = c("Coord", "avg", "num0", "Min", "Max")
colnames(dat3) = c("Coord", "avg", "num0", "Min", "Max")
head(dat1)
head(dat2)
head(dat3)

a=plot_count(dat1, xlab = "GRCh38Y")
b=plot_count(dat2, xlab = "T2T-CHM13v2.0Y")
c=plot_count(dat3, xlab = "T2T-CHM13v2.0Y Only")

data_plot <- arrangeGrob(a, b, c, nrow = 3)
data_plot
ggsave(file = "output/refseq_contam.pdf", data_plot, width = 100, height = 70, units = "mm", device=cairo_pdf)
ggsave(file = "output/refseq_contam.png", data_plot, width = 100, height = 70, units = "mm")


# Draw length distribution
dat1 = fread("input/grch38-hits-per-query.filt.hsat.category", header=F, nThread=4)
dat2 = fread("input/chm13v2-hits-per-query.filt.hsat.category", header=F, nThread=4)
dat3 = fread("input/chm13v2-hits-per-query.filt.hsat.category.chm13v2_only", header=F, nThread=4)

colnames(dat1) = c("RefSeqID", "Length", "Type")
head(dat1)
colnames(dat2) = c("RefSeqID", "Length", "Type")
head(dat2)
colnames(dat3) = c("RefSeqID", "Length", "Type")
head(dat3)

plot_length = function(dat = NULL, xlab = NULL, lim_y) {
    ggplot(data = dat, aes(x = Length, color = Type, fill = Type)) +
        geom_histogram(bins = 50, alpha=0.7, size=0.2) +
        theme_bw() +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        annotation_logticks(sides = "b") +
        scale_y_continuous(labels = comma, limits = c(0, lim_y)) +
        xlab(xlab) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              axis.line  = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              legend.title = element_blank(),
              legend.text = element_text(size = 5),
              legend.spacing = unit(0.5, 'mm'),
              legend.position = c(0.9, 0.7),
              legend.key.size = unit(2, "mm"),
              legend.margin=margin(c(-2, 0, 0, 0), unit="mm"))
}

p1 = plot_length(dat1, xlab = "GRCh38Y", lim_y = 630)
p1
p2 = plot_length(dat2, xlab = "T2T-CHM13v2.0Y", lim_y = 630)
p2
p3 = plot_length(dat3, xlab = "T2T-CHM13v2.0Y Only", lim_y = 160)
p3

data_plot <- arrangeGrob(p1, p2, p3, nrow = 3)
data_plot
ggsave(file = "output/refseq_contam_len.pdf", data_plot, width = 80, height = 70, units = "mm", device=cairo_pdf)
ggsave(file = "output/refseq_contam_len.png", data_plot, width = 80, height = 70, units = "mm")


ggplot(data = dat1, aes(x = Length, color = Type, fill = Type)) +
    geom_histogram(bins = 50, alpha=0.7, size=0.2) +
    geom_histogram(data=dat2, bins = 50, alpha=0.2, size=0.2) +
    theme_bw() +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    scale_y_continuous(labels = comma, limits = c(0, 750)) +
    theme(axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=5),
          axis.line  = element_line(size=0.1),
          axis.ticks = element_line(size=0.1),
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.spacing = unit(0.5, 'mm'),
          legend.position = c(0.9, 0.8),
          legend.key.size = unit(2, "mm"),
          legend.margin=margin(c(-2, 0, 0, 0), unit="mm"))
