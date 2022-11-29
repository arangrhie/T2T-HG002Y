library(data.table)
library(ggplot2)
library(ggridges)
library(scales)
library(gridExtra)
library(gtable)

frac_min = -3
frac_max = 103
dp_min = -10
#dp_c_max = 800
af_c_max = 500
P_NonRef_Lab = "Alternate Allele (%)"
ref_color=c("#E1812C", "#A94359")

fancy_scientific <- function(d) {
    # turn in to character string in scientific notation
    d <- format(d, scientific = TRUE)
    # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format
    d <- gsub("^(.*)e\\+", "'\\1'%*%10^", d)
    # convert 0x10^00 to 0
    d <- gsub("\\'0[\\.0]*\\'(.*)", "'0'", d)
    # return this as an expression
    parse(text=d)
}

plot_count = function(dat = NULL, y_label = "Alternate Allele (%)", xmax = NULL, median = NULL) {
    ggplot(data = dat, aes(x = DP, y = Alt, color = Reference, fill = Reference)) +
        geom_vline(xintercept = median, size = 0.1, linetype = "dashed", color = "darkgray") +
        geom_point(shape=21, alpha=0.4, size = 0.1) +
        scale_color_manual(values = ref_color) +
        scale_fill_manual(values = ref_color) +
        theme_classic() +
        xlab("Total Read Depth") +
        ylab(y_label) +
        scale_x_continuous(labels = comma, limits = c(dp_min, xmax)) +
        scale_y_continuous(labels = comma, limits = c(frac_min, frac_max)) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              axis.line  = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0)))
}

plot_hist_dp = function(dat = NULL, xmax) {
    ggplot(data = dat, aes(x = DP, color = Reference, fill = Reference)) +
        geom_histogram(binwidth = 5, alpha = 0.4, size = 0.1, position="identity", na.rm = TRUE) +
        scale_color_manual(values = ref_color) +
        scale_fill_manual(values = ref_color) +
        theme_classic() +
        ylab("Count") +
        scale_x_continuous(labels = comma, limits = c(dp_min, xmax)) +
        scale_y_continuous(labels = comma) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.title.x = element_blank(),
              axis.text  = element_text(size=5),
              axis.text.x  = element_blank(),
              axis.line  = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              legend.position="none")
}

plot_hist_p_alt = function(dat = NULL) {
    ggplot(data = dat, aes(x = Alt, color = Reference, fill = Reference)) +
        geom_histogram(binwidth = 2, alpha = 0.4, size = 0.1, position="identity", na.rm = TRUE) +
        scale_color_manual(values = ref_color) +
        scale_fill_manual(values = ref_color) +
        theme_classic() +
        ylab("Count") +
        scale_x_continuous(labels = comma, limits = c(frac_min, frac_max)) +
        scale_y_continuous(labels = comma) +
        coord_flip() +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.title.y = element_blank(),
              axis.text  = element_text(size=5),
              axis.text.y  = element_blank(),
              axis.line  = element_line(size=0.1),
              axis.ticks = element_line(size=0.1),
              legend.position="none")
}


plot_all = function(dat = NULL, outfile = NULL, median = NULL) {
    dp_max = max(s$DP) + 30
    a = plot_count(dat, y_label = P_NonRef_Lab, xmax = dp_max, median = median)
    b = plot_hist_dp(dat, xmax = dp_max)
    c = plot_hist_p_alt(dat)
    l = gtable_filter(ggplotGrob(a), "guide-box")
    data_plot <- arrangeGrob(b, l, a + theme(legend.position="none"), c, nrow = 2, ncol = 2, widths = c(4, 1), heights = c(1, 4))
    
    ggsave(file = paste("output/", outfile, ".png", sep = ""), data_plot, width = 90, height = 90, units = "mm")
    ggsave(file = paste("output/", outfile, ".pdf", sep = ""), data_plot, width = 55, height = 55, units = "mm", device=cairo_pdf)
}


dat1 = fread("input/v1_GRCh38Y_variants.syntenic.noPAR.txt", header=T, nThread=4)
dat2 = fread("input/v2_HG002Y_variants.syntenic.noPAR.txt", header=T, nThread=4)

head(dat1)
max(dat1$HG00116_DP, dat1$HG01130_DP, dat1$HG01885_DP)
# 648

min(dat1$HG00116_DP, dat1$HG01130_DP, dat1$HG01885_DP)

# HG00116 (R1b - should match better GRCh38Y)
dat1_s = dat1[dat1$HG00116_DP>dat1$HG00116_AD_REF,]
dp = dat1_s$HG00116_DP
P_NonRef = 100*(1 - dat1_s$HG00116_AD_REF/dp)
s1 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("GRCh38-Y"))
head(s1)

dat2_s = dat2[dat2$HG00116_DP>dat2$HG00116_AD_REF,]
dp = dat2_s$HG00116_DP
P_NonRef = 100*(1 - dat2_s$HG00116_AD_REF/dp)
s2 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("T2T-Y"))
sample1 = rbind(s1, s2)

summary(s1$DP)[[3]] # Median
#  Min. 1st Qu.  Median   Mean  3rd Qu.   Max. 
# 1.00   14.00   28.00   37.06   48.50  648.00 
summary(s2$DP)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   16.00   20.00   20.52   24.00   75.00

summary(s1$Alt)
summary(s1[s1$Alt > 90,])

plot_all(sample1, "Variants_HG00116_R1b", median = 20)

# HG01130 (J1 - should match better HG002Y)
dat1_s = dat1[dat1$HG01130_DP>dat1$HG01130_AD_REF,]
dp = dat1_s$HG01130_DP
P_NonRef = 100*(1 - dat1_s$HG01130_AD_REF/dp)
s1 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("GRCh38-Y"))
head(s1)

dat2_s = dat2[dat2$HG01130_DP>dat2$HG01130_AD_REF,]
dp = dat2_s$HG01130_DP
P_NonRef = 100*(1 - dat2_s$HG01130_AD_REF/dp)
s2 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("T2T-Y"))
s = rbind(s1, s2)
summary(s2$DP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   11.00   14.00   15.09   19.00   42.00 
plot_all(s, "Variants_HG01130_J1", median = 14)

# HG01885 (E1b)
dat1_s = dat1[dat1$HG01885_DP>dat1$HG01885_AD_REF,]
dp = dat1_s$HG01885_DP
P_NonRef = 100*(1 - dat1_s$HG01885_AD_REF/dp)
s1 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("GRCh38-Y"))
head(s1)

dat2_s = dat2[dat2$HG01885_DP>dat2$HG01885_AD_REF,]
dp = dat2_s$HG01885_DP
P_NonRef = 100*(1 - dat2_s$HG01885_AD_REF/dp)
s2 = data.frame(DP = dp, Alt = P_NonRef, Reference = c("T2T-Y"))
s = rbind(s1, s2)
summary(s2$DP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   12.00   15.00   14.85   18.00   58.00 
plot_all(s, "Variants_HG01885_E1b", median = 15)

