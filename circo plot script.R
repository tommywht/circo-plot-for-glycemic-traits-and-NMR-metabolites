# creating all plots 
# circo plots ----------------------
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

setwd("C:/Users/Tommy Wong/Desktop/glycemic traits MR")

## circo plot
circo_plot <- function(dat, var){start_degree = 90
section_track_height = 0.1

discrete_palette <- c("#00378f", # track 1 colour
                      "#ffc067", # track 2 colour
                      "#894300") # track 3 colour

## section header specifics
section_fill_colour <- "snow2"
section_text_colour <- "black"
section_line_colour <- "grey"
section_line_thickness <- 1.5
section_line_type <- 1

## reference lines that go around the tracks
reference_line_colour <- "deeppink"
reference_line_thickness <- 1
reference_line_type <- 1

## point specifics
point_pch <- 21
point_cex <- 1.25

point_col1 <- discrete_palette[1]
point_bg1 <- "white"
point_col1_sig <- "white"
point_bg1_sig <- discrete_palette[1]

point_col2 <- discrete_palette[2]
point_bg2 <- "white"
point_col2_sig <- "white"
point_bg2_sig <- discrete_palette[2]

point_col3 <- discrete_palette[3]
point_bg3 <- "white"
point_col3_sig <- "white"
point_bg3_sig <- discrete_palette[3]

## confidence intervals
ci_lwd <- 1.5
ci_lty <- 1
ci_col1 <- discrete_palette[1]
ci_col2 <- discrete_palette[2]
ci_col3 <- discrete_palette[3]

## lines specifics
lines_col1 <- discrete_palette[1]
lines_col2 <- discrete_palette[2]
lines_col3 <- discrete_palette[3]

lines_lwd <- 3
lines_lty <- 1

## y axis specifics
y_axis_location <- "left"
y_axis_tick <- FALSE
y_axis_tick_length <- 0
y_axis_label_cex <- 0.75

## label specifics
label_distance <- 1.5 # distance from track 0 to plot labels
label_col <- "black"
label_cex <- 0.6

dat$x <- with(dat, 
              ave(seq_along(Group), Group, FUN = seq_along)) # sequence in each group
y <- as.vector(table(dat$Group)) # no. in each group
for(i in 1:nrow(dat)){
  dat$n[i] <- as.numeric(nrow(subset(dat, 
                                     dat$Group == dat$Group[i])))
  dat$ncat[i] <- dat$x[i] / dat$n[i]
}

dat$section_n <- factor(dat$Group,
                        labels = 1:nlevels(dat$Group))

# gap for axis
start_gap <- 10 # this indicate the gap between end and start of circle
gap <- c(rep(1, nlevels(dat$Group)-1), start_gap)

pdf(file = file.path("circo plot", paste0(var, ".pdf")), width = 21, height = 14)
circos.clear()
graphics::par(xpd = NA, cex = 0.8,
              mar = c(0.8, 0.5, 0.3, 0.5) * 26)

circos.par(cell.padding = c(0, 0.5, 0, 0.5),
           start.degree = start_degree,
           gap.degree = gap,
           points.overflow.warning = FALSE,
           track.height = section_track_height,
           clock.wise = TRUE)

# initialise circle
circos.initialize(sector.width = y,
                  xlim = c(0, 1),
                  factors = dat$section_n)

# create and plot section header (the numbers)
circos.track(factors = dat$section_n,
             track.index = 1,
             x = dat$ncat,
             ylim = c(0, 1),
             track.height = 0.075,
             panel.fun = function(x,y){
               chr = circlize::get.cell.meta.data("sector.index") #dont change as this gathers all of the info you need automatically
               xlim = circlize::get.cell.meta.data("xlim") #dont change as this gathers all of the info you need automatically
               ylim = circlize::get.cell.meta.data("ylim") #dont change as this gathers all of the info you need automatically
               circlize::circos.rect(xlim[1], 0, xlim[2], 1, # n (+ and -) length of track away from centre (low number means smaller) - want it large enough to encompass text
                                     border = NA, 
                                     col = section_fill_colour) #colour of track
               circlize::circos.text(mean(xlim), #text location on x axis
                                     mean(ylim), #text location on y axis
                                     chr, # the actual text i.e. sector number
                                     cex = 0.8, #text size
                                     facing = "outside", #diretion of text
                                     niceFacing = TRUE, #flip so text is readable
                                     col = section_text_colour)
             }, bg.border = NA)
# add labels
circlize::circos.trackText(factors = dat$section_n,
                           track.index = 1, #choose labels based on the track we have just made as you can only plot text once a track has been created
                           x = dat$ncat, #location on the x axis where we will plot the lables
                           y = dat$b_ivw *0 + 1.5, #dictates where the labels are plotted - * 0 to give 0 and then choose how far away from the 0 we want to plot the text (this will be trial and error before you get what works best for your data set)
                           # remember sector, x, and y have to be of the same length
                           labels = dat$Subgroup, #where you are taking the names for the labels from
                           facing = "reverse.clockwise",
                           niceFacing = TRUE, #flip the text so it is readable
                           adj = c(1, 1),
                           col = label_col,
                           cex = 0.7) #size of the text

# add data points (track 2)
ci_ivw_min <- min(dat$lci_ivw)
ci_ivw_max <- max(dat$uci_ivw)
track1_axis_min <- round(ci_ivw_min - (ci_ivw_min * 0.1), 2)
track1_axis_max <- round(ci_ivw_max + (ci_ivw_max * 0.1), 2)
ci_egger_min <- min(dat$lci_egger)
ci_egger_max <- max(dat$uci_egger)
track2_axis_min <- round(ci_egger_min - (ci_egger_min * 0.1), 2)
track2_axis_max <- round(ci_egger_max + (ci_egger_max * 0.1), 2)
ci_wm_min <- min(dat$lci_wm)
ci_wm_max <- max(dat$uci_wm)
track3_axis_min <- round(ci_wm_min - (ci_wm_min * 0.1), 2)
track3_axis_max <- round(ci_wm_max + (ci_wm_max * 0.1), 2)

# first track - IVW results
for(i in 1:nlevels(dat$section_n)){
  data1 = subset(dat, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 2,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$b_ivw, #variable you want to plot
                                   ylim = c(track1_axis_min, track1_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$b_ivw * 0 + data1$lci_ivw + 0.01, # y coordinates for start point
                                                               y1 = data1$b_ivw * 0 + data1$uci_ivw, # y coordinates for end point
                                                               col = ci_col1,
                                                               lwd = 3,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# second track - MR-Egger results
for(i in 1:nlevels(dat$section_n)){
  data1 = subset(dat, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 3,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$b_egger, #variable you want to plot
                                   ylim = c(track2_axis_min, track2_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$b_egger * 0 + data1$lci_egger + 0.01, # y coordinates for start point
                                                               y1 = data1$b_egger * 0 + data1$uci_egger, # y coordinates for end point
                                                               col = ci_col2,
                                                               lwd = 3,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# third track - Weighted median results
for(i in 1:nlevels(dat$section_n)){
  data1 = subset(dat, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 4,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$b_wm, #variable you want to plot
                                   ylim = c(track3_axis_min, track3_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$b_wm * 0 + data1$lci_wm + 0.01, # y coordinates for start point
                                                               y1 = data1$b_wm * 0 + data1$uci_wm, # y coordinates for end point
                                                               col = ci_col3,
                                                               lwd = 2.5,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

## 4. layer on top of the confidence intervals the beta
### a. points not reaching significance
pvalue_adjustment <- 0.0009
circlize::circos.trackPoints(factors = subset(dat, dat$pval_ivw > pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(dat, dat$pval_ivw > pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_ivw > pvalue_adjustment)$b_ivw,
                             cex = 1.25,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[1])
circlize::circos.trackPoints(factors = subset(dat, dat$pval_egger > pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(dat, dat$pval_egger > pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_egger > pvalue_adjustment)$b_egger,
                             cex = 1.25,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[2])
circlize::circos.trackPoints(factors = subset(dat, dat$pval_wm > pvalue_adjustment)$section_n,
                             track.index = 4,
                             x = subset(dat, dat$pval_wm > pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_wm > pvalue_adjustment)$b_wm,
                             cex = 1,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[3])
### b. points reaching significance
circlize::circos.trackPoints(factors = subset(dat, dat$pval_ivw <= pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(dat, dat$pval_ivw <= pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_ivw <= pvalue_adjustment)$b_ivw,
                             cex = 1.25,
                             pch = point_pch,
                             col = discrete_palette[1],
                             bg = "white")
circlize::circos.trackPoints(factors = subset(dat, dat$pval_egger <= pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(dat, dat$pval_egger <= pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_egger <= pvalue_adjustment)$b_egger,
                             cex = 1.25,
                             pch = point_pch,
                             col = discrete_palette[2],
                             bg = "white")
circlize::circos.trackPoints(factors = subset(dat, dat$pval_wm <= pvalue_adjustment)$section_n,
                             track.index = 4,
                             x = subset(dat, dat$pval_wm <= pvalue_adjustment)$ncat,
                             y = subset(dat, dat$pval_wm <= pvalue_adjustment)$b_wm,
                             cex = 1,
                             pch = point_pch,
                             col = discrete_palette[3],
                             bg = "white")

## 5. add axis labels - only in the first sector
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 2,
                       at = c(track1_axis_min, 0,  track1_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 3,
                       at = c(track2_axis_min, 0,  track2_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 4,
                       at = c(track3_axis_min, 0,  track3_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)

# Legend
legend1 <- ComplexHeatmap::Legend(at = c("MR_IVW", "MR_Egger", "MR_Weighted median"),
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c(discrete_palette[1], 
                                                                 discrete_palette[2],
                                                                 discrete_palette[3])), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 19, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

# p-value legend points
legend2 <- ComplexHeatmap::Legend(at = "P val <= 0.0009", # breaks, can be wither numeric or character
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 21, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

## Assign legend section labelling
names <- levels(dat$Group)
names <- paste(1:nlevels(dat$Group), names, sep=". ")
legend3 <- ComplexHeatmap::Legend(at = names,
                                  labels_gp = grid::gpar(fontsize = 9),
                                  ncol = 8,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(10, "mm"),
                                  direction = "horizontal")

## 6.D - Pack lagend together
legend4 <- ComplexHeatmap::packLegend(legend1, legend2, direction = "vertical", gap = grid::unit(0, "mm"))
legend <- ComplexHeatmap::packLegend(legend4, legend3, direction = "horizontal", gap = grid::unit(0, "mm"))
legend_height <- legend@grob[["vp"]][["height"]]
legend_width <- legend@grob[["vp"]][["width"]]

# ## 6.E - Layer legend ontop of plot
grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                  y = grid::unit(0.2, "npc"),
                                  width = legend_width,
                                  height = legend_height,
                                  just = c("center", "top")))
grid::grid.draw(legend)
grid::upViewport()


dev.off() }

#### circo plot fasting glucose===============
dat <- readRDS("mr results/fulldat_fg.rds") %>% 
  dplyr::mutate(outcome = tolower(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
dat <- full_join(dat, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)})) %>%
  dplyr::filter(Subgroup != "Glucose")

circo_plot(dat, var = "testing")

#### circo plot 2hr glucose===============
dat <- readRDS("mr results/fulldat_2hglu.rds") %>%
  dplyr::mutate(outcome = tolower(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
dat <- full_join(dat, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)})) %>%
  dplyr::filter(Subgroup != "Glucose")

circo_plot(dat, var = "2-hour glucose")

#### circo plot hba1c ===============
dat <- readRDS("mr results/fulldat_hba1c.rds") %>%
  dplyr::mutate(outcome = tolower(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
dat <- full_join(dat, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)}))

circo_plot(dat, var = "hba1c")

#### circo plot type 2 diabetes ===============
dat <- readRDS("mr results/fulldat_dm.rds") %>%
  dplyr::mutate(outcome = tolower(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
dat <- full_join(dat, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)}))

circo_plot(dat, var = "type 2 diabetes")

#### circo plot fasting insulin ===============
dat <- readRDS("mr results/fulldat_fi.rds") %>%
  dplyr::mutate(outcome = tolower(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
dat <- full_join(dat, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)}))

circo_plot(dat, var = "fasting insulin")

# circo plot for mvmr -----------------------------------------------------
mvmrout <- readRDS("hba1c_hgb_mvmr/mvmr_out.rds") %>% 
  filter(!grepl("ratio ", .$outcome), !grepl("Ratio ", .$outcome)) %>%
  dplyr::mutate(outcome = tolower(outcome)) %>%
  mutate(outcome = str_replace(outcome, "\\|.*$", "")) %>%
  mutate(outcome = str_trim(outcome))

nmrgroup <- openxlsx::read.xlsx("Nightingale NMR grouping.xlsx", sheet = 1) %>% 
  dplyr::mutate(Title = tolower(Title))
mvmrout <- full_join(mvmrout, nmrgroup, by = c("outcome" = "Title")) %>%
  mutate(across(c(Group, Subgroup), .fns = function(x){factor(x, ordered = FALSE)}))

## circo plot for hgb ----------------------------------------------------
start_degree = 90
section_track_height = 0.1

discrete_palette <- c("#00378f", # track 1 colour
                      "#ffc067", # track 2 colour
                      "#894300") # track 3 colour

## section header specifics
section_fill_colour <- "snow2"
section_text_colour <- "black"
section_line_colour <- "grey"
section_line_thickness <- 1.5
section_line_type <- 1

## reference lines that go around the tracks
reference_line_colour <- "deeppink"
reference_line_thickness <- 1.5
reference_line_type <- 1

## point specifics
point_pch <- 21
point_cex <- 1.5

point_col1 <- discrete_palette[1]
point_bg1 <- "white"
point_col1_sig <- "white"
point_bg1_sig <- discrete_palette[1]

point_col2 <- discrete_palette[2]
point_bg2 <- "white"
point_col2_sig <- "white"
point_bg2_sig <- discrete_palette[2]

point_col3 <- discrete_palette[3]
point_bg3 <- "white"
point_col3_sig <- "white"
point_bg3_sig <- discrete_palette[3]

## confidence intervals
ci_lwd <- 2
ci_lty <- 1
ci_col1 <- discrete_palette[1]
ci_col2 <- discrete_palette[2]
ci_col3 <- discrete_palette[3]

## lines specifics
lines_col1 <- discrete_palette[1]
lines_col2 <- discrete_palette[2]
lines_col3 <- discrete_palette[3]

lines_lwd <- 3
lines_lty <- 1

## y axis specifics
y_axis_location <- "left"
y_axis_tick <- FALSE
y_axis_tick_length <- 0
y_axis_label_cex <- 0.75

## label specifics
label_distance <- 1.5 # distance from track 0 to plot labels
label_col <- "black"
label_cex <- 0.6

mvmrout$x <- with(mvmrout, 
                  ave(seq_along(Group), Group, FUN = seq_along)) # sequence in each group
y <- as.vector(table(mvmrout$Group)) # no. in each group
for(i in 1:nrow(mvmrout)){
  mvmrout$n[i] <- as.numeric(nrow(subset(mvmrout, 
                                         mvmrout$Group == mvmrout$Group[i])))
  mvmrout$ncat[i] <- mvmrout$x[i] / mvmrout$n[i]
}

mvmrout$section_n <- factor(mvmrout$Group,
                            labels = 1:nlevels(mvmrout$Group))

# gap for axis
start_gap <- 10 # this indicate the gap between end and start of circle
gap <- c(rep(1, nlevels(mvmrout$Group)-1), start_gap)

dev.off()
pdf(file = "circo plot/test.pdf", width = 21, height = 14)
circos.clear()
graphics::par(xpd = NA, cex = 0.8,
              mar = c(0.8, 0.5, 0.3, 0.5) * 26)

circos.par(cell.padding = c(0, 0.5, 0, 0.5),
           start.degree = start_degree,
           gap.degree = gap,
           points.overflow.warning = FALSE,
           track.height = section_track_height,
           clock.wise = TRUE)

# initialise circle
circos.initialize(sector.width = y,
                  xlim = c(0, 1),
                  factors = mvmrout$section_n)

# create and plot section header (the numbers)
circos.track(factors = mvmrout$section_n,
             track.index = 1,
             x = mvmrout$ncat,
             ylim = c(0, 1),
             track.height = 0.075,
             panel.fun = function(x,y){
               chr = circlize::get.cell.meta.data("sector.index") #dont change as this gathers all of the info you need automatically
               xlim = circlize::get.cell.meta.data("xlim") #dont change as this gathers all of the info you need automatically
               ylim = circlize::get.cell.meta.data("ylim") #dont change as this gathers all of the info you need automatically
               circlize::circos.rect(xlim[1], 0, xlim[2], 1, # n (+ and -) length of track away from centre (low number means smaller) - want it large enough to encompass text
                                     border = NA, 
                                     col = section_fill_colour) #colour of track
               circlize::circos.text(mean(xlim), #text location on x axis
                                     mean(ylim), #text location on y axis
                                     chr, # the actual text i.e. sector number
                                     cex = 0.8, #text size
                                     facing = "outside", #diretion of text
                                     niceFacing = TRUE, #flip so text is readable
                                     col = section_text_colour)
             }, bg.border = NA)
# add labels
circlize::circos.trackText(factors = mvmrout$section_n,
                           track.index = 1, #choose labels based on the track we have just made as you can only plot text once a track has been created
                           x = mvmrout$ncat, #location on the x axis where we will plot the lables
                           y = mvmrout$hgb_ivw_beta *0 + 1.5, #dictates where the labels are plotted - * 0 to give 0 and then choose how far away from the 0 we want to plot the text (this will be trial and error before you get what works best for your data set)
                           # remember sector, x, and y have to be of the same length
                           labels = mvmrout$Subgroup, #where you are taking the names for the labels from
                           facing = "reverse.clockwise",
                           niceFacing = TRUE, #flip the text so it is readable
                           adj = c(1, 1),
                           col = label_col,
                           cex = 0.7) #size of the text

# add data points (track 2)
ci_ivw_min <- min(mvmrout$hgb_ivw_lci)
ci_ivw_max <- max(mvmrout$hgb_ivw_uci)
track1_axis_min <- round(ci_ivw_min - (ci_ivw_min * 0.1), 2)
track1_axis_max <- round(ci_ivw_max + (ci_ivw_max * 0.1), 2)
ci_egger_min <- min(mvmrout$hgb_egger_lci)
ci_egger_max <- max(mvmrout$hgb_egger_uci)
track2_axis_min <- round(ci_egger_min - (ci_egger_min * 0.1), 2)
track2_axis_max <- round(ci_egger_max + (ci_egger_max * 0.1), 2)

# first track - IVW results
for(i in 1:nlevels(mvmrout$section_n)){
  data1 = subset(mvmrout, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 2,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$hgb_ivw_beta, #variable you want to plot
                                   ylim = c(track1_axis_min, track1_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$hgb_ivw_beta * 0 + data1$hgb_ivw_lci, # y coordinates for start point
                                                               y1 = data1$hgb_ivw_beta * 0 + data1$hgb_ivw_uci, # y coordinates for end point
                                                               col = ci_col1,
                                                               lwd = 4,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# second track - MR-Egger results
for(i in 1:nlevels(mvmrout$section_n)){
  data1 = subset(mvmrout, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 3,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$hgb_egger_beta, #variable you want to plot
                                   ylim = c(track2_axis_min, track2_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$hgb_egger_beta * 0 + data1$hgb_egger_lci, # y coordinates for start point
                                                               y1 = data1$hgb_egger_beta * 0 + data1$hgb_egger_uci, # y coordinates for end point
                                                               col = ci_col2,
                                                               lwd = 3.5,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# third track - Weighted median results - no weighted median results for mvmr

## 4. layer on top of the confidence intervals the beta
### a. points not reaching significance
pvalue_adjustment <- 0.0009
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hgb_ivw_pval > pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(mvmrout, mvmrout$hgb_ivw_pval > pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hgb_ivw_pval > pvalue_adjustment)$hgb_ivw_beta,
                             cex = 1.5,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[1])
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hgb_egger_pval > pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(mvmrout, mvmrout$hgb_egger_pval > pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hgb_egger_pval > pvalue_adjustment)$hgb_egger_beta,
                             cex = 1.25,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[2])

### b. points reaching significance
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hgb_ivw_pval <= pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(mvmrout, mvmrout$hgb_ivw_pval <= pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hgb_ivw_pval <= pvalue_adjustment)$hgb_ivw_beta,
                             cex = 1.5,
                             pch = point_pch,
                             col = discrete_palette[1],
                             bg = "white")
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hgb_egger_pval <= pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(mvmrout, mvmrout$hgb_egger_pval <= pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hgb_egger_pval <= pvalue_adjustment)$hgb_egger_beta,
                             cex = 1.25,
                             pch = point_pch,
                             col = discrete_palette[2],
                             bg = "white")

## 5. add axis labels - only in the first sector
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 2,
                       at = c(track1_axis_min, 0,  track1_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 3,
                       at = c(track2_axis_min, 0,  track2_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)

# Legend
legend1 <- ComplexHeatmap::Legend(at = c("MR_IVW", "MR_Egger"),
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c(discrete_palette[1], 
                                                                 discrete_palette[2])), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 19, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

# p-value legend points
legend2 <- ComplexHeatmap::Legend(at = "P val <= 0.0009", # breaks, can be wither numeric or character
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 21, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

## Assign legend section labelling
names <- levels(mvmrout$Group)
names <- paste(1:nlevels(mvmrout$Group), names, sep=". ")
legend3 <- ComplexHeatmap::Legend(at = names,
                                  labels_gp = grid::gpar(fontsize = 9),
                                  ncol = 8,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(10, "mm"),
                                  direction = "horizontal")

## 6.D - Pack lagend together
legend4 <- ComplexHeatmap::packLegend(legend1, legend2, direction = "vertical", gap = grid::unit(0, "mm"))
legend <- ComplexHeatmap::packLegend(legend4, legend3, direction = "horizontal", gap = grid::unit(0, "mm"))
legend_height <- legend@grob[["vp"]][["height"]]
legend_width <- legend@grob[["vp"]][["width"]]

# ## 6.E - Layer legend ontop of plot
grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                  y = grid::unit(0.2, "npc"),
                                  width = legend_width,
                                  height = legend_height,
                                  just = c("center", "top")))
grid::grid.draw(legend)
grid::upViewport()


dev.off() 

## circo plot for hba1c ----------------------------------------------------
start_degree = 90
section_track_height = 0.1

discrete_palette <- c("#00378f", # track 1 colour
                      "#ffc067", # track 2 colour
                      "#894300") # track 3 colour

## section header specifics
section_fill_colour <- "snow2"
section_text_colour <- "black"
section_line_colour <- "grey"
section_line_thickness <- 1.5
section_line_type <- 1

## reference lines that go around the tracks
reference_line_colour <- "deeppink"
reference_line_thickness <- 1.5
reference_line_type <- 1

## point specifics
point_pch <- 21
point_cex <- 1.25

point_col1 <- discrete_palette[1]
point_bg1 <- "white"
point_col1_sig <- "white"
point_bg1_sig <- discrete_palette[1]

point_col2 <- discrete_palette[2]
point_bg2 <- "white"
point_col2_sig <- "white"
point_bg2_sig <- discrete_palette[2]

point_col3 <- discrete_palette[3]
point_bg3 <- "white"
point_col3_sig <- "white"
point_bg3_sig <- discrete_palette[3]

## confidence intervals
ci_lwd <- 1.5
ci_lty <- 1
ci_col1 <- discrete_palette[1]
ci_col2 <- discrete_palette[2]
ci_col3 <- discrete_palette[3]

## lines specifics
lines_col1 <- discrete_palette[1]
lines_col2 <- discrete_palette[2]
lines_col3 <- discrete_palette[3]

lines_lwd <- 3
lines_lty <- 1

## y axis specifics
y_axis_location <- "left"
y_axis_tick <- FALSE
y_axis_tick_length <- 0
y_axis_label_cex <- 0.75

## label specifics
label_distance <- 1.5 # distance from track 0 to plot labels
label_col <- "black"
label_cex <- 0.6

mvmrout$x <- with(mvmrout, 
                  ave(seq_along(Group), Group, FUN = seq_along)) # sequence in each group
y <- as.vector(table(mvmrout$Group)) # no. in each group
for(i in 1:nrow(mvmrout)){
  mvmrout$n[i] <- as.numeric(nrow(subset(mvmrout, 
                                         mvmrout$Group == mvmrout$Group[i])))
  mvmrout$ncat[i] <- mvmrout$x[i] / mvmrout$n[i]
}

mvmrout$section_n <- factor(mvmrout$Group,
                            labels = 1:nlevels(mvmrout$Group))

# gap for axis
start_gap <- 10 # this indicate the gap between end and start of circle
gap <- c(rep(1, nlevels(mvmrout$Group)-1), start_gap)

dev.off()
pdf(file = "circo plot/hba1c_mvmr.pdf", width = 21, height = 14)
circos.clear()
graphics::par(xpd = NA, cex = 0.8,
              mar = c(0.8, 0.5, 0.3, 0.5) * 26)

circos.par(cell.padding = c(0, 0.5, 0, 0.5),
           start.degree = start_degree,
           gap.degree = gap,
           points.overflow.warning = FALSE,
           track.height = section_track_height,
           clock.wise = TRUE)

# initialise circle
circos.initialize(sector.width = y,
                  xlim = c(0, 1),
                  factors = mvmrout$section_n)

# create and plot section header (the numbers)
circos.track(factors = mvmrout$section_n,
             track.index = 1,
             x = mvmrout$ncat,
             ylim = c(0, 1),
             track.height = 0.075,
             panel.fun = function(x,y){
               chr = circlize::get.cell.meta.data("sector.index") #dont change as this gathers all of the info you need automatically
               xlim = circlize::get.cell.meta.data("xlim") #dont change as this gathers all of the info you need automatically
               ylim = circlize::get.cell.meta.data("ylim") #dont change as this gathers all of the info you need automatically
               circlize::circos.rect(xlim[1], 0, xlim[2], 1, # n (+ and -) length of track away from centre (low number means smaller) - want it large enough to encompass text
                                     border = NA, 
                                     col = section_fill_colour) #colour of track
               circlize::circos.text(mean(xlim), #text location on x axis
                                     mean(ylim), #text location on y axis
                                     chr, # the actual text i.e. sector number
                                     cex = 0.8, #text size
                                     facing = "outside", #diretion of text
                                     niceFacing = TRUE, #flip so text is readable
                                     col = section_text_colour)
             }, bg.border = NA)
# add labels
circlize::circos.trackText(factors = mvmrout$section_n,
                           track.index = 1, #choose labels based on the track we have just made as you can only plot text once a track has been created
                           x = mvmrout$ncat, #location on the x axis where we will plot the lables
                           y = mvmrout$hba1c_ivw_beta *0 + 1.5, #dictates where the labels are plotted - * 0 to give 0 and then choose how far away from the 0 we want to plot the text (this will be trial and error before you get what works best for your data set)
                           # remember sector, x, and y have to be of the same length
                           labels = mvmrout$Subgroup, #where you are taking the names for the labels from
                           facing = "reverse.clockwise",
                           niceFacing = TRUE, #flip the text so it is readable
                           adj = c(1, 1),
                           col = label_col,
                           cex = 0.7) #size of the text

# add data points (track 2)
ci_ivw_min <- min(mvmrout$hba1c_ivw_lci)
ci_ivw_max <- max(mvmrout$hba1c_ivw_uci)
track1_axis_min <- round(ci_ivw_min - (ci_ivw_min * 0.1), 2)
track1_axis_max <- round(ci_ivw_max + (ci_ivw_max * 0.1), 2)
ci_egger_min <- min(mvmrout$hba1c_egger_lci)
ci_egger_max <- max(mvmrout$hba1c_egger_uci)
track2_axis_min <- round(ci_egger_min - (ci_egger_min * 0.1), 2)
track2_axis_max <- round(ci_egger_max + (ci_egger_max * 0.1), 2)

# first track - IVW results
for(i in 1:nlevels(mvmrout$section_n)){
  data1 = subset(mvmrout, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 2,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$hba1c_ivw_beta, #variable you want to plot
                                   ylim = c(track1_axis_min, track1_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$hba1c_ivw_beta * 0 + data1$hba1c_ivw_lci + 0.01, # y coordinates for start point
                                                               y1 = data1$hba1c_ivw_beta * 0 + data1$hba1c_ivw_uci + 0.01, # y coordinates for end point
                                                               col = ci_col1,
                                                               lwd = 3,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# second track - MR-Egger results
for(i in 1:nlevels(mvmrout$section_n)){
  data1 = subset(mvmrout, section_n == i)
  
  circlize::circos.trackPlotRegion(factors = data1$section_n, #we plot the first region based on the column in our data set which we create the sections from
                                   track.index = 3,
                                   x = data1$ncat, #set this as ncat as ncat dictates the location of the variable you want to plot within each section and within the circle as a whole
                                   y = data1$hba1c_egger_beta, #variable you want to plot
                                   ylim = c(track2_axis_min, track2_axis_max), #co-ordinates of the Y axis of the track
                                   track.height = 0.2, #how big is the track as % of circle
                                   
                                   #Set sector background
                                   bg.border = NA,
                                   bg.col = NA,
                                   
                                   #Map values
                                   panel.fun = function(x, y) { #this sets x and y as the above defined variables for the following lines of code
                                     
                                     # plot '0' reference line
                                     circlize::circos.lines(x = x,
                                                            y = y * 0,
                                                            col = reference_line_colour, #set the 0 line colour to something distinctive
                                                            lwd = reference_line_thickness, #set the thickness of the line so it's a bit smaller than your point (looks better)
                                                            lty = reference_line_type)
                                     
                                     # confidence interval
                                     circlize::circos.segments(x0 = data1$ncat, # x coordinates for starting point
                                                               x1 = data1$ncat, # x coordinates for end point
                                                               y0 = data1$hba1c_egger_beta * 0 + data1$hba1c_egger_lci, # y coordinates for start point
                                                               y1 = data1$hba1c_egger_beta * 0 + data1$hba1c_egger_uci, # y coordinates for end point
                                                               col = ci_col2,
                                                               lwd = 2.5,
                                                               lty = ci_lty,
                                                               sector.index = i)})}

# third track - Weighted median results - no weighted median results for mvmr

## 4. layer on top of the confidence intervals the beta
### a. points not reaching significance
pvalue_adjustment <- 0.0009
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hba1c_ivw_pval > pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(mvmrout, mvmrout$hba1c_ivw_pval > pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hba1c_ivw_pval > pvalue_adjustment)$hba1c_ivw_beta,
                             cex = 1.5,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[1])
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hba1c_egger_pval > pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(mvmrout, mvmrout$hba1c_egger_pval > pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hba1c_egger_pval > pvalue_adjustment)$hba1c_egger_beta,
                             cex = 1.25,
                             pch = point_pch,
                             col = "white",
                             bg = discrete_palette[2])

### b. points reaching significance
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hba1c_ivw_pval <= pvalue_adjustment)$section_n,
                             track.index = 2,
                             x = subset(mvmrout, mvmrout$hba1c_ivw_pval <= pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hba1c_ivw_pval <= pvalue_adjustment)$hba1c_ivw_beta,
                             cex = 1.5,
                             pch = point_pch,
                             col = discrete_palette[1],
                             bg = "white")
circlize::circos.trackPoints(factors = subset(mvmrout, mvmrout$hba1c_egger_pval <= pvalue_adjustment)$section_n,
                             track.index = 3,
                             x = subset(mvmrout, mvmrout$hba1c_egger_pval <= pvalue_adjustment)$ncat,
                             y = subset(mvmrout, mvmrout$hba1c_egger_pval <= pvalue_adjustment)$hba1c_egger_beta,
                             cex = 1.25,
                             pch = point_pch,
                             col = discrete_palette[2],
                             bg = "white")

## 5. add axis labels - only in the first sector
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 2,
                       at = c(track1_axis_min, 0,  track1_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)
circlize::circos.yaxis(side = "left",
                       sector.index = 1, #the sector this is plotted in
                       track.index = 3,
                       at = c(track2_axis_min, 0,  track2_axis_max), #location on the y axis as well as the name of the label
                       tick = FALSE, tick.length = 0,
                       labels.cex = 0.75)

# Legend
legend1 <- ComplexHeatmap::Legend(at = c("MR_IVW", "MR_Egger"),
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c(discrete_palette[1], 
                                                                 discrete_palette[2])), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 19, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

# p-value legend points
legend2 <- ComplexHeatmap::Legend(at = "P val <= 0.0009", # breaks, can be wither numeric or character
                                  labels_gp = grid::gpar(fontsize = 10),
                                  ncol = 1,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  type = "points", # type of legends, can be grid, points and lines
                                  pch = 21, # type of points
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(15, "mm"),
                                  direction = "vertical")

## Assign legend section labelling
names <- levels(mvmrout$Group)
names <- paste(1:nlevels(mvmrout$Group), names, sep=". ")
legend3 <- ComplexHeatmap::Legend(at = names,
                                  labels_gp = grid::gpar(fontsize = 9),
                                  ncol = 8,
                                  border = NA, # color of legend borders, also for the ticks in the continuous legend
                                  background = NA, # background colors
                                  legend_gp = grid::gpar(col = c("black")), # graphic parameters for the legend
                                  size = grid::unit(15, "mm"), # size of points
                                  grid_height	= grid::unit(15, "mm"),
                                  grid_width = grid::unit(10, "mm"),
                                  direction = "horizontal")

## 6.D - Pack lagend together
legend4 <- ComplexHeatmap::packLegend(legend1, legend2, direction = "vertical", gap = grid::unit(0, "mm"))
legend <- ComplexHeatmap::packLegend(legend4, legend3, direction = "horizontal", gap = grid::unit(0, "mm"))
legend_height <- legend@grob[["vp"]][["height"]]
legend_width <- legend@grob[["vp"]][["width"]]

# ## 6.E - Layer legend ontop of plot
grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                  y = grid::unit(0.2, "npc"),
                                  width = legend_width,
                                  height = legend_height,
                                  just = c("center", "top")))
grid::grid.draw(legend)
grid::upViewport()


dev.off() 