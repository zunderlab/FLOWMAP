rm(list = ls(all = TRUE))
library(igraph)

GRAPHML_DIR = "."
GRAPHML_PATTERN <- "_modularity.graphml$"
PCTILE_COLOR = c(0.02, 0.98)
SCALE = NULL
NORMALIZE = "global"
NODE_SIZE_SCALE_FACTOR = 1
MIN_NODE_SIZE = 10
MAX_NODE_SIZE = 20
PNG_WIDTH = 2400
PNG_HEIGHT = 1800
TEXT_COLOR = "black"
EDGE.COLOR = "grey"
BARE = FALSE
PALETTE = "bluered"
TEST = FALSE

graphml_files <- list.files(path=GRAPHML_DIR,pattern=GRAPHML_PATTERN,full.names=TRUE)
for (gr_file in graphml_files) {
  
  # set output directory
  out_dir = paste("./",sub(GRAPHML_PATTERN,"_png/",basename(gr_file)),sep="")
  dir.create(out_dir)
  
  
  # read in graph file
  gr <- read.graph(gr_file,format = "graphml")
  
  # get matrix/table of all vertex attribute values
  attrs <- c()
  attrs_colnames <- c()
  if (TEST) {
    a="percenttotal"
    #attrs <- cbind(attrs,as.numeric(get.vertex.attribute(gr,a,index=V(gr))))
    attrs <- matrix(data=as.numeric(get.vertex.attribute(gr,a,index=V(gr))),ncol=1)
    attrs_colnames <- c(attrs_colnames,a)
  } else {
    for (a in list.vertex.attributes(gr)) {
      if (is.numeric(get.vertex.attribute(gr,a,index=V(gr)))) {
      attrs <- cbind(attrs,as.numeric(get.vertex.attribute(gr,a,index=V(gr))))
      attrs_colnames <- c(attrs_colnames,a)
      }
    }
  }
  colnames(attrs) <- attrs_colnames
  rownames(attrs) <- V(gr)
  
  # get x-y graph layout
  graph_l <- matrix(data=c(V(gr)$x,V(gr)$y),nrow=length(V(gr)$x),ncol=2)
  
  # for global normalization 
  boundaries <- NULL
  if (NORMALIZE == "global") {
    boundaries <- c()
    all_attrs <- c()
    for (f in graphml_files) {
      tmp_gr <- read.graph(f,format = "graphml")
      tmp_attrs <- c()
      tmp_attrs_colnames <- c()
      for (a in list.vertex.attributes(tmp_gr)) {
        if (is.numeric(get.vertex.attribute(gr,a,index=V(gr)))) {
          tmp_attrs <- cbind(tmp_attrs,as.numeric(get.vertex.attribute(tmp_gr,
            a,index=V(tmp_gr))))
          tmp_attrs_colnames <- c(tmp_attrs_colnames,a)
        }
      }
      colnames(tmp_attrs) <- tmp_attrs_colnames
      rownames(tmp_attrs) <- V(tmp_gr)
      for (c in colnames(tmp_attrs)) {
        all_attrs[[c]] <- c(all_attrs[[c]], tmp_attrs[, c])
      }
    }
    for (i in seq_along(all_attrs)) {
      boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], 
        probs = PCTILE_COLOR, na.rm = TRUE)
    }
  }
  
  # set up color scale
  if (PALETTE == "jet")
    palette <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
  if (PALETTE == "bluered")
    palette <- colorRampPalette(c("blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red"))
  colorSCALE <- palette(100)
  
  # set up node size
  vsize <- attrs[,"percenttotal"]
  vsize <- (vsize-min(vsize, na.rm = TRUE))/(max(vsize, na.rm = TRUE)^(1/NODE_SIZE_SCALE_FACTOR)) * 
    ((MAX_NODE_SIZE)^0.5/pi) + ((MIN_NODE_SIZE)^0.5/pi)
  vsize[is.na(vsize) | (attrs[,"percenttotal"] == 0)] <- (MIN_NODE_SIZE)^0.5/pi
  
  # print out one pdf for each attribute
  for (name in colnames(attrs)) {
    #get attribute name and data
    #name <- colnames(attrs)[i]
    attr <- attrs[, name]
    
    # set up color boundaries
    ifelse (!is.null(SCALE), 
      boundary <- SCALE,
      ifelse (NORMALIZE == "global",
        boundary <- boundaries[[name]],
          boundary <- quantile(attr, probs = PCTILE_COLOR, na.rm = TRUE)
      )
    )
    ifelse (length(grep("^medians|percent|cvs|Dd|Timepoint", name)), 
      boundary <- c(min(boundary), max(boundary)),
        boundary <- c(-max(abs(boundary)), max(abs(boundary)))
    )
    boundary <- round(boundary, 2)
    if (boundary[1] == boundary[2]) {
      boundary <- c(boundary[1] - 1, boundary[2] + 1)
    }
    grad <- seq(boundary[1], boundary[2], length.out = length(colorSCALE))
    color <- colorSCALE[findInterval(attr, grad, all.inside = TRUE)]
    color[is.na(attr) | (attrs[,"percenttotal"] == 0)] <- "grey"
    if (grepl("^percenttotalratiolog$", name)) {
      color[is.na(attr) & attrs[,"percenttotal"] > 0] <- tail(colorSCALE, 1)
    }
    fill_color <- color
    is.na(fill_color) <- is.na(attr)
    frame_color <- color
    png(filename=paste(out_dir, sub("/","-",name), ".png", sep = ""),
        width = PNG_WIDTH, height = PNG_HEIGHT, units = "px", pointsize = 12,
        bg = "transparent")
    graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 
      2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
    par(mar=c(1.5,0,0,0))
    plot(gr, layout = graph_l, vertex.shape = "circle", 
      vertex.color = fill_color, vertex.frame.color = frame_color, 
      edge.color = EDGE.COLOR, vertex.size = vsize, 
      vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
      asp = graph_aspect)
    if (!BARE) {
      if (length(grep("^medians", name))) 
        name <- sub("medians", "Median of ", name)
      if (length(grep("^fold", name))) 
        name <- sub("fold", "Arcsinh diff. of ", name)
      if (grepl("^percenttotal$", name)) 
        name <- sub("percent", "Percent freq. of ", name)
      if (grepl("^percenttotalratiolog$", name)) 
        name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
      if (grepl("^cvs", name)) 
        name <- sub("cvs", "Coeff. of Variation of ", name)
      if (grepl("_clust$", name)) 
        name <- sub("_clust", "\n(Used for tree-building)", name)
      
      #title(main = paste(strsplit(basename(gr_file), GRAPHML_PATTERN)[[1]][1], sub = name, sep = "\n"),
            #col.main="white")
      title(sub=sub("\\(.*","",name),line=0,col.sub=TEXT_COLOR,cex.sub=2)
      #par(col.lab="white",col.axis="white",mar=c(1,1,2,1))
      #subplot(image(grad, c(1), matrix(1:length(colorSCALE), 
        #ncol = 1), col = colorSCALE, xlab = ifelse(is.null(SCALE), 
        #paste("Range:", PCTILE_COLOR[1], "to", PCTILE_COLOR[2], 
        #"pctile"), ""), ylab = "", yaxt = "n", xaxp = c(boundary, 
        #1),main=sub("\\(.*","",name),col.main="white"), x = "right,bottom",
        #size = c(1, 0.2))
    }
    dev.off()
  }
}
