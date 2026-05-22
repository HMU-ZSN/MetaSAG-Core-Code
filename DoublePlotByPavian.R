

###############Reference the drawing code of pavian in the R package###############
#https://github.com/fbreitwieser/pavian/tree/master/R
#pavian/R/sample-build_sankey_network.R
#pavian/R/datainput-read_report.R
#####################################################

library(getopt)
library(sankeyD3)
spec <- matrix(
  c("Kraken2Report",  "i", 2, "character", "This is kraken report file.",
    "OutPlotHtml", "o", 2, "character",  "This is pavian plot html"),
  byrow=TRUE, ncol=5)

# use getopt
opt <- getopt(spec)
my_report=opt$Kraken2Report
OutHtml=opt$OutPlotHtml

print(my_report)
print(OutHtml)









build_kraken_tree <- function(report) {
  if (nrow(report) == 0 || nrow(report) == 1) {
    # this should only happen if the original input to the function has a size <= 1
    return(list(report))
  }
  
  ## select the current depth as the one of the topmost data.frame row
  sel_depth <- report[,'depth'] == report[1,'depth']
  
  ## partition the data.frame into parts with that depth
  depth_partitions <- cumsum(sel_depth)
  
  ## for each depth partition
  res <- lapply(unique(depth_partitions),
                function(my_depth_partition) {
                  sel <- depth_partitions == my_depth_partition
                  
                  ## return the data.frame row if it is only one row (leaf node, ends recursion)
                  if (sum(sel) == 1)
                    return(report[sel,,drop=F])
                  
                  ## otherwise: take first row as partition descriptor ..
                  first_row <- which(sel)[1]
                  ##  and recurse deeper into the depths with the remaining rows
                  dres <- build_kraken_tree(report[which(sel)[-1],,drop=F])
                  
                  attr(dres,"row") <- report[first_row,,drop=F]
                  dres
                })
  names(res) <- report$name[sel_depth]
  res
}


## Collapse taxonomic taxRanks to only those mentioned in keep_taxRanks
collapse.taxRanks <- function(krakenlist,keep_taxRanks=LETTERS,filter_taxon=NULL) {
  ## input: a list, in which each element is either a
  ##            a list or a data.frame (for the leafs)
  ##   the input has an attribute row that gives details on the current taxRank
  
  ## columns whose values are added to the next taxRank when
  ##  a taxRank is deleted
  cols <- c("taxonReads","n_unique_kmers","n_kmers")
  if (length(krakenlist) == 0 || is.data.frame(krakenlist)) {
    return(krakenlist)
  }
  
  parent_row <- attr(krakenlist,"row")
  all.child_rows <- c()
  
  if (is.null(parent_row)) {
    return(do.call(rbind,lapply(krakenlist,collapse.taxRanks,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)))
  }
  
  ## rm.cladeReads captures the number of cladeReads that are deleted.
  ##  this has to be propagated to the higher taxRank
  rm.cladeReads <- 0
  
  for (kl in krakenlist) {
    if (is.data.frame(kl)) {  ## is a leaf node?
      child_rows <- kl
    } else {                 ## recurse deeper into tree
      child_rows <- collapse.taxRanks(kl,keep_taxRanks,filter_taxon=filter_taxon)
      if ('rm.cladeReads' %in% names(attributes(child_rows))) {
        rm.cladeReads <- rm.cladeReads + attr(child_rows,'rm.cladeReads')
      }
    }
    
    ## check if this taxRank and the taxRanks below should be removed
    delete.taxon <- child_rows[1,'name'] %in% filter_taxon
    if (delete.taxon) {
      rm.cladeReads <- rm.cladeReads + child_rows[1,'cladeReads']
      dmessage(sprintf("removed %7s cladeReads, including %s childs, for %s",child_rows[1,'"cladeReads"'],nrow(child_rows)-1,child_rows[1,'name']))
      
      ## remove all children
      child_rows <- NULL
      
    } else {
      
      ## check if the last (top-most) row should be kept
      keep_last.child <- child_rows[1,'taxRank'] %in% keep_taxRanks
      
      if (!keep_last.child) {
        cols <- cols[cols %in% colnames(parent_row)]
        
        ## save the specified colum information to the parent
        parent_row[,cols] <- parent_row[,cols] + child_rows[1,cols]
        
        ## remove row
        child_rows <- child_rows[-1,,drop=FALSE]
        
        ## decrease depths of rows below child row
        if (nrow(child_rows) > 0)
          child_rows[,'depth'] <- child_rows[,'depth'] - 1
        
      }
    }
    all.child_rows <- rbind(all.child_rows,child_rows)
  }
  
  ## subtract deleted read count from parent row
  parent_row[,'cladeReads'] <- parent_row[,'cladeReads'] - rm.cladeReads
  res <- rbind(parent_row,all.child_rows)
  
  if (parent_row[,'cladeReads'] < 0)
    stop("mistake made in removing cladeReads")
  #if (parent_row[,'"cladeReads"'] == 0)
  #  res <- c()
  
  if (rm.cladeReads > 0)
    attr(res,'rm.cladeReads') <- rm.cladeReads
  return(res)
}

delete_taxRanks_below <- function(report,taxRank="S") {
  del_taxRank <- 0
  do_del <- FALSE
  del_row <- 0
  
  cols <- c("taxonReads","n_unique_kmers","n_kmers")
  sub.sums <- c(0,0,0)
  
  rows_to_delete <- c()
  for (i in seq_len(nrow(report))) {
    if (report[i,'taxRank'] %in% taxRank) {
      del_depth <- report[i,'depth']
      do_del <- TRUE
      del_row <- i
      sub.sums <- c(0,0,0)
    } else {
      if (do_del) {
        if (report[i,'depth'] > del_taxRank) {
          rows_to_delete <- c(rows_to_delete,i)
          sub.sums <- sub.sums + report[i,cols]
        } else {
          report[del_row,cols] <- report[del_row,cols]+sub.sums
          sub.sums <- c(0,0,0)
          do_del <- FALSE
        }
      }
    }
  }
  report[-rows_to_delete,]
}

#' Read kraken or centrifuge-style report
#'
#' @param myfile kraken report file
#' @param collapse  should the results be collapsed to only those taxRanks specified in keep_taxRanks?
#' @param keep_taxRanks taxRanks to keep when collapse is TRUE
#' @param min.depth minimum depth
#' @param filter_taxon filter certain taxon names
#' @param has_header if the kraken report has a header or not
#' @param add_taxRank_columns if TRUE, for each taxRank columns are added
#'
#' @return report data.frame
#' @export
#'
read_report2 <- function(myfile,collapse=TRUE,keep_taxRanks=c("D","K","P","C","O","F","G","S"),min.depth=0,filter_taxon=NULL,
                         has_header=NULL,add_taxRank_columns=FALSE) {
  
  first.line <- readLines(myfile,n=1)
  isASCII <-  function(txt) all(charToRaw(txt) <= as.raw(127))
  if (!isASCII(first.line)) {
    dmessage(myfile," is no valid report - not all characters are ASCII")
    return(NULL)
  }
  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z]",first.line)
  }
  
  if (has_header) {
    report <- utils::read.table(myfile,sep="\t",header = T,
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
    #colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")
    
    ## harmonize column names. TODO: Harmonize them in the scripts!
    colnames(report)[colnames(report)=="clade_perc"] <- "percentage"
    colnames(report)[colnames(report)=="perc"] <- "percentage"
    
    colnames(report)[colnames(report)=="n_reads_clade"] <- "cladeReads"
    colnames(report)[colnames(report)=="n.clade"] <- "cladeReads"
    
    colnames(report)[colnames(report)=="n_reads_taxo"] <- "taxonReads"
    colnames(report)[colnames(report)=="n.stay"] <- "taxonReads"
    
    colnames(report)[colnames(report)=="rank"] <- "taxRank"
    colnames(report)[colnames(report)=="tax_rank"] <- "taxRank"
    
    colnames(report)[colnames(report)=="taxonid"] <- "taxID"
    colnames(report)[colnames(report)=="tax"] <- "taxID"
    
  } else {
    report <- utils::read.table(myfile,sep="\t",header = F,
                                col.names = c("percentage","cladeReads","taxonReads","taxRank","taxID","name"),
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
  }
  
  report$depth <- nchar(gsub("\\S.*","",report$name))/2
  report$name <- gsub("^ *","",report$name)
  report$name <- paste(tolower(report$taxRank),report$name,sep="_")
  
  
  ## Only stop at certain taxRanks
  ## filter taxon and further up the tree if 'filter_taxon' is defined
  kraken.tree <- build_kraken_tree(report)
  report <- collapse.taxRanks(kraken.tree,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)
  
  ## Add a metaphlan-style taxon string
  if (add_taxRank_columns) {
    report[,keep_taxRanks] <- NA
  }
  report$taxLineage = report$name
  rows_to_consider <- rep(FALSE,nrow(report))
  
  for (i in seq_len(nrow(report))) {
    ## depth > 2 correspond to taxRanks below 'D'
    if (i > 1 && report[i,"depth"] > min.depth) {
      ## find the maximal index of a row below the current depth
      idx <- report$depth < report[i,"depth"] & rows_to_consider
      if (!any(idx)) { next() }
      
      current.taxRank <- report[i,'taxRank']
      my_row <- max(which(idx))
      report[i,'taxLineage'] <- paste(report[my_row,'taxLineage'],report[i,'taxLineage'],sep="|")
      
      if (add_taxRank_columns) {
        if (report[my_row,'taxRank'] %in% keep_taxRanks) {
          taxRanks.cp <- keep_taxRanks[seq(from=1,to=which(keep_taxRanks == report[my_row,'taxRank']))]
          report[i,taxRanks.cp] <- report[my_row,taxRanks.cp]
        }
        
        report[i,report[i,'taxRank']] <- report[i,'name']
      }
    }
    rows_to_consider[i] <- TRUE
  }
  
  report <- report[report$depth >= min.depth,]
  
  report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100
  
  for (column in c("taxonReads", "cladeReads")) 
    if (all(floor(report[[column]]) == report[[column]])) 
      report[[column]] <- as.integer(report[[column]])
  
  if ('n_unique_kmers'  %in% colnames(report))
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
  #report$taxRankperc <- 100/taxRank(report$cladeReads)
  
  rownames(report) <- NULL
  
  report
}


build_sankey_network <- function(my_report, taxRanks =  c("D","K","P","C","O","F","G","S"), maxn=10,
                                 zoom = F, title = NULL,
                                 ...) {
  stopifnot("taxRank" %in% colnames(my_report))
  if (!any(taxRanks %in% my_report$taxRank)) {
    warning("report does not contain any of the taxRanks - skipping it")
    return()
  }
  my_report <- subset(my_report, taxRank %in% taxRanks)
  my_report <- plyr::ddply(my_report, "taxRank", function(x) x[utils::tail(order(x$cladeReads,-x$depth), n=maxn), , drop = FALSE])
  
  my_report <- my_report[, c("name","taxLineage","taxonReads", "cladeReads","depth", "taxRank")]
  
  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  #my_report$name <- sub("^-_root.", "", my_report$name)
  
  splits <- strsplit(my_report$taxLineage, "\\|")
  
  ## for the root nodes, we'll have to add an 'other' link to account for all cladeReads
  root_nodes <- sapply(splits[sapply(splits, length) ==2], function(x) x[2])
  
  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]
  
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n=2))), stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  links$value <- my_report[sel,"cladeReads"]
  
  my_taxRanks <- taxRanks[taxRanks %in% my_report$taxRank]
  taxRank_to_depth <- stats::setNames(seq_along(my_taxRanks)-1, my_taxRanks)
  
  
  nodes <- data.frame(name=my_report$name,
                      depth=taxRank_to_depth[my_report$taxRank],
                      value=my_report$cladeReads,
                      stringsAsFactors=FALSE)
  
  for (node_name in root_nodes) {
    diff_sum_vs_all <- my_report[my_report$name == node_name, "cladeReads"] - sum(links$value[links$source == node_name])
    if (diff_sum_vs_all > 0) {
      nname <- paste("other", sub("^._","",node_name))
      #nname <- node_name
      #links <- rbind(links, data.frame(source=node_name, target=nname, value=diff_sum_vs_all, stringsAsFactors = FALSE))
      #nodes <- rbind(nodes, nname)
    }
  }
  
  names_id = stats::setNames(seq_len(nrow(nodes)) - 1, nodes[,1])
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[links$source != links$target, ]
  
  nodes$name <- sub("^._","", nodes$name)
  links$source_name <- nodes$name[links$source + 1]
  
  if (!is.null(links))
    sankeyD3::sankeyNetwork(
      Links = links,
      Nodes = nodes,
      doubleclickTogglesChildren = TRUE,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      NodeGroup = "name",
      NodePosX = "depth",
      NodeValue = "value",
      dragY = TRUE,
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 15,
      linkGradient = TRUE,
      nodeShadow = TRUE,
      nodeCornerRadius = 5,
      units = "cladeReads",
      fontSize = 12,
      iterations = maxn * 100,
      align = "none",
      highlightChildLinks = TRUE,
      orderByPath = TRUE,
      scaleNodeBreadthsByString = TRUE,
      zoom = zoom,
      ...
    )
}




report_table=read_report2(file(my_report))


report_plot=build_sankey_network(report_table)
saveNetwork(report_plot,OutHtml)
