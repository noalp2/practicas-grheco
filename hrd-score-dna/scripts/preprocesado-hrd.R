#' Preprocessing seqz files
#'
#' Preprocesses seqz files
#' @param seg A segmentation file
#' @param ploidy0 ploidy, optional
#' @return Output is
preprocess.sequenza<-function(seg, ploidy0=NULL, chr.in.names=TRUE, outputdir=NULL){
  if (is.null(ploidy0)){
    ploidy01 = seq(1, 5.5, 0.1)
  } else {
    ploidy01= seq(ploidy0-0.5,ploidy0+0.5,0.1)
  }
  
  if (is.null(outputdir)){
    outputdir = getwd()
  }
  
  run_name<-gsub(".*/","",gsub("_small.seqz","",gsub("gz","",seg)))
  if(chr.in.names){
    extract<-sequenza.extract(seg, chromosome.list=paste('chr',c(1:24),sep=''))
  } else {
    extract<-sequenza.extract(seg, chromosome.list=c(1:24))
  }
  extract.fit<-sequenza::sequenza.fit(extract, N.ratio.filter = 10, N.BAF.filter = 1, segment.filter = 3e6, mufreq.treshold = 0.10, ratio.priority = FALSE,ploidy=ploidy01, mc.cores = 1)
  #  sequenza.results(extract, extract.fit, out.dir = getwd(),sample.id =run_name)
  
  seg.tab <- do.call(rbind, extract$segments[extract$chromosomes])
  seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06
  cint <- get.ci(extract.fit)
  cellularity <- cint$max.cellularity
  ploidy <- cint$max.ploidy
  avg.depth.ratio <- mean(extract$gc$adj[, 2])
  info_seg<-c(cellularity,ploidy,avg.depth.ratio)
  names(info_seg)<-c("cellularity","ploidy","avg.depth.ratio")
  write.table(t(info_seg),paste0(outputdir,"/",run_name,"_info_seg.txt"),sep="\t",row.names=F)
  allele.cn <- sequenza:::baf.bayes(Bf = seg.tab$Bf, CNt.max = 20, depth.ratio = seg.tab$depth.ratio, avg.depth.ratio = 1,
                                    cellularity = cint$max.cellularity, ploidy = cint$max.ploidy,
                                    sd.ratio = seg.tab$sd.ratio, weight.ratio = seg.len, sd.Bf = seg.tab$sd.BAF,
                                    weight.Bf = 1, ratio.priority = FALSE, CNn = 2)
  seg.tab$CN <- allele.cn[,1]
  allele.cn <- as.data.table(allele.cn)
  #Making imput file
  seg <- data.frame(SampleID = as.character(run_name), Chromosome = seg.tab$chromosome, Start_position = seg.tab$start.pos,
                    End_position = seg.tab$end.pos, Nprobes = 1, total_cn = allele.cn$CNt, A_cn = allele.cn$B,
                    B_cn = allele.cn$A, ploidy = ploidy)
  seg$contamination <- 1
  seg<-seg[!is.na(seg$A_cn),]
  seg<-seg[!is.na(seg$B_cn),]
  return(seg)
}
sequenza.fit <- function (sequenza.extract, female = TRUE, N.ratio.filter = 10, 
    N.BAF.filter = 1, segment.filter = 3e+06, mufreq.treshold = 0.1, 
    XY = c(X = "X", Y = "Y"), cellularity = seq(0.1, 1, 0.01), 
    ploidy = seq(1, 7, 0.1), ratio.priority = FALSE, method = "baf", 
    priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24, 
    mc.cores = getOption("mc.cores", 2L)) 
{
    if (is.null(chromosome.list)) {
        mut.all <- do.call(rbind, sequenza.extract$mutations)
        mut.all <- na.exclude(mut.all)
        segs.all <- do.call(rbind, sequenza.extract$segments)
    }
    else {
        mut.all <- do.call(rbind, sequenza.extract$mutations[chromosome.list])
        mut.all <- na.exclude(mut.all)
        segs.all <- do.call(rbind, sequenza.extract$segments[chromosome.list])
    }
    segs.len <- segs.all$end.pos - segs.all$start.pos
    avg.depth.ratio <- sequenza.extract$avg.depth.ratio
    if (method == "baf") {
        avg.sd.ratio <- sum(segs.all$sd.ratio * segs.all$N.ratio, 
            na.rm = TRUE)/sum(segs.all$N.ratio, na.rm = TRUE)
        avg.sd.Bf <- sum(segs.all$sd.BAF * segs.all$N.BAF, na.rm = TRUE)/sum(segs.all$N.BAF, 
            na.rm = TRUE)
        segs.all$sd.BAF[segs.all$sd.BAF == 0] <- max(segs.all$sd.BAF, 
            na.rm = TRUE)
        segs.all$sd.ratio[segs.all$sd.ratio == 0] <- max(segs.all$sd.ratio, 
            na.rm = TRUE)
        segs.filt <- segs.all$N.ratio > N.ratio.filter & segs.all$N.BAF > 
            N.BAF.filter
        segs.filt <- segs.len >= segment.filter & segs.filt
        if (female) {
            segs.is.xy <- segs.all$chromosome == XY["Y"]
        }
        else {
            segs.is.xy <- segs.all$chromosome %in% XY
        }
        filt.test <- segs.filt & !segs.is.xy
        seg.test <- segs.all[filt.test, ]
        seg.len.mb <- segs.len[filt.test]/1e+06
        baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio, 
            sd.ratio = seg.test$sd.ratio, weight.ratio = seg.len.mb, 
            sd.Bf = seg.test$sd.BAF, weight.Bf = seg.len.mb, 
            avg.depth.ratio = avg.depth.ratio, cellularity = cellularity, 
            ploidy = ploidy, priors.table = priors.table, mc.cores = mc.cores, 
            ratio.priority = ratio.priority)
    }
    else if (method == "mufreq") {
        mut.filt <- mut.all$F >= mufreq.treshold
        if (female) {
            mut.is.xy <- mut.all$chromosome == XY["Y"]
        }
        else {
            mut.is.xy <- mut.all$chromosome %in% XY
        }
        filt.test <- mut.filt & !mut.is.xy
        mut.test <- mut.all[filt.test, ]
        w.mufreq <- round(mut.test$good.reads, 0)
        mufreq.model.fit(mufreq = mut.test$F, depth.ratio = mut.test$adjusted.ratio, 
            weight.ratio = 2 * w.mufreq, weight.mufreq = w.mufreq, 
            avg.depth.ratio = avg.depth.ratio, cellularity = cellularity, 
            ploidy = ploidy, priors.table = priors.table, mc.cores = mc.cores)
    }
    else {
        stop("The only available methods are \"baf\" and \"mufreq\"")
    }
}
gc.sample.stats<-function(file, col_types = "c--dd----d----",
    buffer = 33554432, verbose = TRUE) {

    con <- gzfile(file, "rb")

    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, col_types) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, col_names = FALSE,
            skip = 0, n_max = Inf, progress = FALSE)
        u_chr <- unique(x[, 1])
        n_chr <- table(x[, 1])
        gc1 <- lapply(split(x[, 2], x[, 4]), table)
        gc2 <- lapply(split(x[, 3], x[, 4]), table)
        if (verbose){
            message(".", appendLF = FALSE)
        }
        list(unique = u_chr, lines = n_chr, gc_nor = gc1, gc_tum = gc2)
    }
    if (verbose){
        message("Collecting GC information ", appendLF = FALSE)
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, col_types = col_types,
        CH.MAX.SIZE = buffer)
    close(con)
    if (verbose){
        message(" done\n")
    }
    unfold_gc(res, stats = TRUE)
}

unfold_gc <- function(x, stats = TRUE) {
    gc_norm <- get_gc(x[, "gc_nor"])
    gc_tum <- get_gc(x[, "gc_tum"])
    if (stats) {
        ord_chrom <- unique(Reduce("c", Reduce("c", x[, "unique"])))
        stats_chrom <- Reduce("c", x[, "lines"])
        stats_chrom <- sapply(splash_table(x[, "lines"]), sum)
        stats_chrom <- stats_chrom[ord_chrom]
        stats_start <- cumsum(c(1, stats_chrom[-length(stats_chrom)]))
        stats_end   <- stats_start + stats_chrom - 1
        stats_chrom <- data.frame(chr = ord_chrom, n_lines = stats_chrom,
            start = stats_start, end = stats_end)

        list(file.metrics = stats_chrom, normal = gc_norm, tumor = gc_tum)
    } else {
        list(normal = gc_norm, tumor = gc_tum)
    }
}

splash_table <- function(lis_obj){
    lis_obj <- Reduce("c", lis_obj)
    split(lis_obj, names(lis_obj))
}

get_gc <- function(gc_col) {
    sort_char <- function(x) {
        as.character(sort(as.numeric(x)))
    }
    all_depths <- splash_table(gc_col)
    all_depths <- lapply(all_depths, FUN = function(x) {
        sapply(splash_table(x), sum)
    })
    names_gc <- sort_char(names(all_depths))
    all_depths <- all_depths[names_gc]
    names_depths <- sort_char(unique(Reduce("c", lapply(all_depths, names))))
    n <- do.call(rbind, lapply(all_depths, FUN = function(x, names_depths) {
            res <- x[names_depths]
            names(res) <- names_depths
            res
        },
        names_depths = names_depths))
    n[is.na(n)] <- 0
    list(gc = as.numeric(names_gc), depth = as.numeric(names_depths), n = n)
}

median_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.median(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

mean_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.mean(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

depths_gc <- function(depth_n, depth_t, gc) {
    gc_nor <- lapply(split(depth_n, gc), table)
    gc_tum <- lapply(split(depth_t, gc), table)
    list(gc_nor = gc_nor, gc_tum = gc_tum)
}
