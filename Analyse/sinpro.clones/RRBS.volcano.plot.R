#'

S4.to.df.diff <- function(S4) {
  df <- data.frame(chr = S4$chr, start = S4$start, end = S4$end, pvalue = S4$pvalue, qvalue = S4$qvalue, meth.diff = S4$meth.diff)
  
  df <- df %>% dplyr::filter(chr == "1" | chr == "2" | chr == "3" | chr == "4" | chr == "5" | chr == "6" | chr == "7" | chr == "8" | chr == "9" | chr == "10" | chr == "11" | chr == "12" | chr == "13" | chr == "14" | chr == "15" | chr == "16" | chr == "17" | chr == "18" | chr == "19" | chr == "20" | chr == "21" | chr == "22" | chr == "X")

Richtung.chr.num <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")

df$chr <- factor(df$chr, levels = Richtung.chr.num)
  
  return(df)
}

gene.name.annotation <- function(df, df.ref) {
  library(data.table)
  
  df.tb <- setDT(df)
  df.ref.tb <- setDT(df.ref)
  
  df.overlap <- as.data.frame(df.tb[df.ref.tb, on = "chr", allow.cartesian = T]) #many-to-many matching
  
  df.overlap_ <- df.overlap %>% dplyr::select(pvalue, qvalue, meth.diff, ENSG)
  #df.overlap_.agg <- aggregate(. ~ ENSG, df.overlap_, mean)
  
  return(df.overlap_)
}

volcano.plotting <- function(S4, df.ref) {
  
  df <- S4.to.df.diff(S4)
  
  df.out <- gene.name.annotation(df, df.ref)
  
  Kolory_volcano <- c("orangered","skyblue2","grey")
names(Kolory_volcano) <- c("up", "down", "no.sig")

  df.sig <- df.out %>% dplyr::mutate(sig = case_when(pvalue < 0.01 & qvalue < 0.01 & meth.diff > 25 ~ "up", pvalue < 0.01 & qvalue < 0.01 & meth.diff < -25 ~ "down", pvalue > 0.01 | qvalue > 0.01 ~ "no.sig")) %>% na.omit()
  
  df.sig$label_de <- NA
df.sig$label_de[df.sig$sig != "no.sig"] <- df.sig$ENSG[df.sig$sig != "no.sig"]

top.list <- df.sig[order(df.sig$pvalue), ] 
top.list.uni <- unique(top.list)
top.list.uni$top <- ""
top.list.uni$top[1:10] <- top.list.uni$ENSG[1:10]
  
  ggplot(df.sig, aes(x = meth.diff, y = -log10(pvalue)))+geom_point(aes(col = sig), alpha = 0.5)+theme_bw()+scale_color_manual(values= Kolory_volcano)+geom_text_repel(data = top.list.uni, aes(label = top), max.overlaps=Inf)+theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12, colour = "black",angle=0,vjust=0.5), axis.title.y = element_text(size=12),axis.text.y = element_text(size = 12, colour = "black"))+xlab("Differential methylation (genomic region)")+ylab("P value (-log10)")
}
