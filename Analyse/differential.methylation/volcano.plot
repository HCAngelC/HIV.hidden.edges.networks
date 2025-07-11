S4.to.df.diff <- function(S4) {
  df <- data.frame(chr = S4$chr, start = S4$start, end = S4$end, pvalue = S4$pvalue, qvalue = S4$qvalue, meth.diff = S4$meth.diff)
  
  df <- df %>% dplyr::filter(chr == "1" | chr == "2" | chr == "3" | chr == "4" | chr == "5" | chr == "6" | chr == "7" | chr == "8" | chr == "9" | chr == "10" | chr == "11" | chr == "12" | chr == "13" | chr == "14" | chr == "15" | chr == "16" | chr == "17" | chr == "18" | chr == "19" | chr == "20" | chr == "21" | chr == "22" | chr == "X")

Richtung.chr.num <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")

df$chr <- factor(df$chr, levels = Richtung.chr.num)
  
  return(df)
}

volcano.plotting <- function(S4) {
  
  df <- S4.to.df.diff(S4)
  
  Kolory_volcano <- c("orangered","skyblue2","grey")
names(Kolory_volcano) <- c("up", "down", "no.sig")

  df.sig <- df %>% dplyr::mutate(sig = case_when(pvalue < 0.01 & qvalue < 0.01 & meth.diff > 25 ~ "up", pvalue < 0.01 & qvalue < 0.01 & meth.diff < -25 ~ "down", pvalue > 0.01 | qvalue > 0.01 ~ "no.sig")) %>% na.omit()
  
  ggplot(df.sig, aes(x = meth.diff, y = -log10(pvalue)))+geom_point(aes(col = sig), alpha = 0.5)+theme_bw()+scale_color_manual(values= Kolory_volcano)+theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12, colour = "black",angle=0,vjust=0.5), axis.title.y = element_text(size=12),axis.text.y = element_text(size = 12, colour = "black"))+xlab("Differential methylation (region)")+ylab("P value (-log10)")
}


#Execution
volcano.plotting(f814.filter.merge.tiles.diff)
