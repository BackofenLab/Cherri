

library(tidyverse)
library(readxl)
library(writexl)

data <- read_xlsx("benchmark-data.xlsx", na = "NA", sheet="benchmark-data") |> 
  separate(col=qRegion, sep = "\\.\\.", into=c("qStart","qEnd"), convert=T) |> 
  separate(col=tRegion, sep = "\\.\\.", into=c("tStart","tEnd"), convert=T) 


results <- list()
  
# IntaRNA results
# results[["IntaRNAhelix"]] <- read_delim("IntaRNAhelix.csv", delim=";")
results[["IntaRNA"]] <- read_delim("IntaRNA.csv", delim=";")

# RIblast results
results[["RIblast"]] <- 
  read_csv("RIblast.out", skip=2) |> 
  rename(id1 = `Target name`, id2= `Query name`) |> 
  filter( str_extract(id1,"\\d+") == str_extract(id2,"\\d+")) |> 
  # reduce to mfe prediction
  group_by(id1,id2) |> 
  slice_min(`Hybridization Energy`, with_ties=F) |> 
  ungroup() |> 
  select( -Id, -ends_with("Length"), -ends_with("Energy")) |> 
  mutate( 
    id= str_c("ug.",str_extract(id1,"\\d+")),
    BasePair = str_extract(BasePair,"\\d.*\\d")
  ) |> 
  separate(BasePair, sep=":", into=c("qRegion","tRegion")) |> 
  separate(qRegion,sep="-", into=c("start2","end2"), convert = T) |> 
  separate(tRegion,sep="-", into=c("start1","end1"), convert = T) 

# RIsearch2 results
results[["RIsearch2"]] <- 
  read_tsv("RIsearch2.out", col_names = c("id2","start2","end2","id1","start1","end1","strand","energy")) |> 
  filter( str_extract(id1,"\\d+") == str_extract(id2,"\\d+"), 
  # drop negative strand rris
          strand == "+") |> 
  # keep only mfe entry
  group_by(id1,id2) |> 
  slice_min(energy,n=1, with_ties = F) |> 
  ungroup()
  

#######################################################
# check if the prediction overlaps with the experimental data
overlap <-
  bind_rows(results, .id="tool") |> 
  mutate(id = str_c("ug.",str_extract(id1,"\\d+"))) |> 
  left_join(data, by = "id") |> 
  mutate( 
    qOverlap = (start2 <= qStart & qStart <= end2)
    | (qStart <= start2 & start2 <= qEnd),
    tOverlap = (start1 <= tStart & tStart <= end1)
    | (tStart <= start1 & start1 <= tEnd),
    overlap = (is.na(qOverlap | qOverlap) & (is.na(tOverlap)|tOverlap))
  ) |> 
   # view() |> 
  select( tool, id, starts_with("start"), starts_with("end"), ends_with("verlap"), qType) |> 
  distinct()

write_csv(overlap, "overlap.csv")


# stats on overlaps per tool
overlap |> 
  mutate( n = unique(id) |> length()) |> 
  group_by(qType) |> 
  mutate( nQtype = unique(id) |> length() ) |> 
  ungroup() |> 
  group_by(tool) |> 
  summarize(
    pred = n(),
    TP.rate = sum(overlap)/n(),
    FP.rate = (n()-sum(overlap))/n()
  ) 

overlap |> 
  group_by(tool, overlap) |> 
  count() |> 
  ungroup() |> 
  group_by(tool) |> 
  mutate(nTool = sum(n), a = ifelse(overlap,0.4,0.2)) |> 
  unite(t, tool, overlap, remove = F) |> 
  ungroup() |> 
  arrange(overlap,tool) |> 
  mutate(t=fct_inorder(t,ordered=T), height=cumsum(n)) -> tmp
  
library(ggnewscale)
tmp |> 
  ggplot(aes(x="", y=n, group=overlap))+
  theme_light()+
  # geom_bar(aes(fill=overlap, group=overlap),alpha=1,width = 1, stat = "identity", col="white",size=4) +
  # geom_bar(aes(fill=tool,group=overlap, col=overlap),
  #          width = 1, stat = "identity", alpha=0.7) +
  geom_col(mapping=aes(fill=overlap))+
  scale_fill_manual(values = c("red","lightblue"))+
  labs(
    fill="correct\nprediction"
  ) +
  new_scale("fill") + 
  geom_col(alpha=0.5, mapping=aes(fill=tool))+
  scale_fill_grey(start = 0.2, end = .8)+
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5), col="black") +
  scale_y_continuous(breaks=cumsum(tmp$n) - tmp$n / 2, labels= tmp$tool) +
  coord_polar("y", start=0) 

ggsave("overlap.png",width=5,height=5)

#######################################################
# get genome positions
input_data <- 
  data |> 
  select(-c(doi, source, qType, ends_with("Notes"), organism)) |> 
  transmute(
    id,
    source = "literature",
    correct = T,
    chrom1 = tHg38Chr,
    start1 = tHg38Start + tStart -1,
    stop1  = tHg38Start + tEnd -1,
    strand1 = tHg38Strand,
    chrom2 = qHg38Chr,
    start2 = qHg38Start + qStart -1,
    stop2  = qHg38Start + qEnd -1,
    strand2 = qHg38Strand
  ) |> 
  drop_na()

input_predictions <- 
  overlap |> 
  left_join(data,by="id") |> 
  transmute(
    id,
    source = tool,
    correct = overlap,
    chrom1 = tHg38Chr,
    start1 = tHg38Start + start1 -1,
    stop1  = tHg38Start + end1 -1,
    strand1 = tHg38Strand,
    chrom2 = qHg38Chr,
    start2 = qHg38Start + start2 -1,
    stop2  = qHg38Start + end2 -1,
    strand2 = qHg38Strand
  )

cherri_input <- 
  bind_rows( input_data, input_predictions) |> 
  write_csv("cherri_input.csv") 


cherri_input |> 
  summarize( n=n(), TP = sum(correct), FP = n-TP, FP.rate = FP/n)
