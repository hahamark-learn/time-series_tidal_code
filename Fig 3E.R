# a minimal seq track
s0 <- tibble::tibble(
  seq_id = c("a", "b"),
  length = c(600, 550)
)

# a minimal gene track
g0 <- tibble::tibble(
  seq_id = c("a", "a", "b"),
  start = c(50, 350, 80),
  end = c(250, 500, 450)
)

# a simple link track
l0 <- tibble::tibble(
  seq_id = c("a", "a"),
  start = c(50, 400),
  end = c(250, 480),
  seq_id2 = c("b", "b"),
  start2 = c(80, 350),
  end2 = c(300, 430)
)

#按照上面格式修改文件得到原始数据

genemap_length <- read.table("genemap_length",header = T)

genemap_genes <- read.table("genemap_genes",header = T)

genelinks <- read.table("genemap_links",header = T)
gggenomes(genes=genemap_genes, seqs=genemap_length,links = genelinks)+
  geom_seq()+
  geom_seq_label()+
  geom_gene()+
  geom_link()



p <- gggenomes(genes=g0, seqs=s0, links=l0)
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene() +        # draw genes as arrow
  geom_link()          # draw some connections between syntenic regions


setwd("~/Desktop/article3/metagenome/hgt/emales/")
seqs <- read_fai("emales.fna")
emale_seqs  <- read_fai("genemap.fna") 
