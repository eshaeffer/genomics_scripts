remove(list=objects())
#chooseCRANmirror()    # Uncomment and run this command and choose a mirror
                       # if any of install commands do not work
if (!require("readr")) {
  BiocManager::install("readr")
}
if (!require("mygene")) {
  BiocManager::install("mygene")
}

library(mygene)

(tdate <- format(Sys.time(), "%b_%d_%Y"))


my_data <- read_csv("/Users/eliasshaeffer/Desktop/somatic_gatk_low/genes.txt");
my_data$symbol
colnames(my_data)
#length(gene_info$response$query)

length(my_data$symbol)

df_resOrdered <- as.data.frame(my_data)
df_resOrdered

#call mygenes:
gene_info <- queryMany(my_data$symbol, species=9606, scopes="symbol", fields=c("summary", "interpro.desc"), returnall=TRUE)
(gene_info$response$query[6])    # query is mygenes equivalent of our "symbol"
class(gene_info$response$interpro[6])
df_gene_info <- as.data.frame(gene_info$response)
#colnames(df_resOrdered)
(cn_resOrdered <- colnames(df_resOrdered))
#class((df_gene_info$))
(c <- paste( unlist(df_gene_info$interpro[1]), collapse='; '))

#Merge the various interpro values into one:
for (i in 1:(length(df_gene_info$interpro))) {
   df_gene_info$interpro[i] = paste( unlist(df_gene_info$interpro[i]), collapse='; ')
 }
#colnames(df_resOrdered)
#colnames(df_gene_info)
#df_gene_info$query
#df_resOrdered$symbol

df_gene_info_three_columns = df_gene_info[,c("query","summary","interpro")]
mergedResOrdered2 <- 
  merge(df_resOrdered, df_gene_info_three_columns, 
                by.x = c("symbol"), by.y = c("query"), 
                y.all=TRUE, sort = FALSE )
#cols <- c("",colnames(mergedResOrdered2)) 
write.table(df_gene_info_three_columns, "/Users/eliasshaeffer/Desktop/somatic_gatk_low/gene_summaries.txt", sep='\t')
