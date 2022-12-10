library("GEOmetadb")
## ----citation-----------------------------------------------------------------
citation("GEOmetadb")

## ----include=FALSE-----------------------------------
library(knitr)
options(width=55)
opts_chunk$set(cache=TRUE,message=FALSE,warning=FALSE)

## ----------------------------------------------------
library(GEOmetadb)

## ----------------------------------------------------
if( !file.exists("GEOmetadb.sqlite") ) {
  demo_sqlfile <- getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz", type = "demo")
} else {
  demo_sqlfile <- "GEOmetadb.sqlite"
}

## Download the full GEOmetadb database:
# geometadbfile <- getSQLiteFile()


## ----------------------------------------------------
file.info(demo_sqlfile)

## ----------------------------------------------------
con <- dbConnect(SQLite(), demo_sqlfile)
dbDisconnect(con)

## ----------------------------------------------------
con <- dbConnect(SQLite(), demo_sqlfile)

## ----------------------------------------------------
geo_tables <- dbListTables(con)
geo_tables

## ----------------------------------------------------
dbListFields(con,'gse')

## ----------------------------------------------------
dbGetQuery(con,'PRAGMA TABLE_INFO(gpl)')

## ----------------------------------------------------
rs <- dbGetQuery(con,'select * from gse limit 5')
rs[,1:7]

## ----------------------------------------------------
rs <- dbGetQuery(con,paste("select gse,title from gse where",
                           "contributor like '%Sean%Davis%'",sep=" "))
rs

## ----------------------------------------------------
rs <- dbGetQuery(con,paste("select gsm,supplementary_file",
                           "from gsm where gpl='GPL96'",
                           "and supplementary_file like '%CEL.gz'"))
dim(rs)

## ----------------------------------------------------
rs <- dbGetQuery(con,paste("select gpl.bioc_package,gsm.gpl,",
                           "gsm,gsm.supplementary_file",
                           "from gsm join gpl on gsm.gpl=gpl.gpl",
                           "where gpl.manufacturer='Affymetrix'",
                           "and gsm.supplementary_file like '%CEL.gz' "))
rs[1:5,]

## ----------------------------------------------------
getTableCounts <- function(tableName,conn) {
  sql <- sprintf("select count(*) from %s",tableName)
  return(dbGetQuery(conn,sql)[1,1])
}
do.call(rbind,sapply(geo_tables,getTableCounts,con,simplify=FALSE))

## ----------------------------------------------------
conversion <- geoConvert('GPL96', out_type = c("gse", "gpl", "gsm", "gds", "smatrix"), 
                         sqlite_db_name = demo_sqlfile)

## ----------------------------------------------------
lapply(conversion, dim)
conversion$gse[1:5,]
conversion$gsm[1:5,]
conversion$gds[1:5,]
conversion$sMatrix[1:5,]

## ----------------------------------------------------
getBiocPlatformMap(con, bioc=c('hgu133a','hgu95av2'))

## ----------------------------------------------------
sql <- paste("SELECT DISTINCT gse.title,gse.gse",
             "FROM",
             "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
             "  JOIN gse ON gse_gsm.gse=gse.gse",
             "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gsm.molecule_ch1 like '%total RNA%' AND",
             "  gse.title LIKE '%breast cancer%' AND",
             "  gpl.organism LIKE '%Homo sapiens%'",sep=" ")
rs <- dbGetQuery(con,sql)
dim(rs)
print(rs[1:5,],right=FALSE)

## ----fig.width=10------------------------------------
library(tm)
library(wordcloud)
gseTitles = dbGetQuery(con,"select title from gse")
corp = VCorpus(VectorSource(gseTitles))
corp <- tm_map(corp, removePunctuation)
corp <- tm_map(corp, content_transformer(tolower))
corp <- tm_map(corp, removeNumbers)
corp <- tm_map(corp, function(x)removeWords(x,stopwords()))
term.matrix <- TermDocumentMatrix(corp)
term.matrix <- as.matrix(term.matrix)
v <- sort(rowSums(term.matrix),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
n = 100
wordcloud(d[1:n,]$word,d[1:n,]$freq)

## ----------------------------------------------------
library(dplyr)
db = src_sqlite(demo_sqlfile)
gse = tbl(db,'gse')
filter(gse,gse=='GSE2553')

## ----------------------------------------------------
dbDisconnect(con)

## ---- eval=TRUE--------------------------------------
# file.remove('GEOmetadb.sqlite')

## ----------------------------------------------------
sessionInfo()



?
  dbGetQuery 
