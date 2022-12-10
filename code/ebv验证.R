tmpppp=merge(x2,dePC_GSE51575,by.x='symbol',by.y='GeneName')
tmpppp=tmpppp[tmpppp$adj.P.Val<0.05,]
tmpppp=tmpppp[tmpppp$symbol %in% m3$symbol,]
tmpppp2=tmpppp[abs(tmpppp$logFC.y)>0.58,]
tmpppp$logFC.y



rnaExpr222222=normalizeBetweenArrays(rnaExpr)
rnaExpr
