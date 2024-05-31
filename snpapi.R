getSNPaddr = function(rsid="rs6060535") {
 ans = httr::GET(
  sprintf("https://clinicaltables.nlm.nih.gov/api/snps/v3/search?terms=%s",
     rsid))
 #ans = httr::content(ans)
 # check that a valid result returned
 # pick out the exact match (or figure out how to use the API to get the exact match)
 # if it is not in the space we want, use liftOver to get it there
 ans
 # liftOver is in rtracklayer, get the chain files from UCSC
 # 0-based indexing vs 1-based indexing -- does API come back with 1-based while dbSNP uses 0 based?
}
