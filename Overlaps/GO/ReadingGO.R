dat <- readRDS("RDC.rds")
cat(dat,sep="\n")

#grep -o 'GO:[0-9]*' fileName

# sort fruit.txt | uniq -c | sort -r > temp.txt
# awk ' { t = $1; $1 = $2; $2 = t; print; } ' temp.txt

