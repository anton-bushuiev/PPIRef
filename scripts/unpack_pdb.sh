indir="data/pdb"
counter=0; total=$(find $indir -type f -name "*.ent.gz" | wc -l); \
find $indir -type f -name "*.ent.gz" -print0 | while IFS= read -r -d '' file; \
do \
    stem=$(basename "$file"); \
    directory=$(dirname "$file"); \
    counter=$((counter+1)); \
    gzip -dc "$file" > "${directory}/${stem:3:4}.pdb" && rm "$file" && \
    echo -ne "Unpacking: $counter/$total\r"; \
done
