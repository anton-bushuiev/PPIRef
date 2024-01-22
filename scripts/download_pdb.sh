# https://www.wwpdb.org/ftp/pdb-ftp-sites
rsync -rlpt -v -z --delete \
    --ignore-existing \
    rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/pdb/ \
    ./data/pdb
