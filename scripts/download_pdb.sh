# You may need to adjust the domain name (rsync.ebi.ac.uk) depending on your location.
# By default the scripts uses the UK web server, suitable for Europe (see
# https://www.wwpdb.org/ftp/pdb-ftp-sites for details).
rsync -rlpt -v -z --delete \
    --ignore-existing \
    rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/pdb/ \
    ./data/pdb
