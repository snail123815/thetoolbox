#! /bin/bash

TARGETS=(HM84 QL86 BE26 ELS4 PG2)
for t in "${TARGETS[@]}"; do
    rsync -a --info=progress2 Ilis:/vol/local/streptomyces-collection/"$t" ./ --exclude "rawdata"
done

