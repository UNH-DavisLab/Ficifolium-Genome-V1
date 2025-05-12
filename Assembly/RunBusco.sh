#!/bin/bash

docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.3.2_cv1 \
busco -i ordered_new_ficifolium_corrected.fa -o UnmaskedCorrected2024embryophyta_odb10 -l embryophyta_odb10/ \
-m genome -c 22
