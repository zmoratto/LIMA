#!/bin/sh

gunzip -c data/RDR_hadley_LOLA.csv.gz > data/RDR_hadley_LOLA.csv
./lidar2img -l data/RDR_hadley_LOLA.csv -i data/AS15-M-1134-hadley.cub --outputImage hadley_alignment.png
echo "Aligned tracks output to hadley_alignment.png."

