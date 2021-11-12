cat *.chr1.tad.bed | cut -f 1-2 |sort -k1 -k2 -n |uniq > whole.samples.left.borders.info.txt
python get.consensus.TAD.py whole.samples.left.borders.info.txt ./ samples.txt whole.consensus.tad.chr1.bed 1
bedtools merge -i whole.consensus.tad.chr1.bed > whole.consensus.tad.chr1.merge.bed
python get.complement.tad.py whole.consensus.tad.chr1.merge.bed chr1.size.txt whole.consensus.tad.chr1.merge.true.bed
