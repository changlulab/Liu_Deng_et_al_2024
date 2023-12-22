for file in *-1_R1_Full_Normalized.bedGraph
do
base="${file%-1_R1_Full_Normalized.bedGraph}"

cut -d$'\t' -f4 $PWD/$base"-1_R1_Full_Normalized.bedGraph" > $PWD/$base"_1_Signal.bed"
cut -d$'\t' -f4 $PWD/$base"-2_R1_Full_Normalized.bedGraph" > $PWD/$base"_2_Signal.bed"
paste $PWD/$base"_1_Signal.bed" $PWD/$base"_2_Signal.bed" > $PWD/$base"_combined_Signal.bed"

awk -F$'\t' 'BEGIN {OFS=FS}{$3=($1+$2)/2; print}' $PWD/$base"_combined_Signal.bed" > $PWD/$base"_average_Signal.bed"

cut -d$'\t' -f3 $PWD/$base"_average_Signal.bed" > $PWD/$base"_Average.bed"
paste ~/Data/Promotors/mm10/GenomeWindowed_sort.bed $PWD/$base"_Average.bed" > $PWD/$base"_Average.bedGraph"

bedGraphToBigWig $PWD/$base"_Average.bedGraph" ~/Data/Promotors/mm10/Genome.bed $PWD/$base"_Average.bw"

rm $PWD/$base"_1_Signal.bed" $PWD/$base"_2_Signal.bed" $PWD/$base"_combined_Signal.bed" $PWD/$base"_average_Signal.bed" $PWD/$base"_Average.bed" $PWD/$base"_Average.bedGraph"

done
