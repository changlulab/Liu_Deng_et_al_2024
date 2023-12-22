for f in *.EPIC2.bed
do

base=`basename "$f" .EPIC2.bed`
sort -r -g -k 5,5 < $base".EPIC2.bed" > $base"_temp.out"
lncnt=$(wc -l < $base"_temp.out")
percent_linecount_infloat=$(echo "$lncnt*.3" | bc)
float2Int=$(printf %.0f "$percent_linecount_infloat")
head -"$float2Int" $base"_temp.out" > $base"_top03.EPIC2.bed"
rm $PWD/$base"_temp.out"
done
