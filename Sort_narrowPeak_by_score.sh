for f in *_peaks.narrowPeak
do

base=`basename "$f" _peaks.narrowPeak`
sort -r -g -k 5,5 < $base"_peaks.narrowPeak" > $base"_temp.out"
lncnt=$(wc -l < $base"_temp.out")
percent_linecount_infloat=$(echo "$lncnt*.3" | bc)
float2Int=$(printf %.0f "$percent_linecount_infloat")
head -"$float2Int" $base"_temp.out" > $base"_top03.narrowPeak"
rm $PWD/$base"_temp.out"
done
