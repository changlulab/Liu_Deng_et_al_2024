cat *.bed > M12_K27me3_glia_combined.bed
sort -k1,1 -k2,2n M12_K27me3_glia_combined.bed > M12_K27me3_glia_sorted.bed
bedtools merge -i M12_K27me3_glia_sorted.bed > M12_K27me3_glia_all_peaks.bed
