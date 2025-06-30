process find_mapped_and_unmapped_regions_per_sampleChrom {
  label 'medium_job'

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      path chrom_sizes_f
    output:
        tuple val(chrom), val(sample_id), path("${sample_id}_unmapped_regions.${chrom}.bed"), path("${sample_id}_mapped_regions.${chrom}.bed")
    script:
    """
    mapped_regions_f="${sample_id}_mapped_regions.${chrom}.bed"
    unmapped_regions_f="${sample_id}_unmapped_regions.${chrom}.bed"

    samtools view -h ${bam} ${chrom} | bedtools bamtobed -i - | bedtools merge > \$mapped_regions_f
    grep -w ${chrom} ${chrom_sizes_f} | bedtools complement -i \$mapped_regions_f -g stdin > \$unmapped_regions_f
    """
}

process acrossSamples_mapped_unmapped_regions_perChr {
  label 'mini_job'

  input:
      tuple val(chrom), val(sample_ids), path(unmapped_beds), path(mapped_beds)
    output:
        tuple val(chrom), path("split_points_${chrom}.bed"), emit: unmapped_bed
        tuple val(chrom), path("mapped_regions_${chrom}.bed"), emit: mapped_bed
    script:
    """
    mapped_beds=(${mapped_beds.join(' ')})
    unmapped_beds=(${unmapped_beds.join(' ')})

    mapped_output_f=mapped_regions_"${chrom}".bed
    unmapped_output_f=split_points_"${chrom}".bed

    count_all_files="\${#mapped_beds[@]}";

    cat "\${mapped_beds[@]}" | sort -k1,1V -k2,2n | bedtools merge -i stdin > \$mapped_output_f;
    bedtools multiinter -i \${unmapped_beds[@]} | awk -v count_all_files="\${count_all_files}" '\$4==count_all_files' | cut -f1,2,3 > \$unmapped_output_f;
    """
}


process suggest_splits_binarySearch {
    label 'mini_job'
    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"

    input:
      tuple val(chrom), val(sample_ids), path(bams), path(bais), path(unmapped_regions_bed)
      val chunks
      path chrom_sizes_f
    output:
      tuple val(chrom), path("suggested_splits_onebased_coords.${chrom}.bed"), path("suggested_splits_onebased_coords.${chrom}.list"), path("suggested_splits.${chrom}.bed")
    script:
    """
    printf "${bams.join('\n')}" > bams.txt
    python ${baseDir}/dev/split_chr.py  -c ${chunks} -b bams.txt -r "${chrom}" -s ${unmapped_regions_bed} -z ${chrom_sizes_f} -o suggested_splits."${chrom}".bed
    awk -F "\t" '{print \$1"\t"(\$2+1)"\t"\$3}' suggested_splits."${chrom}".bed > suggested_splits_onebased_coords."${chrom}".bed
    awk -F "\t" '{print \$1":"\$2"-"\$3}' suggested_splits_onebased_coords."${chrom}".bed  > suggested_splits_onebased_coords."${chrom}".list
    rm bams.txt
    """
}


process split_bams_perChunk {
  label 'micro_multithread_job'
  input:
    tuple val(chrom), val(sample_ids), path(bams), path(bais), val(formattedRegion), val(programmaticRegion)
    val(suffix)
  output:
    tuple val(chrom), val(sample_ids), path(outputBams), path(outputBais), val(formattedRegion), val(programmaticRegion),path("${programmaticRegion}_total_count.csv")
  script:

  outputBams = sample_ids.collect { sample_id -> sample_id+".${programmaticRegion}.${suffix}.bam" }
  outputBais = sample_ids.collect { sample_id -> sample_id+".${programmaticRegion}.${suffix}.bam.bai" }

  """
  #setting local bash variables
  bams=(${bams.join(' ')})
  sample_ids=(${sample_ids.join(' ')})
  programmaticRegion="${programmaticRegion}"
  formattedRegion="${formattedRegion}"
  chrom="${chrom}"
  suffix="${suffix}"

  total_count=0
  num_samples="\${#sample_ids[@]}"
  for i in \$(seq 0 \$((\$num_samples-1)) ); do
    sample_id="\${sample_ids[\$i]}"
    input_bam="\${bams[\$i]}"
    output_bam="\${sample_id}.\${programmaticRegion}.\${suffix}.bam"

    samtools view -@ ${task.cpus} -h "\${input_bam}" "\${formattedRegion}" |\
    awk -v chrom="\${chrom}" '{ if (\$1=="@SQ") {if(\$2=="SN:"chrom) {print \$0}} else{print \$0} }' |\
    samtools view -@ ${task.cpus} -h -bo "\${output_bam}" -
    samtools index -@ ${task.cpus} "\${output_bam}";
    count=\$(samtools view -@ ${task.cpus} -c \$output_bam)
    total_count=\$((\$total_count+\$count))
  done;
  echo "\${formattedRegion},\${total_count}" > "\${programmaticRegion}_total_count.csv"
  """
}
