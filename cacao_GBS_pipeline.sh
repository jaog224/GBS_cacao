mkdir loc_realigned
mkdir recalibrationç
mkdir final_bam
mkdir calls
mkdir callability
for((i=#;i<=#;i++)); do
        mySM=$(cat list_samples |sed -n ''$i'p' |awk '{print $1}')
        myID=$(cat list_samples |sed -n ''$i'p' |awk '{print $2}')
        myfile1=$(cat list_samples |sed -n ''$i'p' |awk '{print $3}')
        myfile2=$(cat list_samples |sed -n ''$i'p' |awk '{print $4}')
        echo "I have a sample "$mySM" with a unique ID "$myID
        echo "my file 1 its in "$myfile1
        echo "my file 2 its in "$myfile2
        bwa mem -B 6 -R '@RG\tID:'$myID'\tSM:'$mySM'\tPL:Illumina' /path_to_reference_genome/ $myfile1 $myfile2 > $myID.sam
        samtools view -ubhSt /path_to_reference_genome/ $myID.sam -o $myID.bam
        java -Xmx6G -jar /path_to_software/picard.jar CleanSam INPUT=$myID.bam OUTPUT=$myID.clean.bam VALIDATION_STRINGENCY=SILENT
        java -Xmx6G -jar /path_to_software/picard.jar FixMateInformation INPUT=$myID.clean.bam OUTPUT=$myID.clean.fix.bam VALIDATION_STRINGENCY=SILENT
        java -Xmx6G -jar /path_to_software/picard.jar ValidateSamFile INPUT=$myID.clean.fix.bam OUTPUT=$myID.validation VALIDATION_STRINGENCY=LENIENT
        java -Xmx6G -jar /path_to_software/picard.jar CollectAlignmentSummaryMetrics INPUT=$myID.clean.fix.bam OUTPUT=$myID.metrics VALIDATION_STRINGENCY=LENIENT
        java -Xmx4G -jar /path_to_software/picard.jar SortSam INPUT=$myID.clean.fix.bam OUTPUT=$myID.F.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
        /path_to_software/samtools index $myID.F.bam
        /path_to_software/samtools idxstats $myID.F.bam > $myID.idxstats
        /path_to_software/bamtools stats -in $myID.F.bam > $myID.STATS
        rm $myID.bam $myID.clean.bam $myID.clean.fix.bam
		java -Xmx6g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -T RealignerTargetCreator -I $myID.F.bam -o loc_realigned/$myID.forIndelRealigner.intervals
        java -Xmx6g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -I $myID.F.bam -T IndelRealigner -targetIntervals loc_realigned/$myID.forIndelRealigner.intervals -o loc_realigned/$myID.realigned.bam
		/path_to_software/samtools index loc_realigned/$myID.realigned.bam
		java -Xmx6g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -T BaseRecalibrator -knownSites /path_to_known_sites_vcf/ -I loc_realigned/$myID.realigned.bam -o recalibration/$myID.recal_data.grp --covariate ContextCovariate --covariate RepeatLengthCovariate --covariate CycleCovariate --validation_strictness SILENT 
		java -Xmx6g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -T PrintReads -I loc_realigned/$myID.realigned.bam -BQSR recalibration/$myID.recal_data.grp -o final_bam/$myID.recal.bam --validation_strictness SILENT --allow_potentially_misencoded_quality_scores
		java -Xmx12g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -T UnifiedGenotyper -I /path_to_recal.bam_files/ -o calls/multiple_samples.vcf --output_mode EMIT_ALL_SITES --heterozygosity 0.006 --indel_heterozygosity 1.25E-3 -dcov 200
	 	java -Xmx6g -jar /path_to_software/GenomeAnalysisTK.jar -R /path_to_reference_genome/  -T CallableLoci -I annotate/$myID.annotated.bam -summary annotate/$myID.summaries.txt -o callability/$myID.callable_status.bed
		/path_to_software/vcftools --vcf /vcf_from_UnifiedGenotyper/ --max-alleles 2 -min-alleles 2 --max-missing 1 --out final.vcf --recode
done
