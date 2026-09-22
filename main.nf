#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
        Pasteurella multocida LPS analysis pipeline
========================================================================================
 #### Documentation
 #https://github.com/julianzaugg/LPS_typing_Illumina
 #### Authors
 Valentine Murigneux <v.murigneux@uq.edu.au>
 Julian Zaugg <j.zaugg@uq.edu.au>
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	Pasteurella multocida LPS analysis pipeline v${workflow.manifest.version}
	=========================================
	Usage:
	nextflow main.nf --samplesheet /path/to/samplesheet.csv --outdir /path/to/outdir/

	Required arguments:
		--outdir				Path to the output directory to be created
		--samplesheet				Path to the samplesheet file (fastq paths are defined in its short_fastq_1/short_fastq_2 columns)
    
	Optional parameters:
		--threads				Number of threads (default=4)


    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

process fastp {
	cpus "${params.threads}"
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: '*trimmed.fastq.gz', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(reads1), path(reads2)
	output:
		tuple val(sample), path(reads1), path(reads2), path("R1_trimmed.fastq.gz"), path("R2_trimmed.fastq.gz"),  emit: trimmed_fastq
		path("fastp.log")
		path("*fastq.gz")
	when:
	!params.skip_fastp
	script:
	"""
	fastp -i ${reads1} -I ${reads2} -o R1_trimmed.fastq.gz -O R2_trimmed.fastq.gz
	cp .command.log fastp.log
	"""
}

process fastqc {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.zip'
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.html', saveAs: { filename -> "${sample}_$filename" }
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed)
        output:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), emit: reads_qc
		path("fastqc.log")
		path("*fastqc.zip"), emit: fastqc_zip
                path("*fastqc.html")
        when:
        !params.skip_fastqc
        script:
        """
        fastqc -o \$PWD ${reads1_trimmed} 
	fastqc -o \$PWD ${reads2_trimmed} 
        mv R1_trimmed_fastqc.zip ${sample}_R1_trimmed_fastqc.zip
	mv R2_trimmed_fastqc.zip ${sample}_R2_trimmed_fastqc.zip
	cp .command.log fastqc.log
        """
}

process summary_fastqc {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*html'
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*txt'
	input:
		path(fastqc_files)
	output:
		path("2_Illumina_multiqc_report.html"), emit: fastqc_summary
		path("2_Illumina_multiqc_general_stats.txt"), emit: fastqc_stats
	when:
	!params.skip_summary_fastqc
	script:
	"""
	multiqc --fn_as_s_name .
	cp .command.log summary_fastqc.log
	mv multiqc_report.html 2_Illumina_multiqc_report.html
	mv multiqc_data/multiqc_general_stats.txt 2_Illumina_multiqc_general_stats.txt
	"""
}

process shovill {
        cpus { task.attempt <= 2 ? "${params.shovill_threads}" : 1 }
        tag "${sample}"
        label "cpu"
        label "high_memory"
        errorStrategy 'retry'
        maxRetries 2
	publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: '*fa'
        input:
                tuple val(sample), path(reads1), path(reads2)
        output:
                tuple val(sample), path("*contigs.fa"), emit: assembly_out
		path("*contigs.fa"), emit: assembly_fasta
                path("shovill.log")
		path("*fa")
        script:
        """
        shovill --outdir \$PWD --R1 ${reads1} --R2 ${reads2} --gsize ${params.genome_size} --force --cpus ${task.cpus} ${params.shovill_args} --ram ${task.memory}
	mv contigs.fa ${sample}_contigs.fa
	mv spades.gfa ${sample}_contigs.gfa
        cp .command.log shovill.log
        """
}

process summary_shovill {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(shovill_fasta_files)
	output:
		path("3_Illumina_shovill_stats.tsv"), emit: shovill_summary
	script:
	"""
	echo -e "sample\tassembly_coverage\tnb_contigs\tassembly_size" > 3_Illumina_shovill_stats.tsv
	for file in `ls *contigs.fa`; do
		fileName=\$(basename \$file)
		sample=\${fileName%%_contigs.fa}
		grep "^>" \$file | sed s/len=// | sed s/cov=// > tmp
		total_length=`awk '{total_length+=\$2} END {print total_length}' tmp`
		total_cov=`awk '{total_cov+=\$2*\$3} END {print total_cov}' tmp`
		total_cov_decimal=\$(printf "%.0f" "\$total_cov")
		mean_cov=\$(echo "scale=0; \$total_cov_decimal / \$total_length" | bc -l)
		nb_contigs=`grep "^>" \$file | wc -l`
		echo -e \$sample\\\t\$mean_cov\\\t\$nb_contigs\\\t\$total_length  >> 3_Illumina_shovill_stats.tsv
	done
	"""
}

process quast {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
	publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: '*tsv'
	input:
                tuple val(sample), path(assembly)
        output:
		path("*report.tsv"), emit: quast_results
                path("quast.log")
        when:
        !params.skip_quast
        script:
        """
	quast.py ${assembly} --threads ${params.threads} -o \$PWD
	sed "s/_contigs\$//" report.tsv > ${sample}_report.tsv
        rm transposed_report.tsv report.tsv
	cp .command.log quast.log
        """
}

process summary_quast {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(quast_files)
	output:
		path("4_Illumina_quast_report.tsv"), emit: quast_summary	
	when:
	!params.skip_quast
	script:
	"""
	for file in `ls *report.tsv`; do cut -f2 \$file > \$file.tmp.txt; cut -f1 \$file > rownames.txt; done
	paste rownames.txt *tmp.txt > 4_Illumina_quast_report.tsv
	"""
}

process download_checkm_db {
	publishDir "${projectDir}/databases/",  mode: 'copy'
	output:
		path("checkm_data_2015_01_16"), emit: checkm_db_folder
	when:
	params.download_checkm_db
	script:
	"""
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
	mkdir checkm_data_2015_01_16
	tar -xvzf checkm_data_2015_01_16.tar.gz -C checkm_data_2015_01_16
	"""
}

process checkm {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly), path(checkm_db_folder)
        output:
                path("checkm.log")
		path("*checkm_lineage_wf_results.tsv"), emit: checkm_results
        when:
        !params.skip_checkm
        script:
        """
        export CHECKM_DATA_PATH=${params.checkm_db}
        checkm data setRoot ${params.checkm_db}
        checkm lineage_wf --reduced_tree `dirname ${assembly}` \$PWD --threads ${params.threads} --pplacer_threads ${params.threads} --tab_table -f checkm_lineage_wf_results.tsv -x fa
        mv checkm_lineage_wf_results.tsv ${sample}_checkm_lineage_wf_results.tsv
	cp .command.log checkm.log
        """
}

process summary_checkm {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(checkm_files)
	output:
		path("5_Illumina_checkm_lineage_wf_results.tsv"), emit: checkm_summary
	when:
	!params.skip_checkm
	script:
	"""
	echo -e  sampleID\\\tMarker_lineage\\\tNbGenomes\\\tNbMarkers\\\tNbMarkerSets\\\t0\\\t1\\\t2\\\t3\\\t4\\\t5+\\\tCompleteness\\\tContamination\\\tStrain_heterogeneity > header_checkm
	for file in `ls *checkm_lineage_wf_results.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_checkm_lineage_wf_results.tsv}; grep -v Bin \$file | sed s/_contigs//  >> 5_checkm_lineage_wf_results.tsv.tmp; done
	cat header_checkm 5_checkm_lineage_wf_results.tsv.tmp > 5_Illumina_checkm_lineage_wf_results.tsv
	"""
}

process sylph_download_db {
        cpus 1
        label "cpu"
        publishDir "$params.outdir/databases/sylph_database", mode: 'copy', pattern: "*.syldb"
        input:
                val(db)
        output:
                path("*.syldb"), emit: sylph_db
        when:
        !params.skip_download_sylph_db
        script:
        """
        echo "${db}"
        for attempt in 1 2 3 4 5; do
            wget -c -T 60 "${db}" && break
            [ \$attempt -lt 5 ] && sleep 30
        done
        """
}

process sylph {
        cpus "${params.sylph_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads1_trimmed), path(reads2_trimmed), path(db_files)
        output:
                tuple val(sample), path("*sylph_profile.tsv"), emit: sylph_profile
        when:
        !params.skip_sylph
        script:
        """
        sylph profile ${db_files.join(' ')} \
        -t ${params.sylph_threads} \
        -1 ${reads1_trimmed} -2 ${reads2_trimmed} \
        --output-file sylph_profile.tsv
        """
}

process sylph_tax_download_metadata {
        cpus 1
        label "cpu"
        publishDir "$params.outdir/databases/sylph_database", mode: 'copy', pattern: "*.gz"
        input:
                val(metadata_file)
        output:
                path("*.gz"), emit: sylph_tax_metadata
        when:
        !params.skip_download_sylph_db
        script:
        """
        for attempt in 1 2 3 4 5; do
            wget -c -T 60 "${metadata_file}" -O \$PWD/\$(basename $metadata_file) && break
            [ \$attempt -lt 5 ] && sleep 30
        done
        """
}

process sylph_tax {
        cpus "${params.sylph_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.sylphmpa", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(sylph_profile), path(metadata_files)
        output:
                tuple val(sample), path("merged_taxonomic_abundance.tsv"), path("merged_sequence_abundance.tsv"), emit: sylph_tax
        when:
        !params.skip_sylph
        script:
        """
        sylph-tax taxprof \
        ${sylph_profile} \
        -o \$PWD/ \
        -t ${metadata_files.join(' ')}

        # Merge taxonomy outputs
        sylph-tax merge \$PWD/*.sylphmpa \
        --column relative_abundance \
        -o \$PWD/merged_taxonomic_abundance.tsv

        sylph-tax merge \$PWD/*.sylphmpa \
        --column sequence_abundance \
        -o \$PWD/merged_sequence_abundance.tsv
        """
}

process sylph_summary_per_sample {
        input:
               tuple val(sample), path(taxonomic_abundances), path(sequence_abundances)
        output:
               tuple val(sample), path("${sample}_sylph_summary.tsv")
        script:
        """
        taxonomic_abundance_top_species=\$(grep "s__" ${taxonomic_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g")
        sequence_abundance_top_species=\$(grep "s__" ${sequence_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g")

        taxonomic_abundance_pasteurella_multocida=\$(grep "s__Pasteurella multocida" ${taxonomic_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g" | awk -F "\t" '{print \$2}')
        sequence_abundance_pasteurella_multocida=\$(grep "s__Pasteurella multocida" ${sequence_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g" | awk -F "\t" '{print \$2}')
        echo -e "${sample}\t\$taxonomic_abundance_top_species\t\$sequence_abundance_top_species\tPasteurella_multocida\t\$taxonomic_abundance_pasteurella_multocida\t\$sequence_abundance_pasteurella_multocida" > ${sample}_sylph_summary.tsv
        """
}

process summary_sylph {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(sylph_summary_files)
        output:
                path("6_Illumina_sylph_summary.tsv"), emit: sylph_summary
        when:
        !params.skip_sylph
        script:
        """
        echo -e "sample\ttop_species_by_taxonomic_abundance\ttaxonomic_abundance_for_top_species\ttop_species_by_sequence_abundance\tsequence_abundance_for_top_species\tPasteurella_multocida\ttaxonomic_abundance_for_pasteurella_multocida\tsequence_abundance_for_pasteurella_multocida" > 6_Illumina_sylph_summary.tsv
        for file in ${sylph_summary_files.join(' ')}; do
            cat \$file >> 6_Illumina_sylph_summary.tsv
        done
        """
}

process kaptive3 {
        errorStrategy 'ignore'
	cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*fna'
	publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly)
        output:
		tuple val(sample), path("*kaptive_results.tsv"), emit: kaptive_results
		path("*kaptive_results.tsv"),  emit: kaptive_tsv
		path("*fna")
                path("kaptive_v3.log")
        when:
        !params.skip_kaptive3
        script:
        """
	kaptive assembly ${params.kaptive_db_9lps} ${assembly} -f \$PWD -o kaptive_results.tsv
	mv kaptive_results.tsv ${sample}_kaptive_results.tsv
	sed s/_contigs// ${sample}_contigs_kaptive_results.fna > ${sample}_kaptive_results.fna
	rm ${sample}_contigs_kaptive_results.fna
	cp .command.log kaptive_v3.log
        """
}

process summary_kaptive {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(kaptive_files)
	output:
		path("7_Illumina_kaptive_results.tsv"), emit: kaptive_summary
	when:
	!params.skip_kaptive3
	script:
	"""
	echo -e sampleID\\\tBest match locus\\\tBest match type\\\tMatch confidence\\\tProblems\\\tIdentity\\\tCoverage\\\tLength discrepancy\\\tExpected genes in locus\\\tExpected genes in locus, details\\\tMissing expected genes\\\tOther genes in locus\\\tOther genes in locus, details\\\tExpected genes outside locus\\\tExpected genes outside locus, details\\\tOther genes outside locus\\\tOther genes outside locus, details\\\tTruncated genes, details\\\tExtra genes, details >  header_kaptive3
	for file in `ls *_kaptive_results.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_kaptive_results.tsv}; grep -v Assembly \$file | sed s/_contigs//  >> 7_kaptive_results.tsv.tmp; done
	cat header_kaptive3 7_kaptive_results.tsv.tmp > 7_Illumina_kaptive_results.tsv
	"""
}

process snippy {
	errorStrategy 'ignore'
        cpus "${params.snippy_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: 'snps*', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: '*tab'
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), path(kaptive_report)
	output:
                tuple val(sample), path("*snps.tab"), path("*snps.high_impact.tab"), path("snps.raw.vcf"), path("snps.filt.vcf"),  path("snps.bam"), path("snps.bam.bai"),  emit: snippy_results
		path("snippy.log")
		tuple path("*snps.tab"), path("*snps.high_impact.tab"), emit: snippy_impact_tab
        when:
        !params.skip_snippy && !params.skip_kaptive3
        script:
        """
	locus=\$(tail -1 "${kaptive_report}" | cut -f3)
	ref_gb=\$(grep \${locus:0:2} "${params.reference_LPS_directory}/reference_LPS.txt" | cut -f2)
	ref_gb="${params.reference_LPS_directory}/\$ref_gb"
	snippy --cpus ${params.snippy_threads} --force --outdir \$PWD --ref \${ref_gb} --R1 ${reads1_trimmed} --R2 ${reads2_trimmed} ${params.snippy_args}
        egrep "^CHROM|frameshift_variant|stop_gained" snps.tab > snps.high_impact.tab
	mv snps.high_impact.tab ${sample}_snps.high_impact.tab
	mv snps.tab ${sample}_snps.tab
	cp .command.log snippy.log
        """
}

process petg_blast {
        cpus "${params.petg_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/13_petG",  mode: 'copy', pattern: '*_petG_*'
        publishDir "$params.outdir/$sample/13_petG",  mode: 'copy', pattern: 'petg_blast.log', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(assembly), path(petg_reference)
        output:
                path("${sample}_petG_hits.fasta"), emit: petg_hits
                path("${sample}_petG_summary.tsv"), emit: petg_summary
                path("${sample}_petG_blast.tsv")
                path("${sample}_petG_blast.filtered.tsv")
                path("petg_blast.log")
        when:
        !params.skip_petg
        script:
        """
        makeblastdb -in ${assembly} -dbtype nucl -parse_seqids -out assembly_db
        blastn \\
                -query ${petg_reference} \\
                -db assembly_db \\
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
                -num_threads ${params.petg_threads} \\
                -out ${sample}_petG_blast.tsv

        awk -F '\\t' -v OFS='\\t' -v min_len="${params.petg_min_length}" -v min_ident="${params.petg_min_identity}" '
        {
                start = \$9
                end = \$10
                strand = "plus"
                suffix = "+"
                if (start > end) {
                        start = \$10
                        end = \$9
                        strand = "minus"
                        suffix = "-"
                }
                span = end - start + 1
                if (\$3 >= min_ident && span > min_len) {
                        print \$0, start, end, strand, suffix
                }
        }
        ' ${sample}_petG_blast.tsv > ${sample}_petG_blast.filtered.tsv

        touch ${sample}_petG_hits.fasta
        echo -e "SAMPLE\\tPETG_PRESENT" > ${sample}_petG_summary.tsv
        if [[ -s ${sample}_petG_blast.filtered.tsv ]]; then
                while IFS=\$'\\t' read qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore hit_start hit_end strand suffix; do
                        blastdbcmd -db assembly_db -entry "\$sseqid" -range "\${hit_start}-\${hit_end}" -strand "\$strand" | \\
                        awk -v header=">${sample}|\${sseqid}:\${hit_start}-\${hit_end}\${suffix}" 'BEGIN { print header } /^>/ { next } { print }' >> ${sample}_petG_hits.fasta
                done < ${sample}_petG_blast.filtered.tsv
                echo -e "${sample}\\tyes" >> ${sample}_petG_summary.tsv
        else
                echo -e "${sample}\\t" >> ${sample}_petG_summary.tsv
        fi

        cp .command.log petg_blast.log
        """
}

process empty_mlst_report_input {
        output:
                path("empty_mlst_report_input.txt"), emit: empty_mlst
        script:
        """
        touch empty_mlst_report_input.txt
        """
}

process empty_petg_report_input {
        output:
                path("empty_petg_report_input.txt"), emit: empty_petg
        script:
        """
        touch empty_petg_report_input.txt
        """
}

process report {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'	
	input:
		path(snippy_files)
		path(kaptive_summary)
		path(petg_summaries)
		path(mlst_files)
	output:
		tuple path("8_Illumina_snippy_snps.tsv"), path("8_Illumina_snippy_snps.high_impact.tsv"), path("10_Illumina_subtype_report.tsv"), emit: subtype_report	
	when:
	!params.skip_snippy && !params.skip_kaptive3
	script:
	"""
	echo -e sampleID\\\tCHROM\\\tPOS\\\tTYPE\\\tREF\\\tALT\\\tEVIDENCE\\\tFTYPE\\\tSTRAND\\\tNT_POS\\\tAA_POS\\\tEFFECT\\\tLOCUS_TAG\\\tGENE\\\tPRODUCT > header_snippy
	for file in `ls *_snps.high_impact.tab`; do fileName=\$(basename \$file); sample=\${fileName%%_snps.high_impact.tab}; grep -v EVIDENCE \$file | sed s/^/\${sample}\\\t/  >> 8_snippy_snps.high_impact.tsv.tmp; done
	cat header_snippy 8_snippy_snps.high_impact.tsv.tmp > 8_Illumina_snippy_snps.high_impact.tsv
	for file in `ls *_snps.tab`; do fileName=\$(basename \$file); sample=\${fileName%%_snps.tab}; grep -v EVIDENCE \$file | sed s/^/\${sample}\\\t/  >> 8_snippy_snps.tsv.tmp; done
	cat header_snippy 8_snippy_snps.tsv.tmp > 8_Illumina_snippy_snps.tsv

	echo -e "SAMPLE\\tPETG_PRESENT" > petg_lookup.tsv
	find . -maxdepth 1 -name '*_petG_summary.tsv' | sort | while read petg_file; do
		tail -n +2 "\$petg_file" >> petg_lookup.tsv
	done

	echo -e "SAMPLE\\tMLST" > mlst_lookup.tsv
	find . -maxdepth 1 -name '*_mlst.csv' | sort | while read mlst_file; do
		fileName=\$(basename "\$mlst_file")
		sample=\${fileName%%_mlst.csv}
		mlst_st=\$(awk -F',' 'NR == 1 {print \$3; exit}' "\$mlst_file")
		echo -e "\${sample}\\t\${mlst_st}" >> mlst_lookup.tsv
	done

	awk -F'\t' -v OFS='\t' '
	NR == 1 {
		for (i = 1; i <= NF; i++) {
			header[\$i] = i
		}
		sample_col = ("sampleID" in header) ? header["sampleID"] : 1
		locus_col = ("Best match locus" in header) ? header["Best match locus"] : 2
		confidence_col = ("Match confidence" in header) ? header["Match confidence"] : 4
		next
	}
	{
		if (\$confidence_col == "Typeable") {
			split(\$locus_col, a, "-")
			gsub("LPS", "L", a[1])
			print \$sample_col, a[1]
		} else {
			print \$sample_col, "untypeable"
		}
	}' "${kaptive_summary}" > kaptive_tmp

	subtype_db="${params.reference_LPS_directory}/LPS_subtype_database_v2.txt"
	phenotype_lookup="${params.reference_LPS_directory}/phenotype_lookup.tsv"
	if [[ -f "\$phenotype_lookup" ]]; then
		phenotype_lookup_input="\$phenotype_lookup"
	else
		touch phenotype_lookup_empty.tsv
		phenotype_lookup_input="phenotype_lookup_empty.tsv"
	fi

	awk -F '\\t' -v OFS='\\t' -v subtype_db="\$subtype_db" -v phenotype_lookup="\$phenotype_lookup_input" -v petg_lookup="petg_lookup.tsv" -v mlst_lookup="mlst_lookup.tsv" '
	function set_header_fields(    i) {
		for (i = 1; i <= NF; i++) {
			header[\$i] = i
		}
		db_type_col = ("TYPE" in header) ? header["TYPE"] : 1
		db_subtype_col = ("SUBTYPE" in header) ? header["SUBTYPE"] : 2
		db_isolate_col = ("ISOLATE" in header) ? header["ISOLATE"] : 3
		db_chrom_col = ("CHROM" in header) ? header["CHROM"] : 4
		db_pos_col = ("POS" in header) ? header["POS"] : 5
		db_vartype_col = ("VARTYPE" in header) ? header["VARTYPE"] : 6
		db_ref_col = ("REF" in header) ? header["REF"] : 7
		db_alt_col = ("ALT" in header) ? header["ALT"] : 8
		db_gene_col = ("GENE" in header) ? header["GENE"] : 9
		db_pheno_default_col = ("PHENOTYPE_DEFAULT" in header) ? header["PHENOTYPE_DEFAULT"] : 0
		db_pheno_multi_col = ("PHENOTYPE_MULTIPLE_SUBTYPES" in header) ? header["PHENOTYPE_MULTIPLE_SUBTYPES"] : 0
		db_note_col = ("NOTE" in header) ? header["NOTE"] : 10
	}
	function field_value(col) {
		return (col > 0 && col <= NF) ? \$col : ""
	}
	function clean_phenotype(value) {
		return (value == "NA") ? "" : value
	}
	function phenotype_from_rule(type, rule, parts) {
		if (rule == "" || rule == "NA") {
			return ""
		}
		split(rule, parts, "_")
		return (parts[2] == "") ? rule : type "_" parts[2]
	}
	function choose_phenotype(sample, type, default_phenotype, multi_phenotypes, rules, i, parts) {
		default_phenotype = clean_phenotype(default_phenotype)
		if (multi_phenotypes != "" && multi_phenotypes != "NA") {
			rule_count = split(multi_phenotypes, rules, ";")
			for (i = 1; i <= rule_count; i++) {
				split(rules[i], parts, "_")
				if (parts[1] != "" && ((sample SUBSEP parts[1]) in sample_subtype)) {
					return phenotype_from_rule(type, rules[i])
				}
			}
		}
		return default_phenotype
	}
	FILENAME == subtype_db {
		if (FNR == 1) {
			set_header_fields()
			next
		}
		key = field_value(db_chrom_col) OFS field_value(db_pos_col) OFS field_value(db_ref_col) OFS field_value(db_alt_col)
		db_count[key]++
		idx = key SUBSEP db_count[key]
		db_type[idx] = field_value(db_type_col)
		db_subtype[idx] = field_value(db_subtype_col)
		db_isolate[idx] = field_value(db_isolate_col)
		db_chrom[idx] = field_value(db_chrom_col)
		db_pos[idx] = field_value(db_pos_col)
		db_vartype[idx] = field_value(db_vartype_col)
		db_ref[idx] = field_value(db_ref_col)
		db_alt[idx] = field_value(db_alt_col)
		db_gene[idx] = field_value(db_gene_col)
		db_pheno_default[idx] = field_value(db_pheno_default_col)
		db_pheno_multi[idx] = field_value(db_pheno_multi_col)
		db_note[idx] = field_value(db_note_col)
		next
	}
	FILENAME == phenotype_lookup {
		if (FNR > 1 && \$1 != "") {
			phenotype_description[\$1] = \$2
		}
		next
	}
	FILENAME == petg_lookup {
		if (FNR > 1 && \$1 != "") {
			petg_present[\$1] = \$2
		}
		next
	}
	FILENAME == mlst_lookup {
		if (FNR > 1 && \$1 != "") {
			mlst_st[\$1] = \$2
		}
		next
	}
	FNR > 1 {
		sample = \$1
		key = \$2 OFS \$3 OFS \$5 OFS \$6
		for (i = 1; i <= db_count[key]; i++) {
			idx = key SUBSEP i
			row_count++
			row_sample[row_count] = sample
			row_idx[row_count] = idx
			sample_subtype[sample SUBSEP db_subtype[idx]] = 1
		}
	}
	END {
		print "SAMPLE", "MLST", "TYPE", "SUBTYPE", "VARTYPE", "ISOLATE_DATABASE", "CHROM", "POS", "REF", "ALT", "GENE", "PREDICTED_PHENOTYPE", "PREDICTED_PHENOTYPE_DESCRIPTION", "PETG_PRESENT", "NOTE"
		for (i = 1; i <= row_count; i++) {
			sample = row_sample[i]
			idx = row_idx[i]
			phenotype = choose_phenotype(sample, db_type[idx], db_pheno_default[idx], db_pheno_multi[idx])
			description = (phenotype in phenotype_description) ? phenotype_description[phenotype] : ""
			print sample, mlst_st[sample], db_type[idx], db_subtype[idx], db_vartype[idx], db_isolate[idx], db_chrom[idx], db_pos[idx], db_ref[idx], db_alt[idx], db_gene[idx], phenotype, description, petg_present[sample], db_note[idx]
		}
	}
	' "\$subtype_db" "\$phenotype_lookup_input" petg_lookup.tsv mlst_lookup.tsv 8_Illumina_snippy_snps.tsv > 10_Illumina_subtype_report.tsv.tmp
	awk -F'\t' 'NR > 1 {print \$1}' 10_Illumina_subtype_report.tsv.tmp | sort | uniq > list_samples_snippy_exclude
	awk -F'\t' 'FILENAME == ARGV[1] {exclude[\$1] = 1; next} !(\$1 in exclude)' list_samples_snippy_exclude kaptive_tmp > kaptive_to_keep
	awk -F'\t' -v OFS='\t' '
	FILENAME == "mlst_lookup.tsv" {
		if (FNR > 1 && \$1 != "") {
			mlst_st[\$1] = \$2
		}
		next
	}
	FILENAME == "petg_lookup.tsv" {
		if (FNR > 1 && \$1 != "") {
			petg_present[\$1] = \$2
		}
		next
	}
	{
		print \$1, mlst_st[\$1], \$2, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "", "", petg_present[\$1], ""
	}
	' mlst_lookup.tsv petg_lookup.tsv kaptive_to_keep > kaptive_to_keep.tsv
	cat 10_Illumina_subtype_report.tsv.tmp kaptive_to_keep.tsv > 10_Illumina_subtype_report.tsv
	"""
}

process html_report {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*.html'
	input:
		path(report_inputs)
		val(pipeline_version)
		val(skipped_steps)
		val(param_note)
	output:
		path("LPS_typing_report.html"), emit: html_report
	when:
	!params.skip_html_report && !params.skip_snippy && !params.skip_kaptive3
	script:
	"""
	generate_lps_report.py \\
		--report-dir . \\
		--lps-db-dir ${params.reference_LPS_directory} \\
		--out LPS_typing_report.html \\
		--pipeline-version "${pipeline_version}" \\
		--skipped "${skipped_steps}" \\
		--params "${param_note}"
	"""
}

process prepare_mlst_db {
	container 'docker://quay.io/biocontainers/mlst:2.33.1--hdfd78af_0'
	publishDir "${projectDir}/databases/", mode: 'copy'
	input:
		path(mlst_input_dir)
	output:
		path("mlst_db_prepared"), emit: mlst_db_prepared
	script:
	"""
	mkdir -p mlst_db_prepared/blast
	# Copy scheme files, renaming PubMLST .fas allele downloads to .tfa
	for scheme_dir in ${mlst_input_dir}/*/; do
	    scheme=\$(basename "\$scheme_dir")
	    mkdir -p mlst_db_prepared/\$scheme
	    for src in "\$scheme_dir"*; do
	        fname=\$(basename "\$src")
	        if [[ "\$fname" == *.fas ]]; then
	            cp "\$src" "mlst_db_prepared/\$scheme/\${fname%.fas}.tfa"
	        else
	            cp "\$src" "mlst_db_prepared/\$scheme/\$fname"
	        fi
	    done
	done
	# Build BLAST database with scheme-prefixed sequence headers (as mlst expects)
	for scheme_dir in mlst_db_prepared/*/; do
	    scheme=\$(basename "\$scheme_dir")
	    for tfa in "\$scheme_dir"*.tfa; do
	        [ -f "\$tfa" ] || continue
	        sed "s/^>/>\$scheme./" "\$tfa"
	    done
	done > mlst_db_prepared/blast/mlst.fa
	makeblastdb -hash_index -in mlst_db_prepared/blast/mlst.fa -dbtype nucl -title "PubMLST" -parse_seqids
	"""
}

process mlst {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: '*csv'
        input:
                tuple val(sample), path(assembly)
		path mlst_db
        output:
                path("mlst.log")
		path("*_mlst.csv"), emit: mlst_results
        when:
        !params.skip_mlst
        script:
	def custom_db_args = mlst_db ? "--datadir ${mlst_db} --blastdb ${mlst_db}/blast/mlst.fa" : ""
        """
	mlst --scheme ${params.mlst_scheme} ${assembly} --quiet --csv --threads ${params.threads} ${custom_db_args} > mlst.csv
        sed  s/_contigs.fa// mlst.csv > ${sample}_mlst.csv
	cp .command.log mlst.log
        """
}

process summary_mlst {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*csv'
	input:
		path(mlst_files)
	output:
		path("9_Illumina_mlst.csv"), emit: mlst_summary
	when:
	!params.skip_mlst
	script:
	"""
	for file in `ls *_mlst.csv`; do fileName=\$(basename \$file); sample=\${fileName%%_mlst.csv}; cat \$file >> 9_Illumina_mlst.csv; done
	"""
}

process download_bakta_db {
    publishDir "${projectDir}/databases/",  mode: 'copy'
    output:
        path("bakta_db"), emit: bakta_db_folder
    when:
    params.download_bakta_db
    script:
    """
    bakta_db download --output bakta_db --type full
    """
}

process bakta {
	cpus "${params.bakta_threads}"
	tag "${sample}"
	publishDir "$params.outdir/$sample/11_bakta",  mode: 'copy', pattern: "*.log"
	publishDir "$params.outdir/$sample/11_bakta",  mode: 'copy', pattern: '*bakta*'
	input:
		tuple val(sample), path(assembly), path(bakta_db_folder)
	output:
		path("*bakta*")
		path("bakta.log")
	when:
	!params.skip_bakta
	script:
	"""
	bakta --db ${params.bakta_db} --threads ${params.bakta_threads} --prefix ${sample}_bakta --output \$PWD/ ${params.bakta_args} ${assembly}
	cp .command.log bakta.log
	"""
}

process download_amrfinder_db {
	publishDir "${projectDir}/databases/amrfinderplus",  mode: 'copy'
	output:
		path("amrfinderplus_db"), emit: amrfinder_db_folder
	when:
	params.download_amrfinder_db
	script:
	"""
	amrfinder_update -d amrfinderplus_db
	"""
}

process amrfinder {
	tag "${sample}"
	publishDir "$params.outdir/$sample/12_amrfinder",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/12_amrfinder",  mode: 'copy', pattern: '*tsv'
	input:
		tuple val(sample), path(assembly), path(amrfinder_db_folder)
	output:
		path("*.tsv"), emit: amrfinder_results
		path("amrfinder.log")
	when:
	!params.skip_amrfinder
	script:
	"""
	amrfinder -n ${assembly} -d ${params.amrfinder_db} -o \$PWD/${sample}_amrfinder.tsv --name ${sample} --threads ${params.threads} --plus ${params.amrfinder_args}
	cp .command.log amrfinder.log
	"""
}

process summary_amrfinder {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(amrfinder_files)
	output:
		path("12_Illumina_amrfinder.tsv"), emit: amrfinder_summary
	when:
	!params.skip_amrfinder
	script:
	"""
	echo -e Name\\\tProtein id\\\tContig id\\\tStart\\\tStop\\\tStrand\\\tElement symbol\\\tElement name\\\tScope\\\tType\\\tSubtype\\\tClass\\\tSubclass\\\tMethod\\\tTarget length\\\tReference sequence length\\\t% Coverage of reference\\\t% Identity to reference\\\tAlignment length\\\tClosest reference accession\\\tClosest reference name\\\tHMM accession\\\tHMM description > header_amrfinder
	for file in ${amrfinder_files}; do 
		tail -n +2 "\$file" >> 12_amrfinder.tsv.tmp
	done
	cat header_amrfinder 12_amrfinder.tsv.tmp > 12_Illumina_amrfinder.tsv
	"""
}

workflow {
	Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
	.splitCsv(header:true, sep:',')
	.map { row -> tuple(row.sample_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
	.set { ch_samplesheet_illumina }
	fastp(ch_samplesheet_illumina)
	// When skip_fastp is true, pass raw reads through as if they were trimmed
	ch_trimmed = params.skip_fastp
		? ch_samplesheet_illumina.map { sample, r1, r2 -> tuple(sample, r1, r2, r1, r2) }
		: fastp.out.trimmed_fastq
	fastqc(ch_trimmed)
	summary_fastqc(fastqc.out.fastqc_zip.collect())
	ch_shovill_reads = ch_trimmed.map { sample, r1, r2, r1t, r2t -> tuple(sample, r1t, r2t) }
	shovill(ch_shovill_reads)
	summary_shovill(shovill.out.assembly_fasta.collect())
	quast(shovill.out.assembly_out)
	summary_quast(quast.out.quast_results.collect())
	if (!params.skip_sylph) {
		// Use the paired trimmed reads for taxonomy classification
		ch_sylph_reads = ch_trimmed.map { sample, r1, r2, r1t, r2t -> tuple(sample, r1t, r2t) }
		if (!params.skip_download_sylph_db) {
			// Download the Sylph reference databases
			ch_sylph_db = Channel.of("${params.sylph_db_gtdb_file}", "${params.sylph_db_fungal_file}")
			ch_downloaded_dbs = sylph_download_db(ch_sylph_db).sylph_db
			ch_db_list = ch_downloaded_dbs.collect()

			// Download the Sylph-taxa metadata
			ch_sylph_metadata = Channel.of("${params.sylph_tax_gtdb_metadata}", "${params.sylph_tax_fungal_metadata}")
			sylph_tax_download_metadata(ch_sylph_metadata).collect().set{ sylph_tax_metadata }

			ch_db_list.map { dbs -> tuple([dbs]) }.set { ch_db_tuple }
			sylph_tax_metadata.map { dbs -> tuple([dbs]) }.set { sylph_tax_metadata_tuple }
		} else {
			Channel.fromPath("${params.sylph_db}").collect().map { dbs -> tuple([dbs]) }.set { ch_db_tuple }
			Channel.fromPath("${params.sylph_metadata}").collect().map { dbs -> tuple([dbs]) }.set { sylph_tax_metadata_tuple }
		}

		// Run sylph
		ch_sylph_reads
			.combine(ch_db_tuple)
			.map { sample, r1t, r2t, dbs -> tuple(sample, r1t, r2t, dbs) }
			.set { ch_sylph_input }
		sylph(ch_sylph_input)

		// Run sylph-tax
		sylph.out.sylph_profile.combine(sylph_tax_metadata_tuple).set { ch_sylph_tax_input }
		sylph_tax(ch_sylph_tax_input)

		// Per-sample and aggregated summaries
		sylph_summary_per_sample(sylph_tax.out.sylph_tax)
			.map { sample, summary_file -> summary_file }
			.collect()
			.set { all_sylph_summaries }
		summary_sylph(all_sylph_summaries)
	}
	if (!params.skip_checkm) {
		if (params.download_checkm_db) {
			download_checkm_db()
			checkm(shovill.out.assembly_out.combine(download_checkm_db.out.checkm_db_folder))
		} else {
			checkm(shovill.out.assembly_out.combine(Channel.fromPath("${params.checkm_db}")))
		}
		summary_checkm(checkm.out.checkm_results.collect())
	}
	kaptive3(shovill.out.assembly_out)
	summary_kaptive(kaptive3.out.kaptive_tsv.collect())
	snippy(ch_trimmed.join(kaptive3.out.kaptive_results))
	snippy_tab_ch=snippy.out.snippy_impact_tab.collect()
	kaptive_summary_ch=summary_kaptive.out.kaptive_summary
	if (params.mlst_datadir) {
		prepare_mlst_db(Channel.value(file(params.mlst_datadir)))
		ch_mlst_db = prepare_mlst_db.out.mlst_db_prepared
	} else {
		ch_mlst_db = Channel.value([])
	}
	if (!params.skip_mlst) {
		mlst(shovill.out.assembly_out, ch_mlst_db)
		summary_mlst(mlst.out.mlst_results.collect())
		mlst_report_ch = mlst.out.mlst_results.collect()
	} else {
		empty_mlst_report_input()
		mlst_report_ch = empty_mlst_report_input.out.empty_mlst
	}
	if (!params.skip_petg && !params.skip_kaptive3) {
		ch_petg_reference = Channel.fromPath("${params.reference_LPS_directory}/petG_X73_NZ_CM001580.fasta", checkIfExists: true)
		petg_blast(shovill.out.assembly_out.combine(ch_petg_reference))
		petg_report_ch = petg_blast.out.petg_summary.collect()
	} else {
		empty_petg_report_input()
		petg_report_ch = empty_petg_report_input.out.empty_petg
	}
	report(snippy_tab_ch,kaptive_summary_ch,petg_report_ch,mlst_report_ch)
	if (!params.skip_bakta) {
		if (params.download_bakta_db) {
			download_bakta_db()
			bakta(shovill.out.assembly_out.combine(download_bakta_db.out.bakta_db_folder))
		} else {
			bakta(shovill.out.assembly_out.combine(Channel.fromPath("${params.bakta_db}")))
		}
	}
	if (!params.skip_amrfinder) {
		if (params.download_amrfinder_db) {
			download_amrfinder_db()
			amrfinder(shovill.out.assembly_out.combine(download_amrfinder_db.out.amrfinder_db_folder))
		} else {
			amrfinder(shovill.out.assembly_out.combine(Channel.fromPath("${params.amrfinder_db}")))
		}
		summary_amrfinder(amrfinder.out.amrfinder_results.collect())
	}
	// Build the single combined HTML report from the aggregated 10_report outputs.
	// The LPS database directory is read directly via params.reference_LPS_directory
	// (mounted), so only the aggregated TSV/HTML files need to be staged here.
	html_inputs_ch = report.out.subtype_report.flatten()
		.mix(summary_shovill.out.shovill_summary)
		.mix(summary_quast.out.quast_summary)
		.mix(summary_kaptive.out.kaptive_summary)
		.mix(summary_fastqc.out.fastqc_summary)
	if (!params.skip_checkm) {
		html_inputs_ch = html_inputs_ch.mix(summary_checkm.out.checkm_summary)
	}
	if (!params.skip_sylph) {
		html_inputs_ch = html_inputs_ch.mix(summary_sylph.out.sylph_summary)
	}
	if (!params.skip_mlst) {
		html_inputs_ch = html_inputs_ch.mix(summary_mlst.out.mlst_summary)
	}
	if (!params.skip_amrfinder) {
		html_inputs_ch = html_inputs_ch.mix(summary_amrfinder.out.amrfinder_summary)
	}
	skipped_list = []
	if (params.skip_sylph) skipped_list << 'sylph'
	if (params.skip_checkm) skipped_list << 'checkm'
	if (params.skip_quast) skipped_list << 'quast'
	if (params.skip_mlst) skipped_list << 'mlst'
	if (params.skip_petg) skipped_list << 'petG'
	if (params.skip_amrfinder) skipped_list << 'amrfinder'
	if (params.skip_bakta) skipped_list << 'bakta'
	param_note_str = "MLST scheme: ${params.mlst_scheme}; petG reported present when a hit spans > ${params.petg_min_length} bp at >= ${params.petg_min_identity}% identity"
	html_report(html_inputs_ch.collect(), workflow.manifest.version, skipped_list.join(', '), param_note_str)
}
