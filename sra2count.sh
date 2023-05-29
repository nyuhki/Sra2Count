#!/usr/sh

##########################################################################################
##sra2count.sh
#Usage:
#sh RNA.seq -i sra_file.txt -s mm -d paired -c 5
##########################################################################################
#Defining option
function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -i VALUE    input sra_file.txt
    -s VALUE    Species (Homo Sapience:hs, Mus Musculus:mm)
    -d VALUE    Direction of sequence (pair end:paired, single end:single)
    -c VALUE	CPU core
EOM
    exit 2
}

#Definition
while getopts ":i:s:d:c:h" optKey; do
    case "$optKey" in
		i)
          i=$(echo "${OPTARG}")
          ;;
        s)
          s=$(echo "${OPTARG}")
          ;;
        d)
          d=$(echo "${OPTARG}")
          ;;
		c)
          c=$(echo "${OPTARG}")
          ;;
        '-h'|'--help'|* )
          usage
          ;;
    esac
done

#Environment
echo "augumentation: input=$i, species=$s, direction=$d, core=$c"
mkdir ./result
LIST=$(cat ${i})
#SRA->fastq
for SampleID in `echo ${LIST}`
    do
		DIR=./result/${SampleID}
		if [ -e $DIR/${SampleID}.sra ]; then
			echo "${SampleID}.sra exists."
		else
			prefetch ${SampleID} -p
		fi
		mv ./${SampleID} ./result
		if [ -e $DIR/${SampleID}.fastq ]; then
			echo "${SampleID}.sra have been converted."
		elif [ -e $DIR/${SampleID}_1.fastq ]; then
			echo "${SampleID}.sra have been converted."
		else
			if [ "${d}" = "paired" ]; then
				fasterq-dump $DIR/${SampleID}.sra --outdir $DIR --split-files --threads ${c} --progress
				echo "${SampleID}.sra have been converted."
			else
				fasterq-dump $DIR/${SampleID}.sra --outdir $DIR --threads ${c} --progress
				echo "${SampleID}.sra have been converted."
			fi
		fi
done

#FastQC(1st)>Trimming>FastQC(2nd)>Alignment>Count
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		mkdir $DIR/${SampleID}_fastqc \
		$DIR/${SampleID}_fastqc/1st \
		$DIR/${SampleID}_fastqc/2nd
		if [ "${d}" = "paired" ]; then
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}_1.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}_2.fastq
			./TrimGalore-0.6.10/trim_galore \
				-j ${c} \
				-q 30 \
				--paired \
				-o $DIR \
				$DIR/${SampleID}_1.fastq \
				$DIR/${SampleID}_2.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}_1_val_1.fq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}_2_val_2.fq
			if [ "${s}" = "hs" ]; then
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/hg38 \
					-1 $DIR/${SampleID}_1_val_1.fq \
					-2 $DIR/${SampleID}_2_val_2.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/hg38.gtf \
					-o $DIR/${SampleID}.gtf
			else
					./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/mm10 \
					-1 $DIR/${SampleID}_1_val_1.fq \
					-2 $DIR/${SampleID}_2_val_2.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/mm10.gtf \
					-o $DIR/${SampleID}.gtf
			fi
		else
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}.fastq
			./TrimGalore-0.6.10/trim_galore \
				-j ${c} \
				-q 30 \
				-o $DIR \
				$DIR/${SampleID}.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}.fq
			if [ "${s}" = "hs" ]; then
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/hg38 \
					-U $DIR/${SampleID}_trimmed.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/hg38.gtf \
					-o $DIR/${SampleID}.gtf
			else
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/mm10 \
					-U $DIR/${SampleID}_trimmed.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/mm10.gtf \
					-o $DIR/${SampleID}.gtf
			fi
		fi
done
#PrepforBallgown,DEseq
mkdir ./summary \
	  ./result/ballgown
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		mkdir ./result/ballgown/${SampleID}
		echo $DIR/${SampleID}.gtf \
			>> merge_list.txt
		echo ${SampleID} ./result/ballgown/${SampleID}/${SampleID}.bg.gtf \
			>> ballgown_PATH.txt
done
#Merge
if [ "${s}" = "hs" ]; then
	./stringtie/stringtie --merge -p ${c} \
		-G ./index/hg38.gtf \
		-o ./summary/merge.gtf \
		merge_list.txt
else
	./stringtie/stringtie --merge -p ${c} \
		-G ./index/mm10.gtf \
		-o ./summary/merge.gtf \
		merge_list.txt
fi
#Ballgown file genration
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		./stringtie/stringtie -p ${c} -e -B \
			$DIR/${SampleID}.bam \
			-G ./summary/merge.gtf \
			-o ./result/ballgown/${SampleID}/${SampleID}.bg.gtf \
			-A ./summary/${SampleID}.tab
done
#2to3 -w ./stringtie/prepDE.py
python ./stringtie/prepDE.py3 -i ballgown_PATH.txt
mv gene_count_matrix.csv ./result/gene_raw_count.csv
mv ./result/transcript_count_matrix.csv ./result/transcript_raw_count.csv
#Estimation of annotation accuracy
if [ "${s}" = "hs" ]; then
	./gffcompare/gffcompare \
		-r ./index/hg38.gtf \
		-G \
		-o ./summary/gffcompare ./summary/merge.gtf
else
	./gffcompare/gffcompare \
		-r ./index/mm10.gtf \
		-G \
		-o ./summary/gffcompare ./summary/merge.gtf
fi

rm merge_list.txt ballgown_PATH.txt
mv ./summary ./result
mv sra_file.txt ./result
exit=0
