#!/usr/bin/env bash
#This requires fqgrep (https://github.com/indraniel/fqgrep)
usage () {
  echo "usage: $0 [-i <IR seq>] [-a <assembly>] <pfx> "
  echo "Required parameters:"
  echo "-i     The sequence of the transposon end sequence remaining (for junction authentication)"
  echo "-a     The name of the assembly you're using (PAO1)"
  echo "-m     The number of mismatches/indels you want to tolerate during search"

  echo ""
  echo "The required parameters must precede the file prefix for your sequence file:"
  echo "  (e.g. if your sequence file is named condition1.fastq,"
  echo "   the prefix is \"condition1\")"
  echo ""
  echo "Example:"
  echo "$0 -i TATAAGAGTCAG -a PAO1 -m 1 condition1"
}

# Read in the important options
while getopts ":i:a:m:" option; do
  case "$option" in
  	i)  IR="$OPTARG" ;;
  	a)  ASSEMBLY="$OPTARG" ;;
  	m)  MISMATCHES="$OPTARG" ;;
    h)  # it's always useful to provide some help 
        usage
        exit 0 
        ;;
    :)  echo "Error: -$option requires an argument" 
        usage
        exit 1
        ;;
    ?)  echo "Error: unknown option -$option" 
        usage
        exit 1
        ;;
  esac
done    
shift $(( OPTIND - 1 ))

# Do some error checking to make sure parameters are defined
if [ -z "$IR" ]; then
  echo "Error: you must specify the Tn end sequence using -i"
  usage
  exit 1
fi

if [ -z "$ASSEMBLY" ]; then
  echo "Error: you must specify an assembly using -a"
  usage
  exit 1
fi

if [ -z "$MISMATCHES" ]; then
  echo "Error: you must specify a number of mismatches using -m"
  usage
  exit 1
fi

# Give the usage if there aren't enough parameters
if [ $# -lt 1 ] ; then
  echo "you must provide a file prefix for analysis"
  usage
  exit 1
fi

PREFIX=$1
R1=${PREFIX}
BOWTIEREF=$REFGENOME$ASSEMBLY/$ASSEMBLY

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt
echo "Total sequences: " >> $PREFIX-TnSeq.txt
egrep -c '^@HWI|^@M|^@NS|^@SRR|^@NB5|^@M0' $R1.fastq >> $PREFIX-TnSeq.txt

# IRs
echo "$PREFIX: Searching for reads with an IR in right location..."
fqgrep -m $MISMATCHES -r -p $IR $R1.fastq | awk -F "\t" '(($8 >= 33 && $8 <= 44) || $1=="read name")' | trimmer --5-prime > $PREFIX-IR-clip.fastq
~/.local/bin/cutadapt --nextseq-trim=20 -a G{16} -m 20 $PREFIX-IR-clip.fastq > $PREFIX-IR-clip.trim.fastq 2> $PREFIX-cutadapt-report.txt
mv $PREFIX-IR-clip.trim.fastq $PREFIX-IR-clip.fastq
IRSFOUND=$(egrep -c '^@HWI|^@M|^@NS|^@SRR|^@NB5|^@M0' $PREFIX-IR-clip.fastq)
echo "Processed sequences:" >> $PREFIX-TnSeq.txt
echo $IRSFOUND >> $PREFIX-TnSeq.txt

# Map and convert - feel free to change bowtie2 parameters yourself
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
bowtie2 --end-to-end -p 16 -x $BOWTIEREF -U $PREFIX-IR-clip.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt
grep '^@' $PREFIX.sam > $PREFIX-mapped.sam
cat $PREFIX.sam | grep -v '^@' | awk -F "\t" '((and($2, 0x4) != 0x4) && ($5 > 39))' >> $PREFIX-mapped.sam
echo "Number of reads mapping at high enough MAPQ:" >> $PREFIX-TnSeq.txt
grep -v '^@' $PREFIX-mapped.sam | wc -l >> $PREFIX-TnSeq.txt
echo "$PREFIX: Tallying mapping results..."
grep -v '^@' $PREFIX-mapped.sam | awk -F "\t" '{if (and($2, 0x10) != 0x10) print $4; else print $4+length($10)+2}' | grep '[0-9]' | sort | uniq -c | sort -n -r > $PREFIX-sites.txt
echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
wc -l $PREFIX-sites.txt >> $PREFIX-TnSeq.txt
echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
head $PREFIX-sites.txt >> $PREFIX-TnSeq.txt

# Sort output, cleanup
echo "$PREFIX: Cleaning up..."
mkdir $PREFIX-results
mv $PREFIX-TnSeq.txt $PREFIX-results/
rm $PREFIX-IR-clip.fastq
mv $PREFIX.sam $PREFIX-results/
mv $PREFIX-mapped.sam $PREFIX-results/
mv $PREFIX-sites.txt $PREFIX-results/
mv $PREFIX-cutadapt-report.txt $PREFIX-results/
