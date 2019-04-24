#!/usr/bin/env bash
###input: sites file with number of reads mapping to insertion in first column and genome location in second column
###corrects sites file to be tab deliminated and then runs slippage analysis.
###slippage analysis collapses reads mapping within 1 bp to local maximum
usage () {
  echo "usage: $0 <pfx> "
  echo "Required parameters:"
  echo "The required parameters must precede the file prefix for your sequence file:"
  echo "  (e.g. if your sequence file is named condition1-sites.txt,"
  echo "   the prefix is \"condition1\")"
  echo ""
  echo "Example:"
  echo "$0 condition1"
}

PREFIX=$1
R1=${PREFIX}-sites

echo "Performing slippage analysis on $PREFIX..."
echo "slippage stats for $PREFIX" > $PREFIX-slippage-stats.txt
echo "Original_sites: " >> $PREFIX-slippage-stats.txt
wc -l $R1.txt >> $PREFIX-slippage-stats.txt

sed -i 's/       //g' $R1.txt
sed -i 's/      //g' $R1.txt
sed -i 's/     //g' $R1.txt
sed -i 's/    //g' $R1.txt
sed -i 's/   //g' $R1.txt
sed -i 's/  //g' $R1.txt
sed -i 's/ /\t/g' $R1.txt

python slippage_ed.py $R1.txt > $R1-slippage.txt

echo "Corrected_sites: " >> $PREFIX-slippage-stats.txt
wc -l $R1-slippage.txt >> $PREFIX-slippage-stats.txt
