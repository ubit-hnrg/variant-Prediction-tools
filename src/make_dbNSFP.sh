
#!/bin/sh -e

# These variables need to be updated with every version
#base="dbNSFP4.0b1a"                                	# Database version
dbZip=$1   						#full path to zip file
BASENAM=$(basename $dbZip)
base=${BASENAM%.zip}
#path=$(readlink -f $(dirname $1))
#base=$(basename $dbZip)
db=$base.txt          	                       		 # Output file
db_hg19=$base.hg19.txt

echo $dbZip
echo $base
echo $db

# Check dbNSFP
if [ ! -e "$dbZip" ]
then
	echo "ERROR: Expected dbNSFP zip file '$dbZip' not found"
	exit 1
fi

#---
# Create DB
#---
##echo Create file $db
unzip -p $dbZip $base\_variant.chr1.gz |zcat|head -n1 > $db
unzip -p $dbZip $base\_variant.chr1.gz \
		$base\_variant.chr2.gz \
		$base\_variant.chr3.gz \
		$base\_variant.chr4.gz \
		$base\_variant.chr5.gz \
		$base\_variant.chr6.gz \
		$base\_variant.chr7.gz \
		$base\_variant.chr8.gz \
		$base\_variant.chr9.gz \
		$base\_variant.chr10.gz \
		$base\_variant.chr11.gz \
		$base\_variant.chr12.gz \
		$base\_variant.chr13.gz \
		$base\_variant.chr14.gz \
		$base\_variant.chr15.gz \
		$base\_variant.chr16.gz \
		$base\_variant.chr17.gz \
		$base\_variant.chr18.gz \
		$base\_variant.chr19.gz \
		$base\_variant.chr20.gz \
		$base\_variant.chr21.gz \
		$base\_variant.chr22.gz \
		$base\_variant.chrX.gz \
		$base\_variant.chrY.gz \
	| zgrep -v "^#" \
	>> $db


echo BGZIP $db
bgzip $db -@ 4

echo TABIX $db.gz
tabix -s 1 -b 2 -e 2 $db.gz


#cat $db | ../src/dbNSFP_sort.pl 8 9 >$db_hg19

#echo BGZIP $db_hg19
#bgzip $db_hg19

#echo TABIX $db_hg19.gz
#tabix -s 1 -b 2 -e 2 $db_hg19.gz

