echo "what is the name of your InSilico Digest file?"
read foo < /dev/tty
mkdir $foo.forR/
echo "how many scaffolds or contigs or chromosomes are you digesting?"
read moo < /dev/tty
cat $foo |grep -A1 '# Sequence: scaffold' > mylist
cat mylist |grep 'Sequence' |awk -F ":" '{print $2,$4}' |sed -e 's/from  //g' > lengthContigs.tab
cat lengthContigs.tab |awk '{print $1}' > contigsNames.txt
cat mylist |grep 'HitCount' |awk -F ":" '{print $2}' > NumHits.txt
for ((i=1;i<=$moo;i++)); do
	myscaffold=$(sed -n ''$i'p' contigsNames.txt)
	NumHits=$(sed -n ''$i'p' NumHits.txt)
	mynumHits=$(($NumHits + 13))
	mynumHits2=$(($NumHits + 1))
	echo $myscaffold $mynumHits $mynumHits2
	cat $foo |grep -A $mynumHits -w ''$myscaffold'' |tail -n $mynumHits2 |awk '{print "'$myscaffold'",$0}' |awk '!x{x=sub("'$myscaffold'","CHR")}7' > $foo.forR/$myscaffold.tab
done
rm contigsNames.txt lengthContigs.tab mylist NumHits.txt
