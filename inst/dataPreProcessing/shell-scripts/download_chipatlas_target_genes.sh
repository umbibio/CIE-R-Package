for ag in `cat antigens.txt`; do 
    wget -c http://dbarchive.biosciencedbc.jp/kyushu-u/hg19/target/${ag}.5.tsv
done

