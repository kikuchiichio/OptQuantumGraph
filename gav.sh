for t  in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5
do 
#	echo $t
	./a.out $t|grep GENE
done