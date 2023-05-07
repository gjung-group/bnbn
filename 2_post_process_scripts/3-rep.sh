 tail -n +9 BNBN.cart > 1.txt
 sed -n 1,16p BLBL.mol > title.txt 
 tail -n +17 BLBL.mol > 3.txt
 awk '{print $1"\t" $2"\t" }' 3.txt > 4.txt
 awk '{print $8"\t" $9"\t" $10"\t"}' 3.txt > 5.txt
 paste 4.txt 1.txt 5.txt >> para.txt
 cat title.txt para.txt >> twist.mol
 sed -i '12s/12.010700226/10.80000191/g' twist.mol
 sed -i '13s/12.010700226/14.00000000/g' twist.mol
 #cp h-bn-BL.mol twist.mol
