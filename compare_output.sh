
echo "Diffing $1 $2"
sort $1 > temp1.out
sort $2 > temp2.out
diff temp1.out temp2.out