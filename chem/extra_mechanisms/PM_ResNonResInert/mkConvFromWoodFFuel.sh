
for i in PM*
do
  n=$(echo $i| sed "s/WoodFFuelInert/ResNonRes/")
  echo $n
  sed "s/ffuel/nonRes/g" $i | sed "s/wood/Res/g" > $n
done
