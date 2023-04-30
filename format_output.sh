cd output || exit
cp output.txt predictions.csv
zip -r predictions.zip predictions.csv
