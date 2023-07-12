# step1: download the Hasel et al. (2021) Nature Neuroscience 24:1475–1487 data
cd path/to/project
wget -O GSE148822_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167494&format=file'

# step2: unzip the data
gzip GSE148822_RAW.tar
tar -xvzf GSE148822_RAW.tar.gz

## Remove unneeded data
rm *SOX*
rm GSE167492_RAW.tar.gz

## Check the data
gzip -cd GSM5106098_D1-1_barcodes.tsv.gz

## Put individual data into individual directory
for i in $(seq 1 17)
do
 mkdir Donor$i
 mv *D$i* Donor$i
done


## Edit the file names
for dir in .../Step9.IntegrationOfPublicData/Sadick/Don*
do
cd $dir
for i in *.gz
do 
 newName=$(echo "$i" | tr -cd "barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz")
 mv $i $newName
done
done


### end of the code



