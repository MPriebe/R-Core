
# Loop which takes accession numbers entered and carries out DGEA
while read p; do

 # create a directory with the Results
 mkdir GEO_${p}
 
 # make a name file with the accession id
 echo $p > GEO_${p}.txt

 # run Rscript to carry out Differential Gene Expression Analysis
 	  Rscript DGEA.R GEO_${p}.txt

done < accession.id.txt



