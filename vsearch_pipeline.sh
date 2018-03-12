outlocation=$(mktemp -d /media/GalaxyData/files/XXXXXX)
#outlocation = $(mktemp -d /media/GalaxyData/files/XXXXXX)
#$1=-i
#$2=-t
#$3=-cluster_id
#$4=-cluster_size
vsearch_pipeline.py -i $1 -o $outlocation -t $2 -cluster_id $3 -cluster_size $4
mv $outlocation"/adminlog.log" $5
mv $outlocation"/workfiles/otutab.txt" $6
mv $outlocation"/workfiles/filtered_otu.fa" $7
#rm -rf $outlocation 