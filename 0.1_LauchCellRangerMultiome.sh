## data for the H1X sample
# cellranger-arc count --id=H1X_OLD \
#                        --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                        --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1X_old.csv \
#                        --localcores=56 \
#                        --localmem=200


# # data for the H1_3 sample
# nohup cellranger-arc count --id=H1_3_multiome \
#                        --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                        --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1_3.csv \
#                        --localcores=56 \
#                        --localmem=200 &

# data for the H1_WT sample
nohup cellranger-arc count --id=H1_WT_multiome \
                       --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1_WT.csv \
                       --localcores=56 \
                       --localmem=200 > nohup_H1_WT_multiome & 

# # data for the H1_WT sample
# nohup cellranger-arc count --id=H1_WT_RUN1_multiome \
#                        --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                        --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1_WT_RUN1.csv \
#                        --localcores=56 \
#                        --localmem=200 > nohup_H1_WT_multiome_RUN1 & 



# # data for the H1X_H1_3 sample
# nohup cellranger-arc count --id=H1X_H1_3_multiome \
#                        --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                        --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1X_H1_3.csv \
#                        --localcores=56 \
#                        --localmem=200 > nohup_H1X_H1_3_multiome &


# # data for the H1X sample
# nohup cellranger-arc count --id=H1X_multiome \
#                        --reference=/home/sevastopol/data/gserranos/CellRangerReferences/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                        --libraries=/home/sevastopol/data/gserranos/H1_kaust/Data/Library_Files/library_file_H1X.csv \
#                        --localcores=56 \
#                        --localmem=200 > nohup_H1X_multiome &

