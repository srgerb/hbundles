#In this round I used the segments python script to generate helical bundles that match up with the base of Peilong's designs, and extend and flower out on the side opposite the helical repeats.

#I then used Peilong's design protocol to design and close the loops.

#In the future I would like to remodel the ends of the loops, to make the connection better, since some of the curvature of the loops look bad.

#in order to retreive the scores, make a headernames file that lists all the 8-digit names, then run this line (make sure the second argument goes to destination file)
while read line; do cp /projects/boinc-results/$line/*.sc.bz2 runs/scores; done < headernames

rename the files:
sh rename

#To make plot file:
ls scores_fixed/fold* > fold; ls scores_fixed/relax_score_* > relax; paste fold relax | sed "s/^/python plot_ff\.py /g" > plot.sh; rm fold; rm relax

#Unzip all the things!
bzip2 -d *.sc

#if there are failures and weird W00 things:
sed -i "/FAILURE/d" *.sc.out
sed -i "/W_000/d" *.sc.out

#to plot **make sure you have /plots/$folder_with_same_name_as_score_file_folder** for the plots to go into:
sh plot.sh

To extract pdb files:
get silent files from results:
while read line; do cp /projects/boinc-results/$line/*.out.bz2 runs/scores; done < headernames

Get top 10 scores for a file:
sort -n -k2 GCBETAB2_fold_SAVE_ALL_OUT_527272_0.sc | head -n 10 | awk '{print $2 "\t" $31 "\t" $35}' > score_rmsd.dat

source software:
source /software/rosetta/setup.sh

extract pdbs:
extract_pdbs -in:file:silent_struct_type binary -in:file:fullatom -fail_on_bad_hbond false -silent_read_through_errors -in:file:tags specific_decoy_tag -in:file:silent silent_file_example.out
