module add conda-modules-py37
module add bedtools-2.26.0
module load emboss-6.5.7

P4_arm1="P4_arm1::chrY_hg002:19356690-19546842.fasta"
P4_arm2="P4_arm2::chrY_hg002:19586485-19776673.fasta"
P5_arm1="P5_arm1::chrY_hg002:18362320-18857804.fasta"
P5_arm2="P5_arm2::chrY_hg002:18861249-19356682.fasta"
P6_arm1="P6_arm1::chrY_hg002:17066083-17176115.fasta"
P6_arm2="P6_arm2::chrY_hg002:17222312-17332348.fasta"
P7_arm1="P7_arm1::chrY_hg002:16781324-16790045.fasta"
P7_arm2="P7_arm2::chrY_hg002:16802674-16811395.fasta"
P8_arm1="P8_arm1::chrY_hg002:14891095-14926306.fasta"
P8_arm2="P8_arm2::chrY_hg002:14929707-14964929.fasta"


#align arm1 with a reverse of arm2
stretcher -sreverse2 -asequence ${P4_arm1} -bsequence ${P4_arm2} -outfile P4.hg002.selfalign &
stretcher -sreverse2 -asequence ${P5_arm1} -bsequence ${P5_arm2} -outfile P5.hg002.selfalign &
stretcher -sreverse2 -asequence ${P6_arm1} -bsequence ${P6_arm2} -outfile P6.hg002.selfalign &
stretcher -sreverse2 -asequence ${P7_arm1} -bsequence ${P7_arm2} -outfile P7.hg002.selfalign &
stretcher -sreverse2 -asequence ${P8_arm1} -bsequence ${P8_arm2} -outfile P8.hg002.selfalign &

wait

echo "Done."