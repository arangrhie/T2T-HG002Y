module add conda-modules-py37
module add bedtools-2.26.0
module load emboss-6.5.7

P4_arm1="P4_arm1::chrY_hg002:19356701-19546840.fasta"
P4_arm2="P4_arm2::chrY_hg002:19586495-19776671.fasta"
P5_arm1="P5_arm1::chrY_hg002:18362331-18857802.fasta"
P5_arm2="P5_arm2::chrY_hg002:18861259-19356680.fasta"
P6_arm1="P6_arm1::chrY_hg002:17066094-17176113.fasta"
P6_arm2="P6_arm2::chrY_hg002:17222322-17332346.fasta"
P7_arm1="P7_arm1::chrY_hg002:16781335-16790043.fasta"
P7_arm2="P7_arm2::chrY_hg002:16802684-16811393.fasta"
P8_arm1="P8_arm1::chrY_hg002:14891106-14926304.fasta"
P8_arm2="P8_arm2::chrY_hg002:14929718-14964927.fasta"


#align arm1 with a reverse of arm2
stretcher -sreverse2 -asequence ${P4_arm1} -bsequence ${P4_arm2} -outfile P4.hg002.selfalign &
stretcher -sreverse2 -asequence ${P5_arm1} -bsequence ${P5_arm2} -outfile P5.hg002.selfalign &
stretcher -sreverse2 -asequence ${P6_arm1} -bsequence ${P6_arm2} -outfile P6.hg002.selfalign &
stretcher -sreverse2 -asequence ${P7_arm1} -bsequence ${P7_arm2} -outfile P7.hg002.selfalign &
stretcher -sreverse2 -asequence ${P8_arm1} -bsequence ${P8_arm2} -outfile P8.hg002.selfalign &

wait

echo "Done."

