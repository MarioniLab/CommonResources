from Bio import SeqIO
input_handle = open("RepeatMaskerLib.embl", "rU")
output_handle = open("RepeatMaskerLib.fasta", "w")
 
sequences = SeqIO.parse(input_handle, "embl")
count = SeqIO.write(sequences, output_handle, "fasta")
   
output_handle.close()
input_handle.close()

