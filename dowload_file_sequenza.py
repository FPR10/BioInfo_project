import wget
from Bio import SeqIO

#Lettura del file + dowload su "filename"
url = 'https://siloe.dimes.unical.it/~fab/MedicinaTD/group_1_seq.fasta'
filename = wget.download(url)

'''
#Apertura e print del file
with open(filename, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        print(record.seq)
'''
