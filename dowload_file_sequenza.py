import wget
from Bio import SeqIO

# Lettura del file + download su "filename"
url = 'https://siloe.dimes.unical.it/~fab/MedicinaTD/group_1_seq.fasta'
filename = wget.download(url)

# Apertura e calcolo della lunghezza della sequenza
with open(filename, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        sequence_length = len(record.seq)
        print("Lunghezza della sequenza:", f"{sequence_length:,}")
