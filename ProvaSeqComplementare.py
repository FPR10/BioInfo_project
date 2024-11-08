import Bio.SeqIO
from Bio.Seq import Seq

# Sequenza di esempio
sequenza = Seq("ATCG")
sequenza_complementare = sequenza.complement()
sequenza_complementare_inversa = sequenza.reverse_complement()
print("Sequenza originale:", sequenza)
print("Sequenza complementare:", sequenza_complementare)
print ("Sequenza complementare letta in senso contrario:", sequenza_complementare_inversa)

