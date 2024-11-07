#####################################################
''' 
Identificativo del gruppo
group_id = 1

Partecipanti
- Mariaida Perri - mat. 239943
- Andreea Roxana Manolache - mat. 245906
- Francesco Pio Ruffo - mat. 240044
'''
#####################################################

# Importazione moduli
import numpy as np
import pandas as pd
import os
import regex
import Bio
import Bio.SeqIO
from bisect import bisect_left
import pickle


# CARICAMENTO DELLA SEQUENZA FASTA
# scaricare il file
# siloe.dimes.unical.it/~fab/group_<group_id>_seq.fasta
os.system("wget -q https://siloe.dimes.unical.it/~fab/MedicinaTD/group_1_seq.fasta")  

# RICERCA DELLE ORF
# orf_iter restituisce un iteratore sulle orf
# le modalità d'uso sono equivalenti al find_iter di regex

class Match:
	def __init__(self, string, start, stop):
		self.string = string
		self.start_pos = start
		self.stop_pos = stop
		self._group0 = string[start:stop]
	def group(self):
		return self._group0
	def span(self):
		return (self.start_pos, self.stop_pos)
	def __repr__(self):
		return f"<Match object; span={self.span()}, match='{self.group()}'>"

def orf_iter(sequence):
	# Cerco tutti i codoni di stop e salvo il punto di partenza
	stops = []
	for stop in regex.finditer(r"TAA|TAG|TGA", sequence):
		stops.append(stop.span()[0])
	# Itero sugli start
	for start in regex.finditer(r"ATG", sequence):
		# cerco il primo stop con valore di distanza multiplo di 3
		pos_start = start.span()[0]
		pos_stop = bisect_left(stops, pos_start + 3)
		while (pos_stop < len(stops)) and ((stops[pos_stop] - pos_start) % 3 != 0):
			pos_stop = pos_stop + 1
		if pos_stop < len(stops):
			yield Match(sequence, pos_start, stops[pos_stop])
   
   
################################################################################
#                         ALGORITMO DI GENE PREDICTION                         #
# Scopo del progetto è, data la sequenza scaricata, salvare in un dataframe    #
# la posizione delle ORF presenti nella sequenza e l'indicazione se si tratta  #
# di regione codificante o non codificante e se è stata verificata su BLAST    #
################################################################################

#
# CARICARE LA SEQUENZA DAL FILE FASTA
#
seq = next(Bio.SeqIO.parse(f"group_1_seq.fasta", "fasta")).seq

#
# PREPARAZIONE DEL DATAFRAME CON I RISULTATI
#
results = pd.DataFrame(columns=["ORFpos","Coding","BLAST"])
# per aggiungere una riga in posizione i, inizializzando a False i campi:
# results.loc[i] = [orf.span(), False, False]
# per modificare ad esempio il campo Coding della riga i a True:
# results.loc[i, 'Coding'] = True

#
# SCANSIONE DELLE ORF NELLO STRAND DIRETTO
#
for orf in orf_iter(str(seq)):
  # FLITRAGGIO PRELIMINARE
  results.loc[len(results)] = [orf.span()[0], False, False]
  # INSERIRE QUI LE PROPRIE RIGHE DI CODICE
  pass


# SALVATAGGIO DEI RISULTATI SU FILE
with open(f"group_{group_id}_results.pickle", "wb") as f:
	pickle.dump(results, f)
