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
  # FILTRAGGIO PRELIMINARE per lunghezza (>150): se la lunghezza della ORF è minore della soglia non la consideriamo
  #altrimenti la aggiunge al risultato
  lunghezza = orf.span()[1] - orf.span()[0]
  if lunghezza > 250:
	  results.loc[len(results)] = [orf.span()[0], False, False]
pass


# SALVATAGGIO DEI RISULTATI SU FILE
with open(f"group_{group_id}_results.pickle", "wb") as f:
	pickle.dump(results, f)


###########################################################################################################################################################################################

# INSERIRE QUI LE PROPRIE RIGHE DI CODICE

# SENSORI DI SEGNALE
def tata_box(sequenza):
    regex_tata = r"(TATA)(A|T)(A)(A|T)"
    ret = regex.search(regex_tata, sequenza)
    return ret.span()

def iniziatore(sequenza):
    regex_in = r"(C|T){2}A(A|C|G|T)(A|T)(C|T){2}"
    ret = regex.search(regex_in, sequenza)
    return ret.span()

def isole_CpG(sequenza):
    # NUMERO CpG OSSERVATI (numero di occorrenze CpG nella sequenza
    CpG_osservati = len(regex.findall(r"CG", sequenza))

    # NUMERO DI C E G NELLA SEQUENZA
    numeroC = sequenza.count('C')
    numeroG = sequenza.count('G')

    lunghezza = len(sequenza)

    # FREQUENZA ATTESA DI CpG
    # E(CpG) = (numero di C / lunghezza) * (numero di G / lunghezza) * (lunghezza - 1)
    CpG_atteso = (numeroC / lunghezza) * (numeroG / lunghezza) * (lunghezza - 1)

    # RAPPORTO OSSERVATO/ATTESO
    CpG_osservato_atteso = CpG_osservati / CpG_atteso

    # CONTENUTO CG
    contenuto_GC = cg_content(sequenza)

    return CpG_osservato_atteso, contenuto_GC

    # NB. CONTENUTO CG > 0.5; RAPPORTO CpG OSSERVATO/ATTESO > 0.6

def sequenza_Kozak(sequenza):
    regex_K = r"(GCC)?(GCC)(A|G)(CC)(ATG)(G)"
    ret = regex.search(regex_K, sequenza)
    return ret.span()

# SENSORI DI CONTENUTO
def cg_content(sequenza):
    numeroC = sequenza.count('C')
    numeroG = sequenza.count('G')
    risultato = (numeroC + numeroG) / len(sequenza)
    return risultato

# PROVE
# gen = orf_iter(str(seq))
# print(next(gen))
# print(next(gen))

# print(results)
prova = "TATAAATCGTAATCTAGTCC"
print(iniziatore(str(prova)))

prova_CpG = "ATGCGATACGCGTAA"
print(isole_CpG(sequenza_dna))

prova_K1 = "CATGCCGCCGCCATGGTTT"
print(sequenza_Kozak(prova_K1))

prova_K2 = "ATCGCCACCATGGATG"
print(sequenza_Kozak(prova_K2))
