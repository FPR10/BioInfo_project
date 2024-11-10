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
import Bio.Blast
from Bio.Blast import NCBIWWW, NCBIXML
import Entropia
from Bio.Seq import Seq

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
sequenza_diretta = next(Bio.SeqIO.parse("group_1_seq.fasta", "fasta")).seq

#SEQUENZA INVERSA COMPLEMENTARE
seq = sequenza_diretta.reverse_complement()

#
# PREPARAZIONE DEL DATAFRAME CON I RISULTATI
#
results = pd.DataFrame(columns=["ORFpos", "Coding", "checkBLAST", "Entropia"])
# per aggiungere una riga in posizione i, inizializzando a False i campi:
# results.loc[i] = [orf.span(), False, False]
# per modificare ad esempio il campo Coding della riga i a True:
# results.loc[i, 'Coding'] = True

########################################################################################################################
#DICHIARAZIONE FUNZIONI UTILI

# SENSORI DI SEGNALE
def tata_box(sequenza, posA):
    regex_tata = r"(TATA)(A|T)(A)(A|T)"
    substring = sequenza[posA-35:posA-2]
    ret = regex.search(regex_tata, substring)
    if ret is not None:
        return True
    return False

def iniziatore(sequenza):
    regex_in = r"(C|T){2}A(A|C|G|T)(A|T)(C|T){2}"
    ret = regex.search(regex_in, sequenza)
    if ret is not None:
        return True
    return False

def posizione_iniziatore(sequenza):
    regex_in = r"(C|T){2}A(A|C|G|T)(A|T)(C|T){2}"
    ret = regex.search(regex_in, sequenza)
    return ret.span()

def sequenza_Kozak(sequenza):
    regex_K = r"(GCC)?(GCC)(A|G)(CC)(ATG)(G)"
    ret = regex.search(regex_K, sequenza)
    if ret is not None:
        return True
    return False

# SENSORI DI CONTENUTO
def cg_content(sequenza):
    numeroC = sequenza.count('C')
    numeroG = sequenza.count('G')
    risultato = (numeroC + numeroG) / len(sequenza)
    return risultato

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

    if contenuto_GC >= 0.5 and CpG_osservato_atteso >= 0.6:
        return True
    return False

    # NB. CONTENUTO CG > 0.5; RAPPORTO CpG OSSERVATO/ATTESO > 0.6


# FUNZIONE PER LA RICERCA BLAST SEQUENZA PER SEQUENZA

def ricerca_blast(sequenza):
    res = []
    # Esegui la ricerca BLAST
    x = NCBIWWW.qblast('blastn', 'nt', sequenza, service='megablast')

    # Parla direttamente la risposta senza salvare su file
    blast_records = NCBIXML.parse(x)

     # Itera sui risultati e raccogli i dati principali
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                res.append({
                    "Sequenza considerata": seq,
                    "Nome allineamento": alignment.title,
                    "Lunghezza allineamento": alignment.length,
                    "E-value": hsp.expect,
                    "Score": hsp.score,
                    "Inizio query": hsp.query_start,
                    "Fine query": hsp.query_end,
                })
    x.close()
    return res


# FUNZIONE PER COSTRUIRE IL DATAFRAME CON I DATI RICAVATI DALLA RICERCA SU BLAST

def costruisci_dataFrame_risultati_blast(sequenza):
    df = pd.DataFrame(ricerca_blast(sequenza))
    return df

df_codice_orf = pd.DataFrame(columns=["ORF", "Indice posizione"])

########################################################################################################################

#DIZIONARIO PER MEMORIZZARE LE ORF E I RELATIVI PUNTEGGI

raccolta_orf = {}

for orf in orf_iter(str(seq)):
    # FILTRAGGIO PRELIMINARE per lunghezza (>250): se la lunghezza della ORF è minore della soglia non la consideriamo
    # altrimenti la aggiunge al risultato
    lunghezza = orf.span()[1] - orf.span()[0]
    if lunghezza > 250:
        raccolta_orf[orf] = 0
        if sequenza_Kozak((str)(orf)):
            raccolta_orf[orf] += 1
        if iniziatore((str)(seq[:orf.span()[0]])):
            raccolta_orf[orf] += 1
            pos = posizione_iniziatore((str)(seq[:orf.span()[0]]))[0]+2
            if tata_box((str)(seq), pos):
                raccolta_orf[orf] +=1
        if isole_CpG((str)(orf)):
            raccolta_orf[orf] += 1

#FILTRAGGIO ORF CON PUNTEGGIO PIÙ ALTO E INSERIMENTO NEL DATAFRAME FINALE

for orf in raccolta_orf:
    if raccolta_orf[orf] > 1:
        df_codice_orf.loc[len(results)] = [orf, len(results) ]
        results.loc[len(results)] = [orf.span()[0], True, False, Entropia.calcola_entropia(str(orf))] #SE HANNO UN PUNTEGGIO > 1, SUPPONIAMO CHE SIANO CODIFICANTI
        if not costruisci_dataFrame_risultati_blast(orf.group()).empty:
            results.loc[len(results)-1, 'checkBLAST'] = True #SE SI TROVANO CORRISPONDENZE SU BLAST, VIENE AGGIORNATO IL PARAMETRO


# SALVATAGGIO DEI RISULTATI SU FILE
with open(f"group_1_results.pickle", "wb") as f:
    pickle.dump(results, f)


#Visualizzazione del dataframe dei risultati
print(df_codice_orf)
print(results)
