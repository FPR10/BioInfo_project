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
import regex
import Bio
import Bio.SeqIO
from bisect import bisect_left
import pickle
import Bio.Blast
from Bio.Blast import NCBIWWW, NCBIXML
import Entropia


# Ricerca delle ORF: orf_iter restituisce un iteratore sulle orf
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

seq = next(Bio.SeqIO.parse(f"group_1_seq.fasta", "fasta")).seq

#DataFrame contente le sole ORF codificanti
results = pd.DataFrame(columns=["ORFpos", "Coding", "checkBLAST", "Entropia"])

#DataFrame contenente tutte le ORF (applicando il solo parametro della lunghezza)
final_results = pd.DataFrame(columns=["ORFpos", "Coding", "BLAST"])



#####################################      FUNZIONI UTILI        #######################################################


# SENSORI DI SEGNALE
def iniziatore_tata(orf):
    sequenza = (str)(seq[orf.span()[0]-300:orf.span()[0]])
    reg = r"(TATA)(A|T)(A)(A|T)(A|T|C|G){25,35}(C|T){2}(A)(A|C|G|T)(A|T)(C|T){2}"
    ret = regex.search(reg,sequenza)
    if ret is not None:
        return True
    return False

#Gli indici fanno riferimento a orf.span()[0]
def sequenza_Kozak(posInizio):
    sequenza = (str)(seq[posInizio-9: posInizio+4])
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

def isole_CpG(orf):
    sequenza = (str)(seq[orf.span()[0]-3000:orf.span()[0]])
    for i in range(len(sequenza)-300):
        subseq = sequenza[i:i+300]
        # Numero CpG osservati (numero di occorrenze CpG nella sequenza
        CpG_osservati = subseq.count("CG")

        # Numero di C E G nella sequenza
        numeroC = subseq.count('C')
        numeroG = subseq.count('G')

        lunghezza = len(subseq)
        if lunghezza == 0:
            return False

        # FREQUENZA ATTESA DI CpG
        # E(CpG) = (numero di C / lunghezza) * (numero di G / lunghezza) * (lunghezza - 1)
        CpG_atteso = (numeroC / lunghezza) * (numeroG / lunghezza) * (lunghezza-1)
        if CpG_atteso == 0:
            return False

        # Rapporto osservato atteso
        CpG_osservato_atteso = CpG_osservati / CpG_atteso

        # CONTENUTO CG
        contenuto_GC = cg_content((str)(orf))

        #CONTENUTO CG > 0.5; RAPPORTO CpG OSSERVATO/ATTESO > 0.6
        if contenuto_GC >= 0.5 and CpG_osservato_atteso >= 0.6: 
            return True
    return False

    

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

#Costruisce un dataframe partendo dai risultati ottenuti dalla ricerca su blast
def costr_dataframe_ris_blast(sequenza):
    df = pd.DataFrame(ricerca_blast(sequenza))
    return df


###############################   ITERAZIONE SU ORF + AGGIORNAMENTO DATAFRAME   #######################################

#DIZIONARIO PER MEMORIZZARE LE ORF E I RELATIVI PUNTEGGI
raccolta_orf = {}

for orf in orf_iter(str(seq)):
    # FILTRAGGIO PRELIMINARE per lunghezza (>250): se la lunghezza della ORF è minore della soglia non la consideriamo
    # altrimenti la aggiunge al risultato
    id_orf = orf.span()[0]
    lunghezza = orf.span()[1] - orf.span()[0]
    if lunghezza > 250:
        final_results.loc[len(final_results)] = [orf.span()[0], False, False]
        raccolta_orf[id_orf] = 0
        
        if sequenza_Kozak(orf.span()[0]):
           raccolta_orf[id_orf] += 1
           
        if iniziatore_tata(orf):
            raccolta_orf[id_orf] += 1
            
        if isole_CpG(orf):
            raccolta_orf[id_orf] += 1
            
        if cg_content((str)(orf)) > 0.45 and cg_content((str)(orf)) < 0.65:
            raccolta_orf[id_orf] += 1
        
        #Aggiorniamo il DataFrame. Se hanno punteggio > 1, supponiamo siano codificanti
        if raccolta_orf[id_orf] > 1:
            results.loc[len(results)] = [id_orf, True, False, Entropia.calcola_entropia(str(orf))] 

            #SE SI TROVANO CORRISPONDENZE SU BLAST, VIENE AGGIORNATO IL PARAMETRO
            if not costr_dataframe_ris_blast(orf.group()).empty:
                results.loc[len(results)-1, 'checkBLAST'] = True


for i in range(len(final_results)):
    for posizione in results['ORFpos']:
        if final_results.loc[i]['ORFpos'] == posizione:
            final_results.loc[i, 'Coding'] = True
            final_results.loc[i, 'BLAST'] = True


# SALVATAGGIO DEI RISULTATI SU FILE
with open(f"group_1_results.pickle", "wb") as f:
    pickle.dump(results, f)