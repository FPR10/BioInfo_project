import math

def contaOcc(x, stringa):
    return stringa.count(x)


#Calcolo della probabilità dell'evento
def prob_carattere_singolo(carattere, seq):
    """_summary_

    Args:
        carattere (_type_): base azotata
        seq (_type_): sequenza considerata

    Returns:
        out: probabilità della singola base azotata
    """
    return contaOcc(carattere, seq)/len(seq)


def prob_seq_nucleotidica (seq):
    """_summary_

    Args:
        seq (_type_): sequenza 

    Returns:
        _type_: array contente le probabilità di tutte e 4 le basi
    """
    letters = ['A','T','C','G']
    ret = []
    for i in letters:
        ret.append (prob_carattere_singolo(i,seq))
    return ret # ret = [P(A), P(T), P(C), P(G)]
        
        
#devo passare come parametro l'output del metodo "prob_seq_nucleotidica"
def contInformativo (arrayProb):
    ret = [] # ret = [I(A), I(T), I(C), I (G)]
    for elem in arrayProb:
        ret.append (-math.log(elem,2))
    return ret
    
'''
#arrayProb e arrayContInfo possono essere rispettivamente presi
def entropia (arrayProb, arrayContInfo):
    res = 0
    for i in range (len (arrayProb)):
        for j in range (len (contInformativo)):
            if i == j:
                res += arrayProb[i]*arrayContInfo[j]
    return res
'''
#Riceve la sequenza e calcola arrayProb e arrayContInfo internamente. Possiamo così invocare il metodo "entropia" direttamente sulla orf
#senza preoccuparci di effettuare operazioni di calcolo preliminari
def entropia (sequenza):
    arrayProb = prob_seq_nucleotidica(sequenza)
    arrayContInfo = contInformativo (arrayProb)
    res = 0
    for i in range (len (arrayProb)):
        for j in range (len (arrayContInfo)):
            if i == j:
                res += arrayProb[i]*arrayContInfo[j]
    out = "L'entropia della sequenza è: " + str(res) 
    return out


#Prova (esempi presi da capitolo 10)
sequenza1 = "ATGCATGCATGCATGCTTTTGGGGCCCC"
print (entropia(sequenza1))

sequenza2 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGC"
print (entropia(sequenza2))