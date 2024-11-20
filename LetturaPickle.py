import pickle


#APERTURA FILE PICKLE PER VISUALIZZARE I RISULTATI SULLA SEQ. PRINCIPALE
with open ('group_1_results.pickle', 'rb') as file:
    dati = pickle.load(file)
    
print (dati)



#APERTURA FILE PICKLE PER VISUALIZZARE I RISULTATI SULLA SEQ. INVERSA COMPLEMENTARE
with open ('group_1_results_seq_inversa.pickle', 'rb') as file:
    dati = pickle.load(file)
    
print (dati)
