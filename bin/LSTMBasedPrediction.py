#!/usr/bin/env python
#Author: Pulkit Anupam Srivastava
#Version: 1.0
#Last Modified: 26 Feb, 2019
#Description: The script annotates a protein sequence as a member of given LPMO family or not. The prediction method
#implemented in the script is long short-term based neural network.
#Input: Path to fasta file. Name of the family of LPMO.
#Output: A file with potential LPMO protein sequences.
from DataCleaning import genCleanFile
from tensorflow import keras
from keras.preprocessing.sequence import pad_sequences
from keras.preprocessing.text import Tokenizer
from Bio import SeqIO
import pickle
import glob
import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
#Pre-processing of fasta file for prediction
def PrepareFile(path, FastaFile):
    str2ryt = "ID,Sequence\n"
    Sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))
    for key in Sequences:
        seq = Sequences[key].seq
        str2ryt+= str(key)+","+str(seq)+"\n"
    #Writing the pre-processed csv file
    PreparedFile = os.path.join(path, "PreparedFile.csv")
    with open(PreparedFile, "w+") as f_prepared:
        f_prepared.write((str2ryt.strip()))
    f_prepared.close()
    return PreparedFile
#To write a file having probability of a sequence to be LPMO or not
def combinePredTest(Y_Pred, ID_Test, predictions):
    str2ryt="ID\tConfidence\n"
    for i in range(len(Y_Pred)):
        if Y_Pred[i] == 1:
            str2ryt+=ID_Test.iloc[i]+"\t"+str(predictions[i][1])+"\n"
    return ((str2ryt.strip()))
#To annotate a sequence as member of given LPMO family
def getPrediction(X_Test, ID_Test, path, family):
    optimizedModel = os.path.join(sys.path[0],'Models',family,'lstm','BestModel.h5')
    new_model = keras.models.load_model(optimizedModel)
    #
    TokenizerFile = os.path.join(sys.path[0],'Models',family,'lstm','tokenizer.pickle')
    with open(TokenizerFile, 'rb') as handle:
        tokenizer = pickle.load(handle)
    X_Test = tokenizer.texts_to_sequences(X_Test)
    if ("AA9" in family):
        X_Test = pad_sequences(X_Test, maxlen=300)
    elif ("AA10" in family):
        X_Test = pad_sequences(X_Test, maxlen=350)
    #
    predictions = new_model.predict(X_Test)
    Y_Pred = np.argmax(predictions, 1)
    #
    PredTest = combinePredTest(Y_Pred, ID_Test, predictions)
    PredTestFile = os.path.join(path,"LSTM_Predictions.txt")
    with open(PredTestFile, "w+") as f_prediction:
        f_prediction.write(PredTest)
    f_prediction.close()
    return PredTest
#Input to be given while execution
path2file = sys.argv[1]
family = sys.argv[2]
path, input_filename = os.path.split(path2file)
#To filter out protein sequence having unrecognized residues
CleanedFastaFile = genCleanFile(path2file)
#For Prediciton
TestFileData = pd.read_csv(PrepareFile(path, CleanedFastaFile))
X_Test = TestFileData['Sequence']
ID_Test = TestFileData['ID']
common_list = getPrediction(X_Test, ID_Test, path, family)
#For extracting Fasta sequence using ID
Sequences = SeqIO.to_dict(SeqIO.parse(CleanedFastaFile, "fasta"))
ID_List = pd.read_csv(os.path.join(path,"LSTM_Predictions.txt"), sep = "\t")
common_list = ID_List['ID']
str2ryt_newFile = ""
for key in common_list:
   str2ryt_newFile += Sequences[key].format("fasta")
outSeqFile = os.path.join(path,"LSTM_Potential"+family+".fasta")
with open(outSeqFile, "w+") as f_potentialLPMO:
    f_potentialLPMO.write(str2ryt_newFile)
f_potentialLPMO.close()
