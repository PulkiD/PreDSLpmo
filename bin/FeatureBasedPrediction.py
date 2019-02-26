#!/usr/bin/env python
#Author: Pulkit Anupam Srivastava
#Version: 1.0
#Last Modified: 26 Feb, 2019
#Description: The script annotates a protein sequence as a member of given LPMO family or not. The prediction method
#implemented in the script is feature based neural network, where descriptor value of protein sequences for each
#fetaure-set is fed into the neural network. Once prediction is made by script for each feature-set, the common
#protein sequences predicted as a member of LPMO family in all feature-set is further labelled as potential LPMO.
#Input: Path to fasta file. Name of the family of LPMO.
#Output: A file with potential LPMO protein sequences.
import numpy as np
import pandas as pd
import tensorflow as tf
import glob
import os
import sys
import subprocess
import copy
from tensorflow import keras
from Bio import SeqIO
from DataCleaning import genCleanFile
#To extract IDs predicted as LPMO
def getID(df):
    ID_List = list()
    for i in range(1,len(df)):
        if (df.iloc[i].Predicted == 0):
            #print ("D")
            ID_List.append(df.iloc[i].ID)
    return ID_List
#To write a file having probability of a sequence to be LPMO or not
def combinePredTest(Y_pred, ID_Test):
    str2ryt = "Predicted\tID\n"
    for i in range(0,len(ID_Test)):
        str2ryt += str(Y_pred[i])+"\t"+ID_Test.iloc[i]+"\n"
    return str2ryt
#To annotate a sequence as member of given LPMO family
def getPrediction(X_Test, ID_Test, path, func, family):
    optimizedModel = os.path.join(sys.path[0],"Models",family,"fbdl",func+'.h5')
    new_model = keras.models.load_model(optimizedModel)
    predictions = new_model.predict(X_Test)
    Y_Pred = np.argmax(predictions, 1)
    #
    PredTest = combinePredTest(Y_Pred, ID_Test)
    PredTestFile = os.path.join(path,func,"FB_Predictions.txt")
    with open(PredTestFile, "w+") as f_prediction:
        f_prediction.write(PredTest)
    f_prediction.close()
#Input to be given while execution
path2file = sys.argv[1]
family = sys.argv[2]
path, input_filename = os.path.split(path2file)
#To filter out protein sequence having unrecognized residues
CleanedFastaFile = genCleanFile(path2file)
#To generate descriptor value of protein sequences for each feature-set 
subprocess.call("Rscript "+os.path.join(sys.path[0],"ExtractFeatures.R ")+path+" "+CleanedFastaFile, shell=True)
#List of feature-set
func2change = ["extractGeary", "extractCTriad","extractDC",
               "extractMoran", "extractMoreauBroto", "extractTC"]
#Iteration over list of feature set for prediction of LPMO family members
for i in func2change:
    CSVFile = os.path.join(path,i,i+".csv")
    df = pd.read_csv(CSVFile,index_col=0, dtype={'label':str})
    df_Test = df.dropna()
    df_Test.reset_index(drop=True, inplace=True)
    X_Test = df_Test.drop(['ID', 'label'], axis = 1)
    ID = df_Test['ID']
    getPrediction(X_Test, ID, path, i, family)
#For extracting common IDs
first_list = second_list = common_list = list()
for i in range(0,len(func2change)):
    CSVFile = os.path.join(path,func2change[i],"FB_Predictions.txt")
    df = pd.read_csv(CSVFile, sep = '\t')
    if (i == 0):
        ID_List_1 = getID(df)
        first_list = copy.copy(ID_List_1)
    else:
        second_list = getID(df)
        common_list = list()
        for j in first_list:
            if j in second_list:
                common_list.append(j)
        first_list = copy.copy(common_list)
#For extracting Fasta sequence using common ID
Sequences = SeqIO.to_dict(SeqIO.parse(CleanedFastaFile, "fasta"))
str2ryt_newFile = ""
for key in common_list:
   str2ryt_newFile += Sequences[key].format("fasta")
outSeqFile = os.path.join(path,"FB_Potential"+family+".fasta")
with open(outSeqFile, "w+") as f_potentialLPMO:
    f_potentialLPMO.write(str2ryt_newFile)
f_potentialLPMO.close()
