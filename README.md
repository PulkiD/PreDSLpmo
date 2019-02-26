# PreDSLpmo
A neural network based prediction tool for functional annotation of Lytic Polysaccharide Monooxygenases

# Requirements
  - Python v3.5 or above
    - Tensorflow v1.12.0
    - Pandas v0.23.4
    - Numpy v1.15.4
    - Biopython v1.73
    - Keras v2.2.4
  - R v3.4.1 or above
    - ProtR

# Description
  **Motivation**:
  Lytic polysaccharide monooxygenases (LPMO), a family of copper-dependent oxidative enzymes, boost the degradation of crystalline
  polysaccharides, such as cellulose and chitin. While multiple experimental methods are used for accurate identification of LPMOs, a
  computational method that can aid experimental methods through fast and accurate classification of sequences into LPMOs families,
  namely AA9, AA10, AA11, AA13, AA14, and AA15, would be an important step towards production of enzymatic mixtures adept at efficiently
  degrading recalcitrant polysaccharides.
  
  **Framework**:
  
  ![Framework of the NN implemented in both approaches.](https://github.com/PulkiD/PreDSLpmo/blob/master/Images/Figure_1.tif)
  
  A) Input layer that was composed of a number of neurons equal to the number of descriptors in each features-set. Finally,  two hidden
  layers with ReLU activation function and an output layer composed of two neurons with sigmoid activation function was implemented to
  built a traditional NN.
  B) The padded sequence was then fed into the embedding layer composed of equal number of neurons as the padding length. The output
  from embedding layer was then fed to a bi-directional long short term memory units which was further connected to fully-connected
  dense layer composed of 100 and 50 neurons in each layer. The output layer had two neurons that provide the probability of sequence to
  belong in either the AA9/AA10 class or the Non-AA9/NonAA10 class.
  
  **Result**:
  In this study, we developed a machine learning based tool called PreDSLpmo that employs two different approaches to functionally
  classify a given protein sequence(s) into major LPMO families (AA9 and AA10).The first approach uses traditional neural network (NN),
  while the second approach employs bi-directional long short-term memory for sequence classification. Our method shows significant
  improvement in predictive power when compared with dbCAN2, an existing HMM-profile-based CAZyme predicting tool, on both validation
  and independent benchmark sets.
  
  ![Results.](https://github.com/PulkiD/PreDSLpmo/blob/master/Images/Figure_2.tif)

# File Description
  - bin
    - Models: This folder contains two subfolders "AA9" and "AA10" which consist of best/optimized models for both feature based neural
    network and long short-term based neural network.
    - Scripts: The files are named on two coding schemes, namely FeatureBasedPrediction and LSTMBasedPrediction, are python code that
    can be executed to get potential LPMO sequences.
  - Example
    - Ex1.fasta: An example file (protein sequences) whose path has to be provided while executing the script. This file contains
    sequences from which user wants to predict potential LPMOs.

# Example
Command to execute the script:

*python3 /path/to/script/ScriptName.py /path/to/fastfile/FastaFile.fasta FamilyYouWantToAnnotate*

**Eg**: *python3 /bin/FeatureBasedPrediction.py /Exmaple/Ex1.fasta AA9*
  
Output File:

*FB_PotentialAA9.fasta*: Among all the files generated by the script, FB_PotentialAA9.fasta conatins protein sequences that are annotated to be potential member of given LPMO family.
