"""
class predictor

  fiels -
    recs - dictioary of input fasta sequences (key - fasta header, value - raw sequence)
    new recs - dictionary of inouts after feature selection (key - fasta header, value - raw sequence)

  methods -
    parse fasta - import seqio module from biopython, parse the records in fasta file, store in recs
    remove mini features - dimension reduction by removing U and X form seqs, skipping seqs that are larger than 1000
    encode seq - one hot encoding seqs
    show results - retrieve predicted probabilities as a list
                    import pretty table
                    create output table
                    generate labels using optimal threshold and probabilities
                    append predictions to tabel as fasta header and class

    predict class - take fasta file as input
                    call parse fasta
                    call remove mini features
                    call encode seq
                    call show results
"""
from Bio import SeqIO
from keras.preprocessing.text import hashing_trick  # one hot encode method
from keras.preprocessing.sequence import pad_sequences  # padding method
from prettytable import PrettyTable # display results as a table
import itertools  # iterate through multiple lists once
import numpy as np
from keras.models import load_model # load saved model

class Predictor:
    def __init__(self, file):
        self.file = file
        self.recs = dict()
        self.new_recs = dict()

    def parseFasta(self):
        for record in SeqIO.parse(self.file, 'fasta'):  #passing fasta records
            header = record.id  # fasta header
            accession = header.split("|")[1]  # accession number
            # key = accession
            value = record.seq
            self.recs[accession] = value

    def removeMiniFeatures(self):
        for id, seq in self.recs.items():
            new_seq = ''
            for aa in seq:
                if aa not in ["U", "X"]:
                    new_seq += aa

            if len(new_seq) > 1000:
                print("The recored of ID : ", id, "exceeds the maximum possible length")
                print("\tskipping record ", id)
                continue

            else:
                new_seq = " ".join(new_seq)
                self.new_recs[id] = new_seq

    def encodeSeq(self):
        voc_size = 20 # no. of amino acids
        max_len = 1000  # max possible length
        hashed = [hashing_trick(seq, n=voc_size, hash_function="md5", split=" ") for seq in self.new_recs.values()] # one hot encoded list of seq
        pad_seq = pad_sequences(hashed, max_len, padding="post")  # padded seqs
        return pad_seq

    def showResults(self, labels):
        table = PrettyTable()
        table.field_names = ['Record', 'Class']
        threshold = 0.31482762  # optimal threshold
        for (header, p_val) in zip(self.new_recs.keys(), labels): # iterate through new recs and predicted labels
          if p_val >= threshold:
            fn = "Ubiquitin transferase"
          else:
            fn = "Non-ubiquitin transferase"

          table.add_row([header, fn]) # append row
        print(table)


    def predictClass(self):
        self.parseFasta() # parsing records
        self.removeMiniFeatures() # dimension reduction
        pad_seq = self.encodeSeq()  # one hot encoding
        model = load_model("rat_model.h5")
        predicted = model.predict(pad_seq)
        predicted = np.reshape(predicted, (-1))
        self.showResults(predicted) # generate table


# method to automate the user interaction and prediction
# take file path as input, create Predictor object
# call predictClass method

def runModel():
  path = input("Enter the path for fasta file : ")
  unknown = Predictor(path)
  unknown.predictClass()

runModel()
