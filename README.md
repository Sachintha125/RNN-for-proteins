# RNN-for-proteins
This repository implements a Bi-LSTM (Bidirectional Long Short-Term Memory) model for classifying protein sequences as Ubiquitin-transferases or non-Ubiquitin transferases. Ubiquitin-transferases play a crucial role in various cellular processes, making their identification a valuable tool in biological research.

## Proposed Approach
The model leverages the power of Bi-LSTM networks, which are well-suited for sequential data like protein sequences. Bi-LSTMs can effectively capture long-range dependencies within the amino acid sequence, leading to accurate classifications.

## Software Functionality
This software offers the following capabilities:

 * Classification: Classifies protein sequences (max length 1000 amino acids) as Ubiquitin-transferases or non-Ubiquitin transferases.
  
 * Input: Accepts input in FASTA format, allowing for multiple protein sequences (fasta records) in a single file.
  
 * Feature Selection & Vectorization: Performs feature selection (details provided in code comments) by selecting most criticalamino acids and converts sequences      into vectors suitable for the Bi-LSTM model.
  
## Model: Utilizes a trained Bi-LSTM model for robust classification.

## Running the Software

### Prerequisites:
 * Python 3.9
 * TensorFlow 2.15.0
 * Biopython
 * PrettyTable

### Steps:
 * Clone this repository.
 * Install required libraries.
 * Place your FASTA file(s) in the same directory as the script.
 * Run the Predictor.py using the command line.
  
### Output:
  The script will output the classification results for each protein sequence in the FASTA file(s).
  
## Contact
Feel free to contact me at akalankasachintha125@gmail.com for any questions or collaboration opportunities.
