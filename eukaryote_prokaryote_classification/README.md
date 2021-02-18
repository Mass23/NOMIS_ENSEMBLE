# NOMIS repository for classifying prokaryotes and eukaryotes from assemblies
- Using machine learning models to train the classification of pro- and eukaryotes
- Current status includes the following:
	- created fasta files for refseq sequences from prokaryotes, eukaryotes and viruses
	- concatenated all the fasta sequences, and removed the new lines, along with all the headers
	- split the sequences at 50 Kbp in each file
	- running the 'kmer' counts for each of these .fasta files
	- training the ml_model on TNFs. To be done: GCs and other features
