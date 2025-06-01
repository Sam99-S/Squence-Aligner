import argparse
# importing necessary packages
from Bio import Align
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


# defining the class
class SequenceAligner:

    # constructor method to create objects
    def __init__(self, header, sequence, seq_type):
        self.header = header
        self.sequence = sequence
        self.seq_type = seq_type

    # method to get the type of sequence
    def get_seqtype(sequence):

        # checking whether the sequence contains 'M'- if yes -> protein
        if "M" in sequence:
            return "protein"

        # checking whether the sequence contains 'U'- if yes -> RNA
        elif "U" in sequence:
            return "RNA"

        else:
            return "DNA"


    # method to extract sequences from a multi FASTA file
    @staticmethod
    def fasta_split(file):
        sequences = []

        # reading the multi FASTA file
        with open(file, 'r') as f:
            header = ''
            sequence = ''
            for line in f:
                line = line.strip()

                # looking for new FASTA sequences
                if line.startswith('>'):
                    if header:
                        # saving the header, sequence and type in a list
                        seq_type = SequenceAligner.get_seqtype(sequence)
                        sequences.append(SequenceAligner(header, sequence, seq_type))
                        sequence = ''
                    header = line[1:]
                else:
                    sequence += line

                # extracting the information of the last FASTA sequence
            if header and sequence:
                seq_type = SequenceAligner.get_seqtype(sequence)
                sequences.append(SequenceAligner(header, sequence, seq_type))
        return sequences


    # method to calculate the pairwise similarity
    def calculate_global_pairwise_similarity(sequences):
        similarity_scores = {}

        # creating a pairwise alignment instance
        aligner = Align.PairwiseAligner()

        # creating a set of seq types in the file
        unique_seq_types = set(seq.seq_type for seq in sequences)

        # retrieving the sequences of same type from multiple fasta sequences
        for seq_type in unique_seq_types:
            sequences_of_same_type = [seq for seq in sequences if seq.seq_type == seq_type]

            # Checking if sequences of the same type were found
            if sequences_of_same_type:

                # iterate over all combinations of sequences within sequences_of_same_type
                for i, seq1 in enumerate(sequences_of_same_type):
                    for j, seq2 in enumerate(sequences_of_same_type):

                        # avoiding duplicate calculations
                        if i < j:
                            # performing the pairwise global alignment
                            alignments = aligner.align(seq1.sequence, seq2.sequence)

                            # retrieving the alignment score
                            alignment_score = alignments[0].score
                            similarity_scores[(seq1.header, seq2.header)] = alignment_score
        return similarity_scores

    # method to create a dotplot
    @staticmethod
    def create_dotmatrix(seqX, seqY, window_size, threshold):
        # initializing a 2D NumPy array with dimensions equal to the lengths of seqX and seqY
        dp = np.zeros((len(seqX), len(seqY)), dtype=int)

        # iterating over all positions in seqA and seqB
        for i in range(len(seqX)):
            for j in range(len(seqY)):

                # computing the window around position i in seqA and position j in seqB
                windowX = seqX[max(i - window_size, 0):min(i + window_size + 1, len(seqX))]
                windowY = seqY[max(j - window_size, 0):min(j + window_size + 1, len(seqY))]

                # counting the number of matching symbols in the window
                matches = sum([1 for x, y in zip(windowX, windowY) if x == y])

                # setting the matrix element to 1 if the number of matches is at least the specified threshold
                if matches >= threshold:
                    dp[i, j] = 1
        return dp

    # method to display the dotplot
    @staticmethod
    def show_dotplot(seq_dict, accession1, accession2, window_size, threshold):
        # seq_dict = {}
        # retrieving the headers of the complimentary accession number
        headers1 = [header for header in seq_dict if accession1 in header]
        headers2 = [header for header in seq_dict if accession2 in header]

        # extracting the sequences for the alignment
        if len(headers1) == 1 and len(headers2) == 1:
            seq1 = seq_dict[headers1[0]]
            seq2 = seq_dict[headers2[0]]

            # generating a dot matrix of the similarity between the two sequences
            dp = SequenceAligner.create_dotmatrix(seq1, seq2, window_size, threshold)

            # creating a scatter plot for the dot plot
            fig, ax = plt.subplots(figsize=(25, 18))
            ax.set_title('Dot Plot for Pairwise Similarity')

            # plotting dots where the dot plot matrix elements are non-zeros
            rows, cols = np.where(dp)
            ax.scatter(cols, rows, marker='.', color='black')

            # setting labels and grid
            ax.set_xlabel(headers1[0])
            ax.set_ylabel(headers2[0])
            ax.grid(False)

            # displaying the sequence on each axis
            ax.set_xticks(np.arange(len(seq2)))
            ax.set_yticks(np.arange(len(seq1)))
            ax.set_xticklabels(list(seq2))
            ax.set_yticklabels(list(seq1)[::-1])

            # saving the plot as an image
            plt.savefig('dotplot.png')

    # method to run an online BLAST
    def BLAST(file):
        xml_files = []

        # splitting the multi FASTA file
        seq = SequenceAligner.fasta_split(file)
        result_handles = []

        # checking for each sequence type and performing BLAST
        for seqs in seq:
            result_handle = None

            # if the type is DNA -> blastn
            if seqs.seq_type == "DNA":
                result_handle = NCBIWWW.qblast("blastn", "nt", seqs.sequence)
            # if the type is protein -> blastp
            elif seqs.seq_type == "protein":
                result_handle = NCBIWWW.qblast("blastp", "nr", seqs.sequence)

            # each blast result is appended to a list with the respective header
            if result_handle:
                result_handles.append((seqs.header, result_handle))

        # reading the result handles in the list
        for header, result_handle in result_handles:
            blast_result = result_handle.read()
            result_handle.close()

            # creating separate xml file for each sequence to write the blast result
            with open(f"{header}_blast_result.xml", "w") as out_handle:
                out_handle.write(blast_result)

            # maintaining a list of xml files created
            xml_files.append(f"{header}_blast_result.xml")

        # iterating through the list of xml files
        for xml_file in xml_files:
            header = xml_file.split('.')[0]
            with open(xml_file, "r") as handle:

                # reading the output blast hits file
                blast_records = NCBIXML.parse(handle)
                best_hit = None

                # setting the e value to filter the blast results
                e_value_threshold = 0.05

                # creating a txt file to store the results
                with open(f"{header}.txt", "w") as out:

                    # iterating through the blast hit records and extracting information
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if hsp.expect < e_value_threshold:
                                    out.write(f'Title : {alignment.title}\n')
                                    out.write(f'Length : {alignment.length}\n')
                                    out.write(f'E-value : {hsp.expect}\n')
                                    out.write(f'Score : {hsp.score}\n')
                                    out.write(f'Subject sequence : {hsp.sbjct}\n')
                                    out.write(f'Length of the hit : {len(hsp.sbjct)}\n')
                                    out.write('\n\n')

                                # finding the best hit by its lowest e value
                                if best_hit is None or hsp.expect < best_hit.expect:
                                    best_hit = hsp

                # output the information of the best hit
                print("\n- Best Hit for ", header, " -")
                if best_hit:
                    print('Title:', alignment.title)
                    print('E-value:', best_hit.expect)
                    print('Score:', best_hit.score)
                    print('Query sequence:', best_hit.query)
                    print('Subject sequence:', best_hit.sbjct)
                    print('Length of the query sequence:', len(best_hit.query))
                    print('Length of the hit:', len(best_hit.sbjct))

    # method to retrieve a sequence when the accession number is given by the user
    @staticmethod
    def retrieve_sequence():

        # defining user inputs
        accession = input("\nEnter the accession number : ")
        seqType = input("Enter the sequence type : ")
        Entrez.email = "sam@gmail.com"

        # specifying which database to be used to retrieve the sequence
        if seqType in {"DNA", "dna", "RNA", "rna"}:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        elif seqType == "protein" or "Protein":
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")

        # reading the result and writing the sequence in a new fasta file
        record = SeqIO.read(handle, "fasta")
        with open(f"{accession}.fasta", "w") as out_handle:
            SeqIO.write(record, out_handle, "fasta")
        print("FASTA file created successfully.")

    # method to calculate the melting temperature of a nucleotide sequence
    def calculate_melting_temp(sequence, seq_type):

        # checking whether the sequence is of nucleotides and calculating the melting temp
        if seq_type == "DNA" or seq_type == "RNA":
            mt = MeltingTemp.Tm_NN(sequence)
            return mt
        else:
            return None

    # method to calculate the isoelectric point of a protein
    def calculate_isoelectric_point(sequence, seq_type):

        # checking whether the sequence is of amino acids and calculating the isoelectric point
        if seq_type == "protein":
            P_ie = ProteinAnalysis(str(sequence))
            ie = P_ie.isoelectric_point()
            return ie
        else:
            return None

    # method to calculate the aromaticity of a protein
    def calculate_aromaticity(sequence, seq_type):

        # checking whether the sequence is of amino acids and calculating the aromaticity
        if seq_type == "protein":
            P_aro = ProteinAnalysis(str(sequence))
            aro = P_aro.aromaticity()
            return aro
        else:
            return None


def analyze_fasta(sequences):
  """
  Provides options for further analysis after parsing a FASTA file.
  """
  while True:
    print("\nFASTA Analysis Options:")
    print("1. Get Sequence Type")
    print("2. Calculate Pairwise Similarity Scores")
    print("3. Generate Pairwise Similarity Dotplot")
    print("4. Back to Main Menu")

    choice = input("Enter your choice (1-4): ")

    if not choice.isdigit() or int(choice) < 1 or int(choice) > 4:
      print("Invalid choice. Please enter a number between 1 and 4.")
      continue

    choice = int(choice)

    if choice == 1:
        print("\n")
        for seq in sequences:
            print(seq.header + " : " + seq.seq_type)

    elif choice == 2:
        similarity_scores = SequenceAligner.calculate_global_pairwise_similarity(sequences)
        with open("Pairwise_similarity_scores.txt", 'w') as file:
            for (header1, header2), score in similarity_scores.items():
                file.write(f"{header1} \t&\t {header2} \t:\t Similarity Score = {score}\n")

    elif choice == 3:
        try:
                seq_dict = {seq.header: seq.sequence for seq in sequences}
                accession1 = input("Enter the first accession number: ")
                accession2 = input("Enter the second accession number: ")
                window_size = int(input("Enter a suitable window_size (nucleotide seq -> 5-20, amino acid seq -> 3-5): "))
                threshold = int(float(input("Enter a suitable threshold value (nucleotide seq -> 7-9, amino acid seq -> 3-5): ")))
                SequenceAligner.show_dotplot(seq_dict, accession1, accession2, window_size, threshold)
        except Exception as e:
                print("\nCouldn't create the dotplot:", str(e))

    elif choice == 4:
      break

def analyze_sequence(sequence, seq_type):
    while True:
        print("\nSequence Analysis Options:")
        print("1. Calculate Melting Temperature(DNA)")
        print("2. Calculate Isoelectric Point(Protein)")
        print("3. Calculate Aromaticity(Protein)")
        print("4. Back to Main Menu")

        choice = input("Enter your choice (1-4): ")

        if not choice.isdigit() or int(choice) < 1 or int(choice) > 4:
          print("Invalid choice. Please enter a number between 1 and 4.")
          continue

        choice = int(choice)

        if choice == 1 and seq_type not in ["dna", "DNA", "rna", "RNA"]:
              print("Error: Melting temperature calculation is only valid for DNA/RNA sequences.")
              continue
        elif (choice == 2 or choice == 3) and seq_type != "protein":
              print("Error: Isoelectric point and aromaticity calculations are only valid for protein sequences.")
              continue

        if choice == 1:
              mt = SequenceAligner.calculate_melting_temp(sequence, seq_type)
              print(f"Melting Temperature: {mt:.2f}")
        elif choice == 2:
              iso = SequenceAligner.calculate_isoelectric_point(sequence, seq_type)
              print(f"Isoelectric Point: {iso:.2f}")
        elif choice == 3:
              aro = SequenceAligner.calculate_aromaticity(sequence, seq_type)
              print(f"Aromaticity: {aro:.2f}")
        else:
          break

def analyze_sequences(sequences):
    while True:
        print("\nSequence Analysis Options:")
        print("1. Calculate Melting Temperature(DNA)")
        print("2. Calculate Isoelectric Point(Protein)")
        print("3. Calculate Aromaticity(Protein)")
        print("4. Back to Main Menu")

        choice = input("Enter your choice (1-4): ")

        if not choice.isdigit() or int(choice) < 1 or int(choice) > 4:
          print("Invalid choice. Please enter a number between 1 and 4.")
          continue

        choice = int(choice)

        for seqs in sequences:
            if choice == 1:
                mt = SequenceAligner.calculate_melting_temp(seqs.sequence, seqs.seq_type)
                if mt is not None:
                    print(f"{seqs.header} : Melting temperature = %.2f " % mt)
            elif choice == 2:
                iso = SequenceAligner.calculate_isoelectric_point(seqs.sequence, seqs.seq_type)
                if iso is not None:
                    print(f"{seqs.header} : Isoelectric point = %.2f " % iso)
            elif choice == 3:
                aro = SequenceAligner.calculate_aromaticity(seqs.sequence, seqs.seq_type)
                if aro is not None:
                    print(f"{seqs.header} : Aromaticity = %.2f " % aro)

        if choice == 4:
            break


def main():
  print("Welcome to the Sequence Aligner Tool!")

  while True:
    # User menu
    print("\nSelect an option:")
    print("1. Parse and Analyze FASTA File")
    print("2. BLAST Search")
    print("3. Retrieve Sequence by Accession Number")
    print("4. Sequence Analysis")
    print("5. Exit")

    choice = input("Enter your choice (1-5): ")

    # Validate user input
    if not choice.isdigit() or int(choice) < 1 or int(choice) > 5:
      print("Invalid choice. Please enter a number between 1 and 5.")
      continue

    choice = int(choice)

    if choice == 1:
      # Parse and Analyze FASTA File
      fasta_file = input("Enter the FASTA file path: ")
      try:
        # Check if file exists
        with open(fasta_file, 'r') as f:
          pass
        sequences = SequenceAligner.fasta_split(fasta_file)
        print("\nFASTA file parsed successfully!")
        analyze_fasta(sequences)
      except FileNotFoundError:
        print("Error: File not found. Please provide a valid FASTA file path.")

    elif choice == 2:
      # BLAST Search
      fasta_file = input("Enter the FASTA file path: ")
      try:
        # Check if file exists
        with open(fasta_file, 'r') as f:
          pass
        SequenceAligner.BLAST(fasta_file)
        print("\nBLAST search completed. Results saved in XML files.")
      except FileNotFoundError:
        print("Error: File not found. Please provide a valid FASTA file path.")

    elif choice == 3:
      # Retrieve Sequence by Accession Number
      SequenceAligner.retrieve_sequence()

    elif choice == 4:
      # Sequence Analysis
      fasta_file = input("Enter the FASTA file path (or type '>' to use sequence input): ")
      if fasta_file == '>':
        sequence = input("Enter your sequence: ")
        seq_type = input("Enter the sequence type (DNA/RNA/protein): ")

        if seq_type not in ["dna", "DNA", "rna", "RNA", "protein"]:
            print("Error: Invalid sequence type. Please enter DNA, RNA, or protein.")
            continue
        analyze_sequence(sequence, seq_type)

      else:
          try:
            # Check if file exists
            with open(fasta_file, 'r') as f:
              pass
            sequences = SequenceAligner.fasta_split(fasta_file)
            analyze_sequences(sequences)
          except FileNotFoundError:
            print("Error: File not found. Please provide a valid FASTA file path.")

    elif choice == 5:
      break



if __name__ == "__main__":
    main()
# if __name__ == '__main__':
#
#     # creating objects of multiple fasta sequences in a single file
#     sequences = SequenceAligner.fasta_split("SAligner.fasta")
#
#     # # displaying the type of each sequence
#     # for seq in sequences:
#     #     # print(seq.header)
#     #     # print(seq.sequence)
#     #     print(seq.header + " : " + seq.seq_type)
#     #
#     # # calculating the pairwise similarity score and displaying in a new text file
#     # similarity_scores = SequenceAligner.calculate_global_pairwise_similarity(sequences)
#     # with open("Pairwise_similarity_scores.txt", 'w') as file:
#     #     for (header1, header2), score in similarity_scores.items():
#     #         file.write(f"{header1} \t&\t {header2} \t:\t Similarity Score = {score}\n")
#     #
#     sequence = SequenceAligner.fasta_split("SAligner.fasta")
#     # calculating the melting temp, isoelectric point, aromaticity of the sequences
#     for seq in sequence:
#         mt = SequenceAligner.calculate_melting_temp(seq.sequence, seq.seq_type)
#         if mt is not None:
#             print(f"\n{seq.header} : Melting temperature = %.2f " % mt)
#         ie = SequenceAligner.calculate_isoelectric_point(seq.sequence, seq.seq_type)
#         if ie is not None:
#             print(f"{seq.header} : Isoelectric point = %.2f " % ie)
#         aro = SequenceAligner.calculate_aromaticity(seq.sequence, seq.seq_type)
#         if aro is not None:
#             print(f"{seq.header} : Aromaticity = %.2f " % aro)
#     #
#     # # print(SequenceAligner.calculate_melting_temp("attgctagggcatta","DNA"))
#     print(SequenceAligner.calculate_isoelectric_point("LCGSHLVEALYLVCGERGFFYTPMSRREVED","protein"))
#     print(SequenceAligner.calculate_aromaticity("LCGSHLVEALYLVCGERGFFYTPMSRREVED","protein"))
#
#     sequences_dict = {}
#     for seqs in sequences:
#         sequences_dict[seqs.header] = seqs.sequence
#     #
#     # # # calling the dotplot creation method
#     # SequenceAligner.show_dotplot(sequences_dict, "NM_000515.5", "NM_008117.3", 5, 8.5)
#     SequenceAligner.show_dotplot(sequences_dict, "P01308", "P01326", 1, 3)
#     #
#     # # calling the BLAST performing method
#     # f = "SAligner2.fasta"
#     # SequenceAligner.BLAST(f)
#
#     # calling the sequence retrieval method
#     # SequenceAligner.retrieve_sequence()
