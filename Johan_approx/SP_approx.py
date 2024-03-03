import os
import sys
import subprocess
import tempfile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import combinations

# Define function to read FASTA file
def read_fasta_file(file_path):
    """
    Reads a FASTA file and returns a dictionary of sequences.
    """
    sequences = {}
    with open(file_path, "r") as fasta_file:
        for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), start=1):
            sequences[f"seq{i}"] = str(record.seq)
    return sequences



# Define functions to generate pairwise combinations and run alignment script
def generate_pairwise_combinations(sequences):
    """
    Generates all possible pairwise combinations of sequences.
    """
    pairwise_combinations = list(combinations(sequences, 2))
    return pairwise_combinations

def run_alignment_script_score(fasta_file, matrix_file):
    """
    Runs alignment script with the given FASTA file and returns the alignment score.
    """
    command = f"python global_align_linear.py {fasta_file} -g 5 -m {matrix_file} --hide-alignments"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode())
        return None
    output = stdout.decode()
    # Extract alignment score from output
    alignment_score = float(output.strip().split()[-1])
    return alignment_score

def run_alignment_script_backtrack(fasta_file, matrix_file):
    """
    Runs alignment script with the given FASTA file and returns the aligned sequences.
    """
    command = f"python global_align_linear.py {fasta_file} -g 5 -m {matrix_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode())
        return None
    output = stdout.decode()
    # Extract alignment from output
    alignment_str = output.strip().split('\n')[1]  # Extracting the aligned sequences only
    # Remove unnecessary characters
    alignment_str = alignment_str.strip("[").strip("]").replace("'", "")
    # Convert the string to a list by splitting on comma and space
    alignment_list = alignment_str.split(', ')
    return alignment_list

def process_alignment_result(alignment_result):
    """
    Transposes the alignment result.
    """
    alignment_transposed = [list(position) for position in zip(*alignment_result)]
    return alignment_transposed

def extend_msa(alignment):
    MA = alignment[0]
    for sublist in alignment[1:]:
        command = ["python", "extend_msa.py", str(MA), str(sublist)]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            print(stderr.decode())
            return None
        MA = eval(stdout.decode())  # Update MA with the result from the script
    return MA
# Main function
def main():
    # Read sequences from input FASTA file
    fasta_file_path = sys.argv[1]
    sequences = read_fasta_file(fasta_file_path)

    # Create a dictionary to store sequences with initial values of 0
    sequence_scores = {seq_id: 0 for seq_id in sequences}
    alignment=[]
    # Generate pairwise combinations
    pairwise_combinations = generate_pairwise_combinations(sequences)

    # Create temporary directory for temporary FASTA files
    temp_dir = "temp"
    os.makedirs(temp_dir, exist_ok=True)

    
  
    
    # Get path to score matrix file
    matrix_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "substitution_matrix_phylip.txt")
    


    # Iterate through pairwise combinations
    for idx, (seq_id1, seq_id2) in enumerate(pairwise_combinations):
    # Create temporary FASTA file for current pair
        temp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta", dir=temp_dir, delete=False)

    # Adjust sequence IDs in the temporary file
        temp_fasta.write(f">seq{1}\n{sequences[seq_id1]}\n".encode())
        temp_fasta.write(f">seq{2}\n{sequences[seq_id2]}\n".encode())
        temp_fasta.close()
    

    
    # Run alignment script on temporary FASTA file to get alignment score 
        alignment_score = run_alignment_script_score(temp_fasta.name, matrix_file)
        sequence_scores[seq_id2]+= alignment_score
        sequence_scores[seq_id1]+= alignment_score
      
 
    
    # Remove temporary FASTA file
        os.remove(temp_fasta.name)
    
    ## Output alignment scores
    
    for i in sequences:
        if i !=min(sequence_scores):
            # Create temporary FASTA file for current pair
            temp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta", dir=temp_dir, delete=False)

            # Adjust sequence IDs in the temporary file
            temp_fasta.write(f">seq{1}\n{sequences[min(sequence_scores)]}\n".encode())
            temp_fasta.write(f">seq{2}\n{sequences[i]}\n".encode())
            temp_fasta.close()


            
            # Run alignment script on temporary FASTA file to get backtrack
            alignment_result = run_alignment_script_backtrack(temp_fasta.name, matrix_file)
         

                    # Transpose the alignment result
            alignment_transposed = process_alignment_result(alignment_result)

            # Append the transposed alignment to the alignment list
            alignment.append(alignment_transposed)
            
 
            
                # Remove temporary FASTA file
            os.remove(temp_fasta.name)
    #print(alignment)
    extracted_lists = extend_msa(alignment)
    print(extracted_lists)

    
    

    


    




if __name__ == "__main__":
    main()


