def read_fasta_file(filename):
    sequences = {}
    with open(filename, 'r') as file:
        current_sequence_name = None
        current_sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence_name is not None:
                    sequences[current_sequence_name] = current_sequence
                current_sequence_name = line[1:]
                current_sequence = ''
            else:
                current_sequence += line
        if current_sequence_name is not None:
            sequences[current_sequence_name] = current_sequence
    return sequences