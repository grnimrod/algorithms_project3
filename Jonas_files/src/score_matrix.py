def initiate_score_matrix(file):
    raw = open(file, "r").read()
    lines_sep = raw.split('\n')
    nested_dict = {}

    for line in lines_sep:
        parsed_line = line.split()
        key = parsed_line[0]
        scores_per_nucl = {}

        for char in range(1, len(parsed_line)):
            # Defining the inner dictionary
            scores_per_nucl[lines_sep[char-1][0]] = int(parsed_line[char])
        
        nested_dict[key] = scores_per_nucl
    
    return nested_dict