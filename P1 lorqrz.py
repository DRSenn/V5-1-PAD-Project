def ParseSeqFile(filename):
    """Read and validate a file with nucleotide sequences."""

    sequences = []

    with open(filename, "r") as file: #open the file in reading mode / with closes the file
        for line in file: #read one line at a time

            if not line.strip(): #move to the next line if line is empty or contains white spaces
                continue #back to the beginning of the loop

            if line[0] != '>':
                raise ValueError("malformed input")

            parts = line[1:].split() #remove >

            if len(parts) < 2:
                raise ValueError("malformed input")

            label = parts[0] #first element is the label
            sequence = "" #each new label gets a new sequence

            for fragment in parts[1:]:
                sequence = sequence + fragment

            for base in sequence:
                if base.upper() not in 'ACGT':
                    raise ValueError('malformed input')

            sequences.append((label, sequence)) #append one pair to the list during each interaction of the outer loop

    return sequences #sequences after entire file has been processed successfully

print(ParseSeqFile('sequences.txt'))