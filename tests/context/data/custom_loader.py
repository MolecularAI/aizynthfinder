def extract_smiles(filename):
    with open(filename, "r") as fileobj:
        for i, line in enumerate(fileobj.readlines()):
            if i == 0:
                continue
            yield line.strip().split(",")[0]
