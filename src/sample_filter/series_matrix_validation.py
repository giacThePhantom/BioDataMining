import json,string

# TCGA is not included here
datasets = ["GSE4107", "GSE4183", "GSE9348", "GSE13067", "GSE13294", "GSE14333", "GSE15960", "GSE17536", "GSE17537", "GSE18105", "GSE18088", "GSE20916", "GSE23878", "GSE33113", "GSE26682", "GSE8671", "GSE21510", "GSE22598", "GSE23194", "GSE32323", "GSE37364", "GSE157982", "GSE38832", "GSE31595", "GSE39084", "GSE39582", "GSE103479"]
datasets.remove("GSE157982") # this is the vit-D signature, the series matrix file has no expression data for some reason

datasets.remove("GSE26682")       # this dataset has actually 2 series matrix files instead of 1
datasets.append("GSE26682-GPL96")
datasets.append("GSE26682-GPL570")

# these are the datasets that have a different set of genes and/or NA values, don't consider them
for d in ["GSE20916", "GSE103479", "GSE26682-GPL96", "GSE26682-GPL570"]: datasets.remove(d)

################################################################################


# Folder with the GEO series_matrix files vvv
GEO_FOLDER = "G:\\CulebraExt\\geo\\"


################################################################################
def json_read(filePath):
    file = open(filePath, "r")
    fileInfo = file.read(); file.close()
    return json.loads(fileInfo)


# stores in a new group of json files the tags (i.e. genes) from each series_matrix
def extract_row_tags():
    for dataset in datasets:
        known_rows = []
        with open(f"{GEO_FOLDER}{dataset}_series_matrix.txt", "r") as file:
            for row in file:
                if row[:26] == "!series_matrix_table_begin":
                    known_rows.clear()
                else:
                    gene = ""
                    for char in row:
                        if char == "\t": break
                        gene += char
                    known_rows.append(gene)
        with open(f"{GEO_FOLDER}genes\\{dataset}.json", "w") as file:
            file.write(json.dumps(known_rows[1:-1]))
        print(dataset)


# compares the gene names extracted by extract_row_tags and displays the differences
# Result: ["GSE20916", "GSE103479", "GSE26682-GPL96", "GSE26682-GPL570"]
def compare_tags():
    folder = GEO_FOLDER + "genes\\"
    alias = {d:string.ascii_letters[i] for i,d in enumerate(datasets)}
    diffs = {}
    for i,d0_name in enumerate(datasets[:-1]):
        d0 = set(json_read(f"{folder}{d0_name}.json"))

        for d1_name in datasets[i + 1:]:
            d1 = set(json_read(f"{folder}{d1_name}.json"))

            if len(d0 | d1) != len(d0 & d1):
                try:                 diffs[alias[d1_name]] += alias[d0_name]
                except KeyError:
                    try:             diffs[alias[d0_name]] += alias[d1_name]
                    except KeyError: diffs[alias[d1_name]]  = alias[d0_name]

    print(*diffs.items(),sep="\n")
    for k,v in alias.items():
        if v in diffs: print(f"{v}: {k}")

# detects non-numeric values inside series_matrix files' gene expression data
def rough_NA_detection():
    for dataset in datasets:
        active = False
        with open(f"{GEO_FOLDER}{dataset}_series_matrix.txt", "r") as file:
            for row in file:
                if row[:26] == "!series_matrix_table_begin":
                    active = True
                elif active:
                    elements = row.split("\t")
                    if elements[0] != "\"ID_REF\"":
                        for e in elements[1:]:
                            try: float(e)
                            except ValueError:
                                print(f"/// '{e}' found in dataset {dataset}")
                                active = False; break

        print(f">>> Dataset {dataset} finished")

# specifically verifies GSE33113 series_matrix, as the original file had some empty spaces
# samples with missing values (observed with Excel): GSM1100477	GSM1100478	GSM1100479	GSM1100480	GSM1100481	GSM1100482
def specific_verify():
    print_count = 0; active = False
    with open(GEO_FOLDER + "GSE33113_series_matrix.txt", "r") as file:
        for row in file:
            if row[:26] == "!series_matrix_table_begin":
                active = True
            elif row[:8] == "\"ID_REF\"":
                num_of_samples = len(row.split("\t"))
                print(">>> Number of samples:", num_of_samples)
            elif active:
                elements = row.split("\t")
                l = len(elements)
                if l != num_of_samples:
                    print(f">>> Gene {elements[0]} has only {l} values")
                    print_count += 1
                if print_count == 10:
                    print("..."); break

# ------------------------------------------------------------------------------
# extract_row_tags()
# compare_tags()
# rough_NA_detection()
specific_verify()
