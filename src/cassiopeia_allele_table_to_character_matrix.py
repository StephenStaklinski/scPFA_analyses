
import sys 
import pandas as pd
import cassiopeia as cas


alleletable_path = sys.argv[1]
charmatrix_path = sys.argv[2]
mutdict_path = sys.argv[3]

# Read in the allele table
allele_table = pd.read_csv(alleletable_path, sep='\t')

# Convert to character matrix
char_matrix, prob_dict, dict_mapping = cas.pp.convert_alleletable_to_character_matrix(allele_table)

# Print raw data stats
initial_ncells, initial_nsites = char_matrix.shape
print(f"Starting with {initial_ncells} cells and {initial_nsites} sites")

# Remove any columns with all (NaN || [None]) since those are uninformative sites
cols_to_drop = char_matrix.columns[char_matrix.isin([0, -1]).all()]
char_matrix = char_matrix.drop(columns=cols_to_drop)
filtered_ncells, filtered_nsites = char_matrix.shape
print(f"After removing uninformative sites, {filtered_ncells} cells and {filtered_nsites} sites remain")

# Remove dropped columns from the mutation dict
cols_to_drop_renamed = [int(name.replace("r", "")) - 1 for name in cols_to_drop]    # Needs to be offset by one to match
dict_mapping = {k: v for k, v in dict_mapping.items() if k not in cols_to_drop_renamed}

# Rename cols to be incremental
char_matrix.columns = [f"r{idx+1}" for idx in range(char_matrix.shape[1])]
dict_mapping = {f"r{idx+1}": indel_dict for idx, (site, indel_dict) in enumerate(dict_mapping.items())}

# Remove any rows with all (NaN || [None]) since those are uninformative cells
rows_to_drop = char_matrix.index[char_matrix.isin([0, -1]).all(axis=1)]
char_matrix = char_matrix.drop(index=rows_to_drop)
filtered_ncells, filtered_nsites = char_matrix.shape
print(f"After removing uninformative cells, {filtered_ncells} cells and {filtered_nsites} sites remain")

# Write the character matrix to a tsv file
char_matrix.to_csv(charmatrix_path, sep='\t')

# Write out the mutation dictionary to a tsv file
with open(mutdict_path, 'w') as f:
    f.write("site\tindel_int\tindel_actual\n")
    for site, indel_dict in dict_mapping.items():
        for indel_int, indel_actual in indel_dict.items():
            f.write(f"{site}\t{indel_int}\t{indel_actual}\n")

