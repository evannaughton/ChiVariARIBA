import pandas as pd

# Load the ARIBA output file
file_path = 'report.tsv'
ariba_output = pd.read_csv(file_path, sep='\t')

# Define the variant-specific columns (existing logic)
variant_columns = ariba_output.columns[14:30]

# Define core columns that should remain identical for each unique gene entry, excluding 'free_text'
core_columns = [col for col in ariba_output.columns if col not in variant_columns and col != 'free_text']

# Define custom names for the 18 columns resulting from the split of 'free_text'
free_text_column_names = [
    'gff_file_of_gene', 'scaffold_name', 'annotation/protein_id', 'gene_name', 'gene_description/annotation',
    'panaroo_gene_family', 'reference_gene_id', 'pident_compared_to_reference', 'length', 'mismatches',
    'gapopen', 'qstart', 'qend', 'sstart', 'send',
    'evalue', 'bitscore', 'btop_mismatches'
]

# Split the 'free_text' column based on '|' and expand it into multiple columns
if 'free_text' in ariba_output.columns:
    free_text_split = ariba_output['free_text'].str.split('|', expand=True)
    
    # Ensure we have exactly 18 columns after splitting
    if free_text_split.shape[1] == 18:
        free_text_split.columns = free_text_column_names  # Assign custom column names
    else:
        raise ValueError(f"Expected 18 columns from free_text split, but got {free_text_split.shape[1]}")

    # Combine the split columns with the original DataFrame
    ariba_output = pd.concat([ariba_output.drop(columns=['free_text']), free_text_split], axis=1)

# Re-define core columns to include the newly named split 'free_text' columns
core_columns = [col for col in ariba_output.columns if col not in variant_columns]

# Group by the core columns and aggregate the variant columns by joining unique values into a single string
aggregated_output = (
    ariba_output.groupby(core_columns, as_index=False)
    .agg({col: lambda x: "; ".join(x.dropna().unique()) for col in variant_columns})
)

# Save the output to a new file with the collapsed rows and split free_text columns at the end
output_path = 'ariba_collapsed_report.tsv'
aggregated_output.to_csv(output_path, sep='\t', index=False)

print(f"Collapsed file with named free_text columns saved to: {output_path}")








