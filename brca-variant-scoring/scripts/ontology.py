# scripts/ontology.py

from alphagenome.models import dna_client

# Pass your key directly for local use
api_key = "AIzaSyBpZVKZljO2x2HesJZH9cCPrR9qwVOhEN4"
dna_model = dna_client.create(api_key)

# Get all metadata
metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

# Search for breast-related tracks
breast = metadata[
    metadata['biosample_name'].str.contains('breast|mammary|MCF', case=False, na=False)
]

print("=== Breast-related tracks ===")
print(breast[['output_type', 'biosample_name', 'ontology_curie']].drop_duplicates().to_string())
print(f"\nUnique ontology terms: {breast['ontology_curie'].unique()}")