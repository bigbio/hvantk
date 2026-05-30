"""Plugin-local constants for the uniprot_ptm skill.

Moved from hvantk/core/ptm_constants.py per issue #122 (clean-core principle).
"""

# UniProt REST API
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_API_FIELDS = "accession,gene_names,ft_mod_res,xref_ensembl,sequence"
UNIPROT_HUMAN_PTM_QUERY = "organism_id:9606 AND reviewed:true AND ft_mod_res:*"
UNIPROT_BATCH_SIZE = 500
