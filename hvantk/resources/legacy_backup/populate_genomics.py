"""
Script to populate the genomics registry with key genomic annotation datasets.
"""
import json
import sys
from pathlib import Path
from datetime import datetime

# Add the parent directory to path
sys.path.append(str(Path(__file__).parent))
from schema_validator import SchemaValidator

def create_genomics_datasets():
    """Create genomics datasets with proper metadata."""

    validator = SchemaValidator()
    genomics_datasets = []

    # 1. dbNSFP - Missense variants prediction scores
    dbnsfp_dataset = {
        "title": "dbNSFP - Database of Human Nonsynonymous SNPs and their Functional Predictions",
        "accession": "dbNSFP_v4.7",
        "description": "Comprehensive database of functional predictions and annotations for human nonsynonymous and splice-site SNVs. Includes prediction scores from SIFT, PolyPhen2, LRT, MutationTaster, MutationAssessor, FATHMM, PROVEAN, VEST4, CADD, DANN, fitCons, PhyloP, PhastCons, GERP++, and many others.",
        "pubmedid": "21520341",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "quarterly",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 84000000,  # Approximate number of variants
        "license": "Academic use",
        "doi": "10.1002/humu.21517",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel"],
        "frequency_data": True,
        "clinical_significance": True,
        "annotation_sources": ["SIFT", "PolyPhen2", "CADD", "GERP++", "PhyloP", "PhastCons", "VEST4", "DANN", "fitCons"],
        "files": [
            {
                "path": "dbNSFP4.7a.txt.gz",
                "format": "txt",
                "size_bytes": 45000000000,  # ~45GB
                "compression": "gzip",
                "description": "Complete dbNSFP database with all prediction scores",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 2. ClinVar annotations
    clinvar_dataset = {
        "title": "ClinVar - Archive of reports of relationships among variants and phenotypes",
        "accession": "ClinVar_latest",
        "description": "Public archive of reports of the relationships among human variations and phenotypes, with supporting evidence. ClinVar facilitates access to and communication about the relationships asserted between human variation and observed health status.",
        "pubmedid": "24234437",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "monthly",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 2500000,  # Approximate number of variants
        "license": "Public Domain",
        "doi": "10.1093/nar/gkt1113",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel", "CNV", "SV"],
        "frequency_data": False,
        "clinical_significance": True,
        "annotation_sources": ["ClinVar", "OMIM", "MedGen"],
        "files": [
            {
                "path": "clinvar.vcf.gz",
                "format": "vcf",
                "size_bytes": 500000000,  # ~500MB
                "compression": "gzip",
                "description": "ClinVar variant calls with clinical significance annotations",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 3. gnomAD annotations
    gnomad_dataset = {
        "title": "gnomAD - Genome Aggregation Database",
        "accession": "gnomAD_v4.1",
        "description": "Large-scale reference dataset of human genetic variation, aggregating exome and genome sequencing data from diverse populations worldwide. Provides allele frequencies, quality metrics, and population-specific annotations.",
        "pubmedid": "32461654",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "annually",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 730000,  # ~730k individuals
        "license": "ODC Open Database License",
        "doi": "10.1038/s41586-020-2308-7",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel", "SV", "CNV"],
        "frequency_data": True,
        "clinical_significance": False,
        "population_data": {
            "populations": ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "SAS"],
            "ancestry_groups": ["African", "Latino/Admixed American", "Ashkenazi Jewish", "East Asian", "Finnish", "Non-Finnish European", "Other", "South Asian"]
        },
        "annotation_sources": ["gnomAD", "ExAC"],
        "quality_filters": {
            "min_depth": 10,
            "min_quality_score": 30,
            "genotype_quality_threshold": 20
        },
        "files": [
            {
                "path": "gnomad.genomes.v4.1.sites.vcf.bgz",
                "format": "vcf",
                "size_bytes": 150000000000,  # ~150GB
                "compression": "gzip",
                "description": "gnomAD genome sites with population frequencies",
                "checksum": "placeholder_checksum"
            },
            {
                "path": "gnomad.exomes.v4.1.sites.vcf.bgz",
                "format": "vcf",
                "size_bytes": 80000000000,  # ~80GB
                "compression": "gzip",
                "description": "gnomAD exome sites with population frequencies",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 4. INSIDER - Protein-protein interaction sites
    insider_dataset = {
        "title": "INSIDER - Interaction Network of Structured Domains in Proteins",
        "accession": "INSIDER_v1.0",
        "description": "Database of protein-protein interaction sites and interfaces derived from structural data. Provides information about interaction residues, binding sites, and functional domains involved in protein-protein interactions.",
        "pubmedid": "29036289",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "irregular",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 50000,  # Approximate number of proteins
        "license": "Academic use",
        "doi": "10.1093/nar/gkx1044",
        "genome_build": "GRCh38",
        "variant_types": ["SNV"],
        "frequency_data": False,
        "clinical_significance": False,
        "annotation_sources": ["PDB", "InterPro", "Pfam"],
        "files": [
            {
                "path": "insider_interaction_sites.tsv",
                "format": "tsv",
                "size_bytes": 100000000,  # ~100MB
                "description": "Protein-protein interaction site annotations",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 5. Ensembl gene annotations
    ensembl_dataset = {
        "title": "Ensembl Gene Annotations",
        "accession": "Ensembl_v110",
        "description": "Comprehensive gene annotations from the Ensembl project including gene models, transcript isoforms, protein sequences, regulatory features, and comparative genomics data across species.",
        "pubmedid": "31691826",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "quarterly",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 70000,  # Approximate number of genes
        "license": "Apache 2.0",
        "doi": "10.1093/nar/gkz966",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel", "SV"],
        "frequency_data": False,
        "clinical_significance": False,
        "annotation_sources": ["Ensembl", "GENCODE", "RefSeq"],
        "files": [
            {
                "path": "Homo_sapiens.GRCh38.110.gtf.gz",
                "format": "gtf",
                "size_bytes": 800000000,  # ~800MB
                "compression": "gzip",
                "description": "Complete gene annotations in GTF format",
                "checksum": "placeholder_checksum"
            },
            {
                "path": "Homo_sapiens.GRCh38.110.gff3.gz",
                "format": "gff",
                "size_bytes": 900000000,  # ~900MB
                "compression": "gzip",
                "description": "Complete gene annotations in GFF3 format",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 6. GeVIR score
    gevir_dataset = {
        "title": "GeVIR - Genome-wide Variant Intolerance Ranking",
        "accession": "GeVIR_v1.0",
        "description": "Genome-wide scores for ranking the pathogenicity potential of human genetic variants. GeVIR integrates multiple genomic features to provide pathogenicity predictions for both coding and non-coding variants.",
        "pubmedid": "31873297",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "irregular",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 700000000,  # Approximate number of scored variants
        "license": "Academic use",
        "doi": "10.1038/s41467-019-14246-7",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel"],
        "frequency_data": False,
        "clinical_significance": True,
        "annotation_sources": ["GeVIR", "CADD", "phyloP", "phastCons"],
        "files": [
            {
                "path": "gevir_scores_GRCh38.tsv.gz",
                "format": "tsv",
                "size_bytes": 20000000000,  # ~20GB
                "compression": "gzip",
                "description": "GeVIR pathogenicity scores for human variants",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    # 7. CCR - Coding-constrained regions
    ccr_dataset = {
        "title": "CCR - Coding-Constrained Regions",
        "accession": "CCR_v2.0",
        "description": "Identification of constrained coding regions in the human genome based on depletion of protein-truncating variants. CCR scores indicate the level of constraint in coding regions, useful for prioritizing variants in disease gene discovery.",
        "pubmedid": "30395236",
        "data_source": "Custom",
        "last_updated": datetime.now().isoformat(),
        "update_frequency": "irregular",
        "organism": "Homo sapiens",
        "tissue_type": "all",
        "sample_count": 20000,  # Approximate number of genes with CCR scores
        "license": "Academic use",
        "doi": "10.1038/s41467-018-07616-0",
        "genome_build": "GRCh38",
        "variant_types": ["SNV", "indel"],
        "frequency_data": False,
        "clinical_significance": True,
        "annotation_sources": ["gnomAD", "ExAC"],
        "files": [
            {
                "path": "ccr_coordinates_v2_GRCh38.bed.gz",
                "format": "bed",
                "size_bytes": 50000000,  # ~50MB
                "compression": "gzip",
                "description": "Coding-constrained regions coordinates and scores",
                "checksum": "placeholder_checksum"
            }
        ]
    }

    genomics_datasets = [
        dbnsfp_dataset,
        clinvar_dataset,
        gnomad_dataset,
        insider_dataset,
        ensembl_dataset,
        gevir_dataset,
        ccr_dataset
    ]

    # Validate all datasets
    validated_datasets = []
    for dataset in genomics_datasets:
        is_valid, errors = validator.validate_dataset(dataset, "genomics_schema")
        if is_valid:
            validated_datasets.append(dataset)
            print(f"✓ Validated: {dataset['accession']}")
        else:
            print(f"✗ Validation failed for {dataset['accession']}: {errors}")

    return validated_datasets

def update_genomics_registry():
    """Update the genomics registry with the new datasets."""

    resources_dir = Path(__file__).parent
    genomics_file = resources_dir / "registry" / "genomics" / "datasets.json"

    # Create genomics datasets
    new_datasets = create_genomics_datasets()

    # Load existing data (should be empty)
    if genomics_file.exists():
        with open(genomics_file, 'r') as f:
            existing_data = json.load(f)
    else:
        existing_data = []

    # Combine and save
    all_datasets = existing_data + new_datasets

    with open(genomics_file, 'w') as f:
        json.dump(all_datasets, f, indent=2)

    print(f"\\n✅ Successfully updated genomics registry!")
    print(f"Added {len(new_datasets)} genomics datasets:")
    for dataset in new_datasets:
        print(f"- {dataset['accession']}: {dataset['title']}")

    return len(new_datasets)

if __name__ == "__main__":
    count = update_genomics_registry()
    print(f"\\nGenomic annotation datasets now available: {count}")
