import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.commands.ucsc_downloader import ucsc_downloader
from hvantk.commands.expression_atlas_downloader import download_experiments
from hvantk.commands.clingen_downloader import clingen_downloader
from hvantk.commands.gencc_downloader import gencc_downloader
from hvantk.commands.hgnc_downloader import hgnc_downloader
from hvantk.commands.clinvar_downloader import clinvar_downloader
from hvantk.commands.uniprot_ptm_downloader import uniprot_ptm_downloader
from hvantk.commands.peptideatlas_phospho_downloader import peptideatlas_phospho_downloader
from hvantk.commands.cptac_phospho_downloader import cptac_phospho_downloader


@click.group("download", context_settings=CONTEXT_SETTINGS)
def download_group():
    """Download external datasets."""


download_group.add_command(ucsc_downloader, "ucsc")
download_group.add_command(download_experiments, "expression-atlas")
download_group.add_command(clingen_downloader, "clingen")
download_group.add_command(gencc_downloader, "gencc")
download_group.add_command(hgnc_downloader, "hgnc")
download_group.add_command(clinvar_downloader, "clinvar")
download_group.add_command(uniprot_ptm_downloader, "uniprot-ptm")
download_group.add_command(peptideatlas_phospho_downloader, "peptideatlas-phospho")
download_group.add_command(cptac_phospho_downloader, "cptac-phospho")
