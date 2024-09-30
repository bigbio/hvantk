import click

from hvantk.settings import CONTEXT_SETTINGS
from hvantk.commands.make_annotation_tables_cli import make_annotation_tables_cli


# Main CLI entry point for the package (hvantk)


@click.group('hvantk',
             help='A python package for gene and variant annotation.',
             context_settings=CONTEXT_SETTINGS)
def cli():
    """A python package for gene and variant annotation."""
    pass


cli.add_command(make_annotation_tables_cli)

def main():
    cli()


if __name__ == '__main__':
    main()
