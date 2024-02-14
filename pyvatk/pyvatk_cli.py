import click

from pyvatk.settings import CONTEXT_SETTINGS
from pyvatk.cli.make_annotation_tables_cli import make_annotation_tables_cli


# Main CLI entry point for the package (pyvatk)


@click.group('pyvatk',
             help='A python package for gene and variant annotation.',
             context_settings=CONTEXT_SETTINGS)
def cli():
    """A python package for gene and variant annotation."""
    pass


cli.add_command(make_annotation_tables_cli)


def pyvatk_main():
    cli()


if __name__ == '__main__':
    pyvatk_main()
