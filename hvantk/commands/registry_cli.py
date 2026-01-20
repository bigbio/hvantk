#!/usr/bin/env python3
"""
CLI command for the new registry system.

This integrates with the existing CLI structure instead of being a standalone script.
"""

import argparse
from pathlib import Path

from hvantk.register import RegistryManager, RegistryConfig


def registry_web_command(args):
    """Handle the registry web generation command."""
    config = RegistryConfig()
    if args.title:
        config.site_title = args.title
    if args.theme:
        config.theme_color = args.theme

    manager = RegistryManager(registry_file=args.input, config=config)
    output_path = manager.generate_web_interface(args.output)

    stats = manager.get_registry_statistics()
    print(f"✅ Web registry generated in {output_path}")
    print(f"📊 Total datasets: {stats['total_datasets']}")
    print(f"📈 Success rate: {stats['success_rate']}%")


def registry_api_command(args):
    """Handle the registry API generation command."""
    config = RegistryConfig()
    manager = RegistryManager(registry_file=args.input, config=config)
    api_path = manager.generate_api_endpoints(args.output_dir)

    print(f"✅ API endpoints generated in {api_path}")


def registry_complete_command(args):
    """Handle the complete registry generation command."""
    config = RegistryConfig()
    if args.title:
        config.site_title = args.title
    if args.theme:
        config.theme_color = args.theme

    manager = RegistryManager(registry_file=args.input, config=config)
    paths = manager.generate_complete_registry(args.output)

    stats = manager.get_registry_statistics()
    print(f"✅ Complete registry generated in {args.output}")
    print(f"🌐 Web interface: {paths['web_interface']}")
    print(f"📡 API endpoints: {paths['api_endpoints']}")
    print(
        f"📊 Total: {stats['total_datasets']} datasets ({stats['success_rate']}% success)"
    )


def main():
    """Main CLI entry point for registry commands."""
    parser = argparse.ArgumentParser(description="HVANTK Registry Management")
    subparsers = parser.add_subparsers(dest="command", help="Registry commands")

    # Web generation command
    web_parser = subparsers.add_parser("web", help="Generate web interface")
    web_parser.add_argument("--input", required=True, help="Input registry JSON file")
    web_parser.add_argument("--output", required=True, help="Output directory")
    web_parser.add_argument("--title", help="Custom site title")
    web_parser.add_argument("--theme", help="Custom theme color")
    web_parser.set_defaults(func=registry_web_command)

    # API generation command
    api_parser = subparsers.add_parser("api", help="Generate API endpoints")
    api_parser.add_argument("--input", required=True, help="Input registry JSON file")
    api_parser.add_argument("--output-dir", required=True, help="Output directory")
    api_parser.set_defaults(func=registry_api_command)

    # Complete generation command
    complete_parser = subparsers.add_parser(
        "complete", help="Generate complete registry"
    )
    complete_parser.add_argument(
        "--input", required=True, help="Input registry JSON file"
    )
    complete_parser.add_argument("--output", required=True, help="Output directory")
    complete_parser.add_argument("--title", help="Custom site title")
    complete_parser.add_argument("--theme", help="Custom theme color")
    complete_parser.set_defaults(func=registry_complete_command)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
