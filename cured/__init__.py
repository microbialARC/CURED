#!/usr/bin/env python3
import argparse
from . import CURED_Main
from . import CURED_FindREs

def main():
    """Main entrypoint for the CURED console command."""
    parser = argparse.ArgumentParser(prog="cured", description="CURED: Classification Using Restriction Enzyme Diagnostics")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # main subcommand
    parser_main = subparsers.add_parser("main", help="Run main workflow")
    CURED_Main.add_arguments(parser_main)  # Reuse argument definitions
    parser_main.set_defaults(func=CURED_Main.main)

    # res subcommand
    parser_res = subparsers.add_parser("res", help="Find restriction enzymes")
    CURED_FindREs.add_arguments(parser_res)
    parser_res.set_defaults(func=CURED_FindREs.main)

    # Parse arguments and dispatch
    args = parser.parse_args()
    args.func(args)
