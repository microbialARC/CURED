#!/usr/bin/env python3
import argparse
from cured import CURED_Main
from cured import CURED_FindREs

def main():
    parser = argparse.ArgumentParser(prog="CURED")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # main subcommand
    parser_main = subparsers.add_parser("main", help="Run main workflow")
    CURED_Main.add_arguments(parser_main)
    parser_main.set_defaults(func=CURED_Main.main)

    # res subcommand
    parser_res = subparsers.add_parser("res", help="Find restriction enzymes")
    CURED_FindREs.add_arguments(parser_res)
    parser_res.set_defaults(func=CURED_FindREs.main)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()