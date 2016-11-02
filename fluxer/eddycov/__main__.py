"""Top-level main() function to be exposed as an executable script

This function is run via an executable script via setuptools.

"""


import argparse
import fluxer.eddycov.db_flux as db_flux


def main():
    """Top-level tool to run full flux analyses"""
    description = "Perform flux analyses, given a configuration file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    args = parser.parse_args()
    db_flux.main(args.config_file)


if __name__ == "__main__":
    main()
