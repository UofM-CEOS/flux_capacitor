"""Top-level main() function to be exposed as executable script

This function is run via an executable script provided by setuptools.

"""
import argparse
from underway import main as uw_main


def main():
    """Top-level tool to run underway pCO2 analyses"""
    description = ("Perform underway pCO2 calculations, " +
                   "given a configuration file.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    args = parser.parse_args()
    uw_main(args.config_file)


if __name__ == "__main__":
    main()
