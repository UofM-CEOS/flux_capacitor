"""
Basic flags and objects required throughout the package

"""

__all__ = ["FLUX_FLAGS"]

#: Flags used during processing
FLUX_FLAGS = ["failed_prep_flag", "open_flag", "closed_flag",
              "sonic_flag", "motion_flag", "bad_navigation_flag",
              "bad_meteorology_flag"]
