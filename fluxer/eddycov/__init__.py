FLUX_FLAGS = ["open_flag", "closed_flag", "sonic_flag",
              "motion_flag", "bad_navigation_flag",
              "bad_meteorology_flag"]


# Exception classes for catching our conditions
class FluxError(Exception):
    """Base class for Exceptions in this module"""
    pass


class SonicError(FluxError):
    """Critical sonic anemometer Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags


class NavigationError(FluxError):
    """Critical navigation Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags


class MeteorologyError(FluxError):
    """Critical meteorology Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags
