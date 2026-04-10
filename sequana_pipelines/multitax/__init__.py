from importlib.metadata import PackageNotFoundError, version

try:
    version = version("sequana-multitax")
except PackageNotFoundError:
    version = "unknown"
