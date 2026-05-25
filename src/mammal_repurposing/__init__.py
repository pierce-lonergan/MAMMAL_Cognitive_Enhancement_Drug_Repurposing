"""MAMMAL-based drug repurposing pipeline for cognitive enhancement targets."""

# Install sys.modules shims BEFORE any submodule imports mammal/fuse.
# (PyTDC eagerly walks cellxgene/tiledbsoma which has no Windows wheel.)
from mammal_repurposing._compat import install_tdc_shims as _install_tdc_shims

_install_tdc_shims()

__version__ = "0.1.0"
