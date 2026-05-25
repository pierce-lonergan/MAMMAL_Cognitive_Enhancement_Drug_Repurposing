"""Runtime compatibility shims.

The MAMMAL example task code (``mammal.examples.dti_bindingdb_kd.pl_data_module``)
unconditionally imports ``from tdc.multi_pred.dti import DTI``. PyTDC's
``tdc.multi_pred.__init__.py`` eagerly imports every sibling module, including
``perturboutcome`` which pulls in ``single_cell`` -> ``cellxgene_census`` ->
``gget`` -> ``tiledbsoma``. ``tiledbsoma`` has no buildable Windows wheel.

We don't use cellxgene or single-cell features for DTI / BBBP / ClinTox inference.
This shim stubs the problematic transitive modules in ``sys.modules`` BEFORE
the mammal import chain executes, so the eager imports succeed without ever
loading the native cellxgene code.

Stubbed module APIs raise ``RuntimeError`` on access — if you accidentally call
into the cell-genes path, you'll see a loud failure rather than silent breakage.
"""

from __future__ import annotations

import sys
import types


def _make_stub(module_name: str, attrs: dict | None = None) -> types.ModuleType:
    mod = types.ModuleType(module_name)
    # Set common dunder defaults so inspect / pickle / debuggers don't trip the
    # raising __getattr__ below.
    mod.__file__ = f"<stub:{module_name}>"
    mod.__path__ = []  # marks it as a package
    mod.__all__ = []
    attrs = attrs or {}
    for k, v in attrs.items():
        setattr(mod, k, v)

    def _raising_getattr(name: str):
        # Let Python's machinery query dunders without raising — those are
        # used by inspect, pickle, debuggers, IDEs.
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        raise RuntimeError(
            f"{module_name}.{name} called via stub. "
            "Real package not installed (tiledbsoma wheel fails on Windows). "
            "This path is not exercised by DTI/BBBP/ClinTox inference; if you "
            "hit it, the call chain has drifted."
        )

    mod.__getattr__ = _raising_getattr  # PEP 562 module-level __getattr__
    return mod


def install_tdc_shims() -> None:
    """Install sys.modules stubs for tdc/cellxgene cascade. Idempotent."""
    stubs = {
        "tiledbsoma": {"__version__": "0.0.0-stub"},
        "cellxgene_census": {"__version__": "0.0.0-stub"},
        "gget": {"__version__": "0.0.0-stub"},
        "tdc.multi_pred.perturboutcome": {"PerturbOutcome": type(
            "PerturbOutcome", (), {"__init__": lambda self, *a, **kw: None}
        )},
        "tdc.multi_pred.single_cell": {"CellXGeneTemplate": type(
            "CellXGeneTemplate", (), {"__init__": lambda self, *a, **kw: None}
        )},
        "tdc.resource.cellxgene_census": {"CensusResource": type(
            "CensusResource", (), {"__init__": lambda self, *a, **kw: None}
        )},
    }
    for name, attrs in stubs.items():
        if name not in sys.modules:
            sys.modules[name] = _make_stub(name, attrs)
