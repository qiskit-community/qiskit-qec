"""Factory for error propagators."""
import logging


class EPSelector:
    """Select and return an ErrorPropagator."""

    def __init__(self):
        """Create a new factory."""
        self.logger = logging.getLogger(__name__)

    def _attempt_import(self, module_name: str):
        """Try to load a module."""
        try:
            __import__(module_name)
        except ImportError as e:
            self.logger.exception(f"__import__({module_name}) failed, raising {e}")
            return False
        else:
            return True

    def get_error_propagator(self, eptype: str = "auto", qreg_size: int = 1, creg_size: int = 1):
        """Return an error propagator.

        auto = try to automatically detect extensions
        c = force the compiled version
        py = force the python version
        """
        if eptype == "auto":
            if self._attempt_import("qiskit_qec.extensions.compiledextension"):
                from qiskit_qec.analysis.cerrorpropagator import CErrorPropagator

                return CErrorPropagator(qreg_size, creg_size)
            else:
                from qiskit_qec.analysis.pyerrorpropagator import PyErrorPropagator

                return PyErrorPropagator(qreg_size, creg_size)
        elif eptype == "c":
            if not self._attempt_import("qiskit_qec.extensions.compiledextension"):
                raise Exception("extensions not available")
            from qiskit_qec.analysis.cerrorpropagator import CErrorPropagator

            return CErrorPropagator(qreg_size, creg_size)
        elif eptype == "py":
            from qiskit_qec.analysis.pyerrorpropagator import PyErrorPropagator

            return PyErrorPropagator(qreg_size, creg_size)
        else:
            raise Exception("unsupported error propagator type")
