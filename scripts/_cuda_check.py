"""Quick CUDA + GPU verification. Used by setup script."""
import sys

import torch

print(f"PyTorch         : {torch.__version__}")
print(f"CUDA available  : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Device          : {torch.cuda.get_device_name(0)}")
    print(f"Compute cap.    : {torch.cuda.get_device_capability(0)}")
    print(f"CUDA runtime    : {torch.version.cuda}")
else:
    print("WARNING: CUDA not available - inference will run on CPU (slow).", file=sys.stderr)
    sys.exit(2)
