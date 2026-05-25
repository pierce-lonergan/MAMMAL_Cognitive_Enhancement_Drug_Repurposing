"""Local equivalent of notebooks/00_smoke_test.ipynb — single inference call.

Used to validate env + GPU + weight download before kicking off the full pipeline.
"""

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mammal_repurposing  # noqa: F401 -- installs sys.modules shims

import torch
from fuse.data.tokenizers.modular_tokenizer.op import ModularTokenizerOp
from mammal.examples.dti_bindingdb_kd.task import DtiBindingdbKdTask
from mammal.model import Mammal

MODEL_ID = "ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd"
NORM_Y_MEAN = 5.79384684128215
NORM_Y_STD = 1.33808027428196

# Reference inputs from HF model card README
target_seq = "NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"
drug_seq = "CC(=O)NCCC1=CNc2c1cc(OC)cc2"

print(f"PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Device: {torch.cuda.get_device_name(0)}")

t0 = time.perf_counter()
print(f"\nLoading model {MODEL_ID} ...")
model = Mammal.from_pretrained(MODEL_ID)
model.eval()
if torch.cuda.is_available():
    model.to("cuda")
tokenizer_op = ModularTokenizerOp.from_pretrained(MODEL_ID)
load_secs = time.perf_counter() - t0
n_params = sum(p.numel() for p in model.parameters())
print(f"Loaded in {load_secs:.1f}s. Params: {n_params/1e6:.1f}M. Device: {model.device}")

print("\nRunning reference inference (HF README pair)...")
t1 = time.perf_counter()
sample = {"target_seq": target_seq, "drug_seq": drug_seq}
sample = DtiBindingdbKdTask.data_preprocessing(
    sample_dict=sample,
    tokenizer_op=tokenizer_op,
    target_sequence_key="target_seq",
    drug_sequence_key="drug_seq",
    norm_y_mean=None,
    norm_y_std=None,
    device=model.device,
)
batch_dict = model.forward_encoder_only([sample])
batch_dict = DtiBindingdbKdTask.process_model_output(
    batch_dict,
    scalars_preds_processed_key="model.out.dti_bindingdb_kd",
    norm_y_mean=NORM_Y_MEAN,
    norm_y_std=NORM_Y_STD,
)
pkd_ref = float(batch_dict["model.out.dti_bindingdb_kd"][0])
infer_secs = time.perf_counter() - t1
print(f"Reference pair pKd = {pkd_ref:.3f}  ({infer_secs*1000:.0f} ms)")

if torch.cuda.is_available():
    vram_mb = torch.cuda.memory_allocated() / 1e6
    peak_mb = torch.cuda.max_memory_allocated() / 1e6
    print(f"\nVRAM allocated: {vram_mb:.0f} MB (peak {peak_mb:.0f} MB)")

print("\nSMOKE TEST OK")
