#!/usr/bin/env bash
set -euo pipefail

# Disposable, no-weight production registration check.  Config metadata is
# fetched, but no model/tokenizer weights are downloaded or allocated.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
VENV=$(mktemp -d "${TMPDIR:-/tmp}/protos-esmc-preflight.XXXXXX")
trap 'rm -rf "$VENV"' EXIT
PYTHON311=${PYTHON311:-python3.11}
command -v "$PYTHON311" >/dev/null || PYTHON311=/data/fast/envs/lambda/bin/python3.11
"$PYTHON311" -m venv --system-site-packages "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install \
  --ignore-installed \
  'transformers @ git+https://github.com/Biohub/transformers.git@ef32577f55da19a4989cd7b22e004dc43a4998cb' \
  'huggingface_hub>=0.20'
"$VENV/bin/python" -m pip install --no-deps 'accelerate>=1.0.0'
PYTHONPATH="$ROOT/src" ROOT="$ROOT" "$VENV/bin/python" - <<'PY'
import json, os
import sys
import torch
from importlib.metadata import version
from transformers import AutoConfig, AutoModelForMaskedLM
import importlib.util, sys
spec = importlib.util.spec_from_file_location("esmc_adapter", os.path.join(os.environ["ROOT"], "src/protos/processing/embedding/esmc_adapter.py"))
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)
ESMCModelProvenance, ESMCLoadPolicy, HuggingFaceESMCBackend = module.ESMCModelProvenance, module.ESMCLoadPolicy, module.HuggingFaceESMCBackend

model = ESMCModelProvenance()
assert os.path.realpath(__import__("transformers").__file__).startswith(os.path.realpath(sys.prefix))
assert module.validate_transformers_runtime() == version("transformers")
config = AutoConfig.from_pretrained(model.model_id, revision=model.model_revision, trust_remote_code=False)
assert config.model_type == "esmc" and config.d_model == 2560 and config.n_layers == 80
model_class = AutoModelForMaskedLM._model_mapping[type(config)]
assert model_class.__name__
kwargs = HuggingFaceESMCBackend.model_load_kwargs(model, ESMCLoadPolicy(), "cpu", torch_dtype="bfloat16")
assert kwargs["revision"] == model.model_revision and kwargs["trust_remote_code"] is False
from transformers.models.esmc.configuration_esmc import ESMCConfig
tiny = ESMCConfig(d_model=32, n_layers=2, vocab_size=32, n_heads=4)
tiny_model = model_class(tiny)
input_ids = torch.tensor([[0, 4, 5, 2], [0, 6, 2, 1]])
attention_mask = torch.tensor([[1, 1, 1, 1], [1, 1, 1, 0]])
with torch.no_grad():
    tiny_output = tiny_model(input_ids=input_ids, attention_mask=attention_mask, return_dict=True)
assert hasattr(tiny_output, "last_hidden_state") and tiny_output.last_hidden_state.shape[:2] == input_ids.shape
assert not hasattr(tiny_output, "embeddings")
print(json.dumps({"python": os.sys.version.split()[0], "transformers": version("transformers"), "accelerate": version("accelerate"), "transformers_module": __import__("transformers").__file__, "host_torch": version("torch"), "model_type": config.model_type, "registered_model": model_class.__name__, "tiny_output_field": "last_hidden_state", "tiny_shape": list(tiny_output.last_hidden_state.shape), "weights": "production weights not loaded", "loader_kwargs": kwargs}, sort_keys=True))
PY
