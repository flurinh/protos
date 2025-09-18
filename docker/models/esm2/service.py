"""ESM-2 model service."""

import os
import torch
import esm_test
from flask import Flask, request, jsonify
import numpy as np
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create Flask app
app = Flask(__name__)

# Global model holder
model = None
alphabet = None
batch_converter = None
device = None


def load_model():
    """Load ESM-2 model."""
    global model, alphabet, batch_converter, device
    
    model_variant = os.environ.get("MODEL_VARIANT", "esm2_t33_650M")
    device_name = os.environ.get("DEVICE", "cuda" if torch.cuda.is_available() else "cpu")
    
    logger.info(f"Loading ESM-2 model: {model_variant} on {device_name}")
    
    # Load model
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_variant)
    batch_converter = alphabet.get_batch_converter()
    
    # Move to device
    device = torch.device(device_name)
    model = model.to(device)
    model.eval()
    
    logger.info("Model loaded successfully")


@app.route("/health", methods=["GET"])
def health():
    """Health check endpoint."""
    if model is None:
        return jsonify({"status": "starting"}), 503
    return jsonify({"status": "ready"}), 200


@app.route("/info", methods=["GET"])
def info():
    """Model information endpoint."""
    return jsonify({
        "model": "esm2",
        "variant": os.environ.get("MODEL_VARIANT", "esm2_t33_650M"),
        "device": str(device) if device else "not loaded",
        "max_length": 1024,
        "embedding_dim": 1280
    })


@app.route("/predict", methods=["POST"])
def predict():
    """Prediction endpoint."""
    if model is None:
        return jsonify({"error": "Model not loaded"}), 503
    
    try:
        data = request.json
        
        # Extract sequences
        if "sequence" not in data:
            return jsonify({"error": "No sequence provided"}), 400
        
        sequences = data["sequence"]
        if isinstance(sequences, str):
            sequences = [sequences]
        
        # Prepare batch
        batch_data = [(f"protein_{i}", seq) for i, seq in enumerate(sequences)]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)
        
        # Get representations
        with torch.no_grad():
            results = model(
                batch_tokens,
                repr_layers=[-1],
                return_contacts=data.get("return_contacts", False)
            )
        
        # Extract embeddings
        embeddings = results["representations"][-1]
        
        # Remove special tokens
        embeddings = embeddings[:, 1:-1, :]
        
        # Prepare response
        response = {
            "embeddings": {
                "_type": "numpy_array",
                "data": embeddings.cpu().numpy().tolist(),
                "shape": list(embeddings.shape),
                "dtype": "float32"
            }
        }
        
        if "contacts" in results:
            contacts = results["contacts"].cpu().numpy()
            response["contacts"] = {
                "_type": "numpy_array",
                "data": contacts.tolist(),
                "shape": list(contacts.shape),
                "dtype": "float32"
            }
        
        return jsonify(response)
        
    except Exception as e:
        logger.error(f"Prediction error: {str(e)}")
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    # Load model on startup
    load_model()
    
    # Run server
    app.run(host="0.0.0.0", port=8001)