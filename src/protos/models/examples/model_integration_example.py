"""Example of integrating and using AI models in Protos.

This example demonstrates:
1. Creating and configuring models
2. Loading models from checkpoints
3. Making predictions on entities
4. Working with datasets
5. Using the model registry
"""

from pathlib import Path
from protos.models import ModelConfig, ModelRegistry, LambdaModel
from protos.io.paths import ProtosPaths


def example_lambda_model_usage():
    """Example of using Lambda model for property prediction."""
    
    # 1. Create model configuration
    config = ModelConfig(
        model_type="lambda",
        model_name="opsin_predictor_v1",
        version="1.0.0",
        input_formats=["structure", "embedding", "grn"],
        output_format="property",
        model_params={
            "hidden_dim": 256,
            "num_encoder_blocks": 4,
            "encoder_type": "GraphTransformer",
            "pooling_type": "multihead_attention",
            "num_heads": 8,
            "targets": ["lambda_max", "photocycle_rate"],
            "tasks": ["regression", "regression"]
        },
        device="cuda" if torch.cuda.is_available() else "cpu"
    )
    
    # 2. Create model instance
    model = LambdaModel(config=config)
    
    # 3. Load pre-trained weights
    checkpoint_path = Path("models/lambda/opsin_predictor_v1/weights/best_model.pth")
    if checkpoint_path.exists():
        model.load_model(checkpoint_path)
    
    # 4. Make predictions on single entity
    prediction = model.predict("BACR_HALSA")
    print(f"Prediction for BACR_HALSA: {prediction.prediction}")
    
    # 5. Create and predict on dataset
    model.create_dataset(
        "microbial_opsins",
        ["BACR_HALSA", "ChR2", "NpHR", "C1C2"],
        metadata={"organism": "microbial", "family": "rhodopsin"}
    )
    
    dataset_predictions = model.predict_dataset("microbial_opsins")
    for pred in dataset_predictions:
        print(f"{pred.entity_name}: {pred.prediction}")


def example_model_registry_usage():
    """Example of using the model registry."""
    
    # 1. Initialize registry
    registry = ModelRegistry()
    
    # 2. Discover available models
    models = registry.discover_models()
    print(f"Found {len(models)} models")
    
    # 3. List models of specific type
    lambda_models = registry.list_models(model_type="lambda")
    print(f"Lambda models: {lambda_models}")
    
    # 4. Load a model from registry
    if lambda_models:
        model = registry.load_model(lambda_models[0])
        print(f"Loaded model: {model.name}")
    
    # 5. Find models compatible with an entity
    compatible_models = registry.get_models_for_entity(
        "BACR_HALSA",
        output_format="property"
    )
    print(f"Models that can process BACR_HALSA: {compatible_models}")


def example_custom_model_integration():
    """Example of integrating a custom model."""
    
    from protos.models import BaseModel, ModelConfig
    import numpy as np
    
    class SimplePropertyPredictor(BaseModel):
        """A simple custom property predictor."""
        
        def load_model(self, checkpoint_path=None):
            """Load a simple linear model."""
            # In practice, load your model here
            self.model = lambda x: np.random.rand()  # Dummy predictor
            self.is_loaded = True
        
        def _preprocess_input(self, input_data):
            """Convert input data to model format."""
            # Extract features from input_data
            if 'sequence' in input_data:
                # Simple: use sequence length as feature
                return len(input_data['sequence'])
            return 0
        
        def _predict_single(self, input_data):
            """Make prediction."""
            # Use preprocessed input
            return {"property_value": self.model(input_data)}
    
    # Register the custom model type
    from protos.models import register_model_type
    register_model_type("simple_predictor", SimplePropertyPredictor)
    
    # Create and use the model
    config = ModelConfig(
        model_type="simple_predictor",
        model_name="test_model",
        version="1.0",
        input_formats=["sequence"],
        output_format="property"
    )
    
    model = SimplePropertyPredictor(config=config)
    model.load_model()
    
    # Make prediction
    result = model.predict("P62988")  # UniProt ID
    print(f"Custom model prediction: {result.prediction}")


def example_multi_format_prediction():
    """Example using multiple input formats."""
    
    # Configure model that uses multiple inputs
    config = ModelConfig(
        model_type="lambda",
        model_name="multimodal_predictor",
        version="2.0",
        input_formats=["structure", "sequence", "grn", "property"],
        output_format="property",
        model_params={
            "fusion_method": "attention",
            "modality_dropout": 0.1
        }
    )
    
    # The model will automatically:
    # 1. Load structure from CifBaseProcessor
    # 2. Load sequence from SeqProcessor  
    # 3. Load GRN from GRNBaseProcessor
    # 4. Load existing properties from PropertyProcessor
    # 5. Combine all inputs for prediction
    
    model = LambdaModel(config=config)
    # model.load_model()  # Load pre-trained weights
    
    # Prediction will use all available formats for the entity
    # prediction = model.predict("complex_protein")


def example_model_adapter_pattern():
    """Example of using adapters for format conversion."""
    
    from protos.models import LambdaModel, ModelConfig
    
    # Define custom adapters
    def sequence_to_graph_adapter(sequence_data):
        """Convert sequence to graph representation."""
        # Custom logic to create graph from sequence
        # This is a simplified example
        import torch
        from torch_geometric.data import Data
        
        seq_len = len(sequence_data) if isinstance(sequence_data, str) else sequence_data.shape[0]
        
        # Create fully connected graph
        edge_index = torch.combinations(torch.arange(seq_len), r=2).t()
        
        # Use sequence embeddings or one-hot encoding as features
        if isinstance(sequence_data, str):
            # One-hot encode sequence
            features = torch.randn(seq_len, 20)  # Simplified
        else:
            features = torch.tensor(sequence_data, dtype=torch.float32)
        
        return Data(x=features, edge_index=edge_index)
    
    # Create model with adapter
    config = ModelConfig(
        model_type="lambda",
        model_name="sequence_graph_model",
        version="1.0",
        input_formats=["sequence"],
        output_format="property"
    )
    
    model = LambdaModel(config=config)
    
    # Register the adapter
    model.register_input_adapter("sequence", sequence_to_graph_adapter)
    
    # Now the model can process sequences by converting them to graphs
    # model.load_model()
    # prediction = model.predict("EGFR_HUMAN")


if __name__ == "__main__":
    import torch  # Import for device check
    
    print("=== Protos Model Integration Examples ===\n")
    
    print("1. Lambda Model Usage")
    print("-" * 40)
    # example_lambda_model_usage()
    
    print("\n2. Model Registry Usage")
    print("-" * 40)
    example_model_registry_usage()
    
    print("\n3. Custom Model Integration")
    print("-" * 40)
    example_custom_model_integration()
    
    print("\n4. Adapter Pattern")
    print("-" * 40)
    example_model_adapter_pattern()
    
    print("\nExamples completed!")