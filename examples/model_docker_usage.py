"""Example of using Protos models with Docker backend.

This example shows how to use models without worrying about
dependency conflicts - each model runs in its own container.
"""

from protos.models import (
    UnifiedModelClient, ModelBackend, 
    SequenceFormat, list_available_models
)


def example_docker_usage():
    """Show basic Docker-based model usage."""
    print("=== Protos Model System - Docker Backend ===\n")
    
    # Create client with Docker backend
    client = UnifiedModelClient(backend=ModelBackend.DOCKER)
    
    # List available models
    models = list_available_models()
    print(f"Available models: {models}\n")
    
    # Example 1: ESM-2 for sequence embeddings
    print("1. Using ESM-2 for embeddings:")
    try:
        result = client.predict(
            model_name="esm2",
            inputs={"sequence": "MKTAYIAKQRQISFVKSHFSRQ"},
            model_variant="esm2_t6_8M"  # Small variant for testing
        )
        
        embeddings = result.get("embeddings")
        print(f"   Embedding shape: {embeddings.shape}")
        print(f"   Per-residue embeddings generated\n")
    except Exception as e:
        print(f"   Error: {e}\n")
    
    # Example 2: Using different model with conflicting dependencies
    print("2. Using Ankh (different PyTorch version):")
    try:
        result = client.predict(
            model_name="ankh",
            inputs={"sequence": "MKTAYIAKQRQISFVKSHFSRQ"}
        )
        
        embeddings = result.get("embeddings")
        print(f"   Embedding shape: {embeddings.shape}")
        print(f"   No dependency conflicts!\n")
    except Exception as e:
        print(f"   Error: {e}\n")
    
    # Example 3: Check service status
    print("3. Service status:")
    if client.service_manager:
        services = client.service_manager.list_services()
        for name, status in services.items():
            print(f"   {name}: {status.value}")
    
    # Cleanup
    client.shutdown()


def example_auto_backend():
    """Show automatic backend selection."""
    print("\n=== Automatic Backend Selection ===\n")
    
    # Client automatically chooses best backend
    client = UnifiedModelClient(backend=ModelBackend.AUTO)
    
    # Try to use a model
    try:
        # If Docker is available, it will use that
        # Otherwise falls back to local (if dependencies installed)
        result = client.predict(
            model_name="esm2",
            inputs={"sequence": "MKTAYIAK"}
        )
        
        info = client.get_model_info("esm2")
        print(f"Model info:")
        print(f"  Name: {info['full_name']}")
        print(f"  Available backends: {info['backend_available']}")
        
    except Exception as e:
        print(f"Error: {e}")
    
    client.shutdown()


def example_entity_prediction():
    """Show prediction using entity data."""
    print("\n=== Entity-based Prediction ===\n")
    
    # This assumes you have entities loaded in Protos
    client = UnifiedModelClient(backend=ModelBackend.DOCKER)
    
    try:
        # Predict for an entity - automatically loads required formats
        result = client.predict(
            model_name="esm2",
            entity_name="BACR_HALSA"  # Bacteriorhodopsin
        )
        
        print(f"Prediction completed for BACR_HALSA")
        print(f"Results: {list(result.keys())}")
        
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure entity data is loaded in Protos")
    
    client.shutdown()


def example_quick_prediction():
    """Show quick prediction API."""
    print("\n=== Quick Prediction API ===\n")
    
    # One-liner prediction
    from protos.models import predict
    
    try:
        # Automatically manages client lifecycle
        result = predict(
            "esm2", 
            inputs={"sequence": "MKTAYIAKQRQISFVKSHFSRQ"},
            backend=ModelBackend.DOCKER
        )
        
        print("Quick prediction completed!")
        print(f"Result keys: {list(result.keys())}")
        
    except Exception as e:
        print(f"Error: {e}")


def show_docker_commands():
    """Show useful Docker commands for model management."""
    print("\n=== Docker Commands for Model Management ===\n")
    
    print("Build all model images:")
    print("  docker-compose -f docker/docker-compose.yml build")
    print()
    
    print("Start specific model:")
    print("  docker-compose -f docker/docker-compose.yml up -d esm2")
    print()
    
    print("Start all models:")
    print("  docker-compose -f docker/docker-compose.yml up -d")
    print()
    
    print("Check running models:")
    print("  docker ps --filter 'label=protos.type=model-service'")
    print()
    
    print("View model logs:")
    print("  docker logs protos-esm2")
    print()
    
    print("Stop all models:")
    print("  docker-compose -f docker/docker-compose.yml down")


def main():
    """Run all examples."""
    print("=" * 60)
    print("PROTOS MODEL SYSTEM - DOCKER BACKEND")
    print("=" * 60)
    print()
    print("This example demonstrates using AI models through Docker")
    print("containers, avoiding dependency conflicts.\n")
    
    # Check Docker availability
    try:
        import docker
        client = docker.from_env()
        client.ping()
        print("✓ Docker is available\n")
    except:
        print("✗ Docker not available. Please install Docker.\n")
        show_docker_commands()
        return
    
    # Run examples
    example_docker_usage()
    example_auto_backend()
    example_entity_prediction()
    example_quick_prediction()
    
    # Show useful commands
    show_docker_commands()
    
    print("\n" + "=" * 60)
    print("Benefits of Docker backend:")
    print("- No dependency conflicts between models")
    print("- Each model runs in isolation")
    print("- Easy to update individual models")
    print("- Consistent environment across machines")
    print("- GPU support with nvidia-docker")
    print("=" * 60)


if __name__ == "__main__":
    main()