"""Test Docker backend functionality."""

from protos.models import (
    UnifiedModelClient, ModelBackend, 
    ServiceStatus, MODEL_SERVICES,
    list_available_models
)


def test_docker_backend():
    """Test basic Docker backend functionality."""
    print("=== Testing Docker Backend ===\n")
    
    # Check Docker availability
    try:
        import docker
        docker_client = docker.from_env()
        docker_client.ping()
        print("✓ Docker is available")
    except Exception as e:
        print(f"✗ Docker not available: {e}")
        print("\nTo use Docker backend:")
        print("1. Install Docker: https://docs.docker.com/get-docker/")
        print("2. Start Docker daemon")
        return
    
    # Create client
    print("\nCreating Docker-based client...")
    client = UnifiedModelClient(backend=ModelBackend.DOCKER)
    
    # List available models
    models = list_available_models()
    print(f"\nAvailable models: {models}")
    
    # Show service configurations
    print("\nModel service configurations:")
    for name, config in MODEL_SERVICES.items():
        print(f"\n{name}:")
        print(f"  Image: {config.image}")
        print(f"  Port: {config.port}")
        print(f"  GPU: {config.gpu}")
        print(f"  Memory: {config.memory_limit}")
    
    # Check service status
    if client.service_manager:
        print("\nChecking service status...")
        services = client.service_manager.list_services()
        for name, status in services.items():
            print(f"  {name}: {status.value}")
    
    # Cleanup
    client.shutdown()


def test_service_config():
    """Test service configuration functionality."""
    print("\n\n=== Testing Service Configuration ===\n")
    
    from protos.models import ServiceConfig
    
    # Example custom service config
    custom_config = ServiceConfig(
        name="custom_model",
        image="myorg/custom-model:latest",
        port=8100,
        environment={
            "MODEL_PATH": "/models/custom",
            "DEVICE": "cuda",
            "BATCH_SIZE": "32"
        },
        volumes={
            "/data/models/custom": "/models"
        },
        gpu=True,
        memory_limit="16g",
        healthcheck_endpoint="/api/health",
        predict_endpoint="/api/predict"
    )
    
    print("Custom service configuration:")
    print(f"  Name: {custom_config.name}")
    print(f"  Port: {custom_config.port}")
    
    # Convert to Docker config
    docker_config = custom_config.to_docker_config()
    print("\nDocker configuration:")
    for key, value in docker_config.items():
        if key != "environment":  # Skip verbose env vars
            print(f"  {key}: {value}")


def test_model_info():
    """Test model information retrieval."""
    print("\n\n=== Testing Model Information ===\n")
    
    client = UnifiedModelClient()
    
    # Get info for each model
    for model_name in ["esm2", "ankh", "lambda"]:
        info = client.get_model_info(model_name)
        print(f"\n{info['full_name']}:")
        print(f"  Description: {info['description']}")
        print(f"  Inputs: {info['input_formats']}")
        print(f"  Output: {info['output_format']}")
        print(f"  Backends: {info['backend_available']}")


def main():
    """Run all tests."""
    print("=" * 60)
    print("PROTOS DOCKER BACKEND TEST")
    print("=" * 60)
    
    test_docker_backend()
    test_service_config()
    test_model_info()
    
    print("\n" + "=" * 60)
    print("Docker backend provides:")
    print("- Isolated environments for each model")
    print("- No dependency conflicts")
    print("- Easy scaling and deployment")
    print("- Consistent behavior across machines")
    print("=" * 60)


if __name__ == "__main__":
    main()