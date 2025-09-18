"""Model service architecture for containerized models.

This module provides a service-based architecture where each model runs
in its own Docker container with the appropriate dependencies.
"""

from typing import Any, Dict, List, Optional, Union
from pathlib import Path
import json
import time
import requests
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum

import numpy as np
import pandas as pd


class ServiceStatus(Enum):
    """Status of a model service."""
    OFFLINE = "offline"
    STARTING = "starting"
    READY = "ready"
    ERROR = "error"
    BUSY = "busy"


@dataclass
class ServiceConfig:
    """Configuration for a model service."""
    name: str
    image: str
    port: int
    environment: Dict[str, str]
    volumes: Dict[str, str]
    gpu: bool = False
    gpu_device: Optional[int] = None
    memory_limit: Optional[str] = None
    cpu_limit: Optional[float] = None
    healthcheck_endpoint: str = "/health"
    predict_endpoint: str = "/predict"
    
    def to_docker_config(self) -> Dict[str, Any]:
        """Convert to Docker configuration."""
        config = {
            "image": self.image,
            "ports": {f"{self.port}/tcp": self.port},
            "environment": self.environment,
            "volumes": self.volumes,
            "restart_policy": {"Name": "unless-stopped"},
            "labels": {
                "protos.model": self.name,
                "protos.type": "model-service"
            }
        }
        
        if self.gpu:
            config["device_requests"] = [{
                "driver": "nvidia",
                "count": 1 if self.gpu_device is None else -1,
                "device_ids": [str(self.gpu_device)] if self.gpu_device else [],
                "capabilities": [["gpu"]]
            }]
        
        if self.memory_limit:
            config["mem_limit"] = self.memory_limit
            
        if self.cpu_limit:
            config["cpu_quota"] = int(self.cpu_limit * 100000)
            config["cpu_period"] = 100000
            
        return config


class ModelServiceInterface(ABC):
    """Abstract interface for model services."""
    
    @abstractmethod
    def predict(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Make a prediction."""
        pass
    
    @abstractmethod
    def health_check(self) -> ServiceStatus:
        """Check service health."""
        pass
    
    @abstractmethod
    def get_info(self) -> Dict[str, Any]:
        """Get service information."""
        pass


class RemoteModelService(ModelServiceInterface):
    """Client for communicating with containerized model services."""
    
    def __init__(self, service_config: ServiceConfig, timeout: int = 30):
        self.config = service_config
        self.base_url = f"http://localhost:{self.config.port}"
        self.timeout = timeout
        
    def predict(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Make a prediction via HTTP API."""
        url = f"{self.base_url}{self.config.predict_endpoint}"
        
        # Serialize numpy arrays and pandas objects
        serialized_inputs = self._serialize_inputs(inputs)
        
        try:
            response = requests.post(
                url, 
                json=serialized_inputs,
                timeout=self.timeout
            )
            response.raise_for_status()
            
            # Deserialize response
            result = response.json()
            return self._deserialize_outputs(result)
            
        except requests.exceptions.ConnectionError:
            raise RuntimeError(f"Model service {self.config.name} is not running")
        except requests.exceptions.Timeout:
            raise RuntimeError(f"Model service {self.config.name} timed out")
        except requests.exceptions.HTTPError as e:
            raise RuntimeError(f"Model service error: {e}")
    
    def health_check(self) -> ServiceStatus:
        """Check if the service is healthy."""
        url = f"{self.base_url}{self.config.healthcheck_endpoint}"
        
        try:
            response = requests.get(url, timeout=5)
            if response.status_code == 200:
                data = response.json()
                return ServiceStatus(data.get("status", "ready"))
            else:
                return ServiceStatus.ERROR
                
        except requests.exceptions.ConnectionError:
            return ServiceStatus.OFFLINE
        except:
            return ServiceStatus.ERROR
    
    def get_info(self) -> Dict[str, Any]:
        """Get service information."""
        url = f"{self.base_url}/info"
        
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
            return response.json()
        except:
            return {"error": "Could not get service info"}
    
    def wait_for_ready(self, max_wait: int = 120, check_interval: int = 5):
        """Wait for the service to be ready."""
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.health_check()
            if status == ServiceStatus.READY:
                return True
            elif status == ServiceStatus.ERROR:
                raise RuntimeError(f"Service {self.config.name} is in error state")
            
            time.sleep(check_interval)
        
        raise TimeoutError(f"Service {self.config.name} did not become ready in {max_wait}s")
    
    def _serialize_inputs(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Serialize complex inputs for JSON transmission."""
        serialized = {}
        
        for key, value in inputs.items():
            if isinstance(value, np.ndarray):
                serialized[key] = {
                    "_type": "numpy_array",
                    "data": value.tolist(),
                    "shape": value.shape,
                    "dtype": str(value.dtype)
                }
            elif isinstance(value, pd.DataFrame):
                serialized[key] = {
                    "_type": "pandas_dataframe",
                    "data": value.to_dict(orient="records"),
                    "columns": list(value.columns),
                    "index": value.index.tolist()
                }
            elif isinstance(value, pd.Series):
                serialized[key] = {
                    "_type": "pandas_series",
                    "data": value.tolist(),
                    "index": value.index.tolist(),
                    "name": value.name
                }
            else:
                serialized[key] = value
                
        return serialized
    
    def _deserialize_outputs(self, outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Deserialize outputs from JSON."""
        deserialized = {}
        
        for key, value in outputs.items():
            if isinstance(value, dict) and "_type" in value:
                if value["_type"] == "numpy_array":
                    arr = np.array(value["data"], dtype=value.get("dtype", "float32"))
                    deserialized[key] = arr.reshape(value["shape"])
                elif value["_type"] == "pandas_dataframe":
                    deserialized[key] = pd.DataFrame(value["data"], columns=value["columns"])
                elif value["_type"] == "pandas_series":
                    deserialized[key] = pd.Series(
                        value["data"], 
                        index=value["index"], 
                        name=value.get("name")
                    )
                else:
                    deserialized[key] = value
            else:
                deserialized[key] = value
                
        return deserialized


# Model-specific service configurations
MODEL_SERVICES = {
    "esm2": ServiceConfig(
        name="esm2",
        image="protos/esm2:latest",
        port=8001,
        environment={
            "MODEL_VARIANT": "esm2_t33_650M",
            "DEVICE": "cuda"
        },
        volumes={
            "/data/models/esm2": "/models"
        },
        gpu=True,
        memory_limit="8g"
    ),
    
    "ankh": ServiceConfig(
        name="ankh",
        image="protos/ankh:latest",
        port=8002,
        environment={
            "MODEL_PATH": "/models/ankh_large",
            "DEVICE": "cuda"
        },
        volumes={
            "/data/models/ankh": "/models"
        },
        gpu=True,
        memory_limit="16g"
    ),
    
    "lambda": ServiceConfig(
        name="lambda",
        image="protos/lambda:v2",
        port=8003,
        environment={
            "MODEL_CONFIG": "/models/config.json",
            "CHECKPOINT": "/models/lambda_checkpoint.pt",
            "DEVICE": "cuda"
        },
        volumes={
            "/data/models/lambda": "/models"
        },
        gpu=True,
        memory_limit="12g"
    ),
    
    "boltz1": ServiceConfig(
        name="boltz1",
        image="protos/boltz:1.0",
        port=8004,
        environment={
            "MODEL_PATH": "/models/boltz1",
            "DEVICE": "cuda",
            "TORCH_VERSION": "2.0.1"
        },
        volumes={
            "/data/models/boltz": "/models"
        },
        gpu=True,
        memory_limit="24g"
    ),
    
    "esmfold": ServiceConfig(
        name="esmfold",
        image="protos/esmfold:latest",
        port=8005,
        environment={
            "MODEL_PATH": "/models/esmfold_v1",
            "DEVICE": "cuda",
            "OPENMM_CPU_THREADS": "4"
        },
        volumes={
            "/data/models/esmfold": "/models"
        },
        gpu=True,
        memory_limit="16g"
    )
}


class ModelServiceManager:
    """Manages Docker containers for model services."""
    
    def __init__(self, docker_host: Optional[str] = None):
        try:
            import docker
            self.docker_client = docker.from_env() if docker_host is None else docker.DockerClient(base_url=docker_host)
        except ImportError:
            raise ImportError("Docker Python SDK not installed. Run: pip install docker")
        except Exception as e:
            raise RuntimeError(f"Could not connect to Docker: {e}")
        
    def start_service(self, service_name: str, wait_ready: bool = True) -> RemoteModelService:
        """Start a model service container."""
        if service_name not in MODEL_SERVICES:
            raise ValueError(f"Unknown service: {service_name}")
            
        config = MODEL_SERVICES[service_name]
        
        # Check if already running
        container = self._get_container(service_name)
        if container and container.status == "running":
            print(f"Service {service_name} is already running")
        else:
            # Start new container
            print(f"Starting service {service_name}...")
            docker_config = config.to_docker_config()
            
            container = self.docker_client.containers.run(
                name=f"protos-{service_name}",
                detach=True,
                **docker_config
            )
            
        # Create service interface
        service = RemoteModelService(config)
        
        if wait_ready:
            print(f"Waiting for {service_name} to be ready...")
            service.wait_for_ready()
            print(f"Service {service_name} is ready!")
            
        return service
    
    def stop_service(self, service_name: str):
        """Stop a model service container."""
        container = self._get_container(service_name)
        if container:
            print(f"Stopping service {service_name}...")
            container.stop()
            container.remove()
        else:
            print(f"Service {service_name} is not running")
    
    def list_services(self) -> Dict[str, ServiceStatus]:
        """List all model services and their status."""
        status = {}
        
        for service_name in MODEL_SERVICES:
            container = self._get_container(service_name)
            if container:
                if container.status == "running":
                    # Check health
                    service = RemoteModelService(MODEL_SERVICES[service_name])
                    status[service_name] = service.health_check()
                else:
                    status[service_name] = ServiceStatus.OFFLINE
            else:
                status[service_name] = ServiceStatus.OFFLINE
                
        return status
    
    def _get_container(self, service_name: str):
        """Get container by service name."""
        try:
            return self.docker_client.containers.get(f"protos-{service_name}")
        except:
            return None