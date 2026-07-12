"""Smoke tests for the supported ProtOS package surface."""

from __future__ import annotations

import importlib.util


REMOVED_MODEL_MODULES = (
    "protos.models.model_manager",
    "protos.models.model_client",
    "protos.models.model_service",
    "protos.models.job_client",
    "protos.models.job_server",
    "protos.models.model_specs",
    "protos.models.embedding_runtime",
    "protos.models.dummy_runtime",
)

REMOVED_MODEL_EXPORTS = (
    "ModelManager",
    "ModelBackend",
    "UnifiedModelClient",
    "predict",
    "ServiceStatus",
    "ServiceConfig",
    "RemoteModelService",
    "ModelServiceManager",
    "MODEL_SERVICES",
)


def test_supported_package_imports() -> None:
    import protos
    import protos.models as models
    from protos.models import BaseModel, ModelDefinition, ModelRegistry
    from protos.processing.embedding import EmbeddingProcessor

    assert protos.get_data_path
    assert BaseModel
    assert ModelDefinition
    assert ModelRegistry
    assert EmbeddingProcessor
    assert hasattr(EmbeddingProcessor, "ingest_from_artifact")
    assert not hasattr(EmbeddingProcessor, "ingest_from_invocation")
    assert set(REMOVED_MODEL_EXPORTS).isdisjoint(models.__all__)


def test_removed_orchestration_modules_and_exports_are_absent() -> None:
    import protos.models as models

    for module_name in REMOVED_MODEL_MODULES:
        assert importlib.util.find_spec(module_name) is None

    for symbol in REMOVED_MODEL_EXPORTS:
        assert not hasattr(models, symbol)
