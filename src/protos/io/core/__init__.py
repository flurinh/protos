"""Core data management infrastructure."""

from .base_processor import BaseProcessor
from .entity_registry import EntityRegistry
from .dataset_manager import DatasetManager
from .input_manager import InputManager
from .registry_health import RegistryHealthChecker
from .conflict_resolver import ConflictResolver


# Singleton access for EntityRegistry
_registry_instance = None


def get_registry() -> EntityRegistry:
    """Get or create EntityRegistry singleton."""
    global _registry_instance
    if _registry_instance is None:
        _registry_instance = EntityRegistry()
    return _registry_instance


def reset_registry(*, backup: bool = True):
    """Reset the singleton registry contents and optionally create a backup."""

    registry = get_registry()
    registry.reset(backup=backup)


__all__ = [
    'BaseProcessor',
    'EntityRegistry',
    'DatasetManager',
    'InputManager',
    'RegistryHealthChecker',
    'ConflictResolver',
    'get_registry',
    'reset_registry'
]
