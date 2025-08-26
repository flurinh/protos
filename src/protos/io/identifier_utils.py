"""
Identifier resolution utilities for the Protos dual ID system.

This module provides utility functions for resolving between human-readable
identifiers and hash-based entity IDs. It serves as a convenience layer
on top of the EntityRegistry.
"""

from typing import Optional, List, Dict, Union, Tuple
from pathlib import Path
import logging

from .data_access import EntityRegistry, generate_entity_id

logger = logging.getLogger(__name__)


def get_entity_registry(registry_path: Optional[str] = None) -> EntityRegistry:
    """
    Get or create the entity registry instance.
    
    Args:
        registry_path: Optional path to registry file
        
    Returns:
        EntityRegistry instance
    """
    if registry_path is None:
        registry_path = Path("data") / "entity_registry.json"
    
    return EntityRegistry(str(registry_path))


def resolve_to_entity_id(
    identifier: str,
    format_type: Optional[str] = None,
    registry: Optional[EntityRegistry] = None,
    create_if_missing: bool = True
) -> str:
    """
    Resolve any identifier to an entity hash ID.
    
    This is the primary function for converting human-readable identifiers
    to hash IDs. It handles:
    - Existing hash IDs (validates and returns them)
    - Human-readable IDs (looks up in registry)
    - Aliases (alternative names for entities)
    - New identifiers (creates hash if create_if_missing=True)
    
    Args:
        identifier: Any identifier (PDB ID, UniProt ID, custom name, or hash)
        format_type: Optional format type for disambiguation
        registry: Optional EntityRegistry instance (creates one if not provided)
        create_if_missing: If True, generates hash for unknown identifiers
        
    Returns:
        Entity hash ID (10-character alphanumeric string)
        
    Examples:
        >>> resolve_to_entity_id("1ubq")
        '36c2c0da93'
        
        >>> resolve_to_entity_id("P62988")  # UniProt ID for same protein
        '36c2c0da93'  # Same hash if aliased
        
        >>> resolve_to_entity_id("36c2c0da93")  # Already a hash
        '36c2c0da93'
    """
    if registry is None:
        registry = get_entity_registry()
    
    # Use registry's resolve_identifier method
    entity_id = registry.resolve_identifier(identifier, format_type)
    
    # If not found and create_if_missing is False, check if it exists
    if not create_if_missing and not registry.entity_exists(entity_id):
        # Try to find by original ID
        found_id = registry.find_entity_by_original_id(identifier, format_type)
        if found_id:
            return found_id
        raise ValueError(f"Entity not found for identifier: {identifier}")
    
    return entity_id


def resolve_to_human_id(
    identifier: str,
    format_type: Optional[str] = None,
    registry: Optional[EntityRegistry] = None,
    prefer_original: bool = True
) -> str:
    """
    Resolve any identifier to a human-readable ID.
    
    This function converts hash IDs back to human-readable identifiers,
    or validates existing human-readable IDs.
    
    Args:
        identifier: Any identifier (hash or human-readable)
        format_type: Optional format type for filtering
        registry: Optional EntityRegistry instance
        prefer_original: If True, returns original_id; if False, may return alias
        
    Returns:
        Human-readable identifier
        
    Raises:
        ValueError: If entity not found or has no human-readable ID
        
    Examples:
        >>> resolve_to_human_id("36c2c0da93")
        '1ubq'
        
        >>> resolve_to_human_id("1ubq")
        '1ubq'
    """
    if registry is None:
        registry = get_entity_registry()
    
    # First resolve to entity ID
    entity_id = resolve_to_entity_id(identifier, format_type, registry, create_if_missing=False)
    
    # Get human ID
    human_id = registry.get_human_id(entity_id)
    if human_id:
        return human_id
    
    # If no original ID and not preferring original, check aliases
    if not prefer_original:
        aliases = registry.get_aliases(entity_id)
        if aliases:
            return aliases[0]
    
    # Last resort - return the entity ID itself
    # This handles cases where entities are registered without human IDs
    return entity_id


def resolve_identifiers_batch(
    identifiers: List[str],
    format_type: Optional[str] = None,
    registry: Optional[EntityRegistry] = None,
    to_format: str = "hash"
) -> Dict[str, str]:
    """
    Batch resolve multiple identifiers.
    
    Args:
        identifiers: List of identifiers to resolve
        format_type: Optional format type for all identifiers
        registry: Optional EntityRegistry instance
        to_format: Target format - "hash" or "human"
        
    Returns:
        Dictionary mapping input identifiers to resolved identifiers
        
    Examples:
        >>> resolve_identifiers_batch(["1ubq", "2gb1", "P62988"])
        {'1ubq': '36c2c0da93', '2gb1': 'ba1837f945', 'P62988': '36c2c0da93'}
    """
    if registry is None:
        registry = get_entity_registry()
    
    results = {}
    
    if to_format == "hash":
        for identifier in identifiers:
            try:
                results[identifier] = resolve_to_entity_id(identifier, format_type, registry)
            except Exception as e:
                logger.warning(f"Failed to resolve {identifier}: {e}")
                results[identifier] = None
    
    elif to_format == "human":
        for identifier in identifiers:
            try:
                results[identifier] = resolve_to_human_id(identifier, format_type, registry)
            except Exception as e:
                logger.warning(f"Failed to resolve {identifier}: {e}")
                results[identifier] = None
    
    else:
        raise ValueError(f"Invalid to_format: {to_format}. Use 'hash' or 'human'")
    
    return results


def is_hash_id(identifier: str) -> bool:
    """
    Check if an identifier appears to be a hash ID.
    
    Args:
        identifier: Identifier to check
        
    Returns:
        True if identifier looks like a hash ID (10 alphanumeric chars)
    """
    return len(identifier) == 10 and identifier.isalnum() and identifier.islower()


def is_registered(
    identifier: str,
    format_type: Optional[str] = None,
    registry: Optional[EntityRegistry] = None
) -> bool:
    """
    Check if an identifier is registered in the entity registry.
    
    Args:
        identifier: Any identifier (hash or human-readable)
        format_type: Optional format type to check
        registry: Optional EntityRegistry instance
        
    Returns:
        True if the entity is registered
    """
    if registry is None:
        registry = get_entity_registry()
    
    try:
        entity_id = resolve_to_entity_id(identifier, format_type, registry, create_if_missing=False)
        return registry.entity_exists(entity_id)
    except:
        return False


def get_entity_formats(
    identifier: str,
    registry: Optional[EntityRegistry] = None
) -> List[str]:
    """
    Get all available formats for an entity.
    
    Args:
        identifier: Any identifier (hash or human-readable)
        registry: Optional EntityRegistry instance
        
    Returns:
        List of format types available for the entity
    """
    if registry is None:
        registry = get_entity_registry()
    
    try:
        entity_id = resolve_to_entity_id(identifier, None, registry, create_if_missing=False)
        return registry.get_entity_formats(entity_id)
    except:
        return []


def get_all_identifiers(
    identifier: str,
    registry: Optional[EntityRegistry] = None
) -> Dict[str, Union[str, List[str]]]:
    """
    Get all identifiers (hash, original, aliases) for an entity.
    
    Args:
        identifier: Any identifier for the entity
        registry: Optional EntityRegistry instance
        
    Returns:
        Dictionary with 'hash', 'original', and 'aliases' keys
    """
    if registry is None:
        registry = get_entity_registry()
    
    try:
        entity_id = resolve_to_entity_id(identifier, None, registry, create_if_missing=False)
        
        return {
            'hash': entity_id,
            'original': registry.get_original_id(entity_id),
            'aliases': registry.get_aliases(entity_id)
        }
    except:
        return {'hash': None, 'original': None, 'aliases': []}


def create_entity_id(
    identifier: str,
    check_existing: bool = True,
    registry: Optional[EntityRegistry] = None
) -> Tuple[str, bool]:
    """
    Create an entity ID for an identifier.
    
    Args:
        identifier: Human-readable identifier
        check_existing: If True, checks if entity already exists
        registry: Optional EntityRegistry instance
        
    Returns:
        Tuple of (entity_id, is_new) where is_new indicates if this is a new entity
    """
    if registry is None:
        registry = get_entity_registry()
    
    # Generate the hash ID
    entity_id = generate_entity_id(identifier)
    
    is_new = True
    if check_existing:
        # Check if it already exists
        if registry.entity_exists(entity_id):
            is_new = False
        else:
            # Also check if the identifier exists as an original ID or alias
            existing_id = registry.find_entity_by_original_id(identifier)
            if existing_id:
                entity_id = existing_id
                is_new = False
    
    return entity_id, is_new