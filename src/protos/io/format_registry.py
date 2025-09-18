"""
Centralized file format registry for Protos.

This module defines all supported file formats, their extensions, handlers,
and processor associations in one place to avoid duplication across the codebase.
"""

from typing import Dict, List, Set, Optional, Any
from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class ProcessorType(Enum):
    """Types of processors in Protos."""
    STRUCTURE = "structure"
    SEQUENCE = "sequence"
    GRN = "grn"
    PROPERTY = "property"
    EMBEDDING = "embedding"
    LIGAND = "ligand"
    GRAPH = "graph"


class FormatCategory(Enum):
    """Categories of file formats."""
    STRUCTURE = "structure"
    SEQUENCE = "sequence"
    TABLE = "table"
    EMBEDDING = "embedding"
    DATASET = "dataset"
    CONFIG = "config"
    ALIGNMENT = "alignment"
    NETWORK = "network"


@dataclass
class FileFormat:
    """Definition of a file format."""
    name: str                           # Format name (e.g., "PDB", "FASTA")
    extensions: List[str]               # File extensions (e.g., [".pdb"])
    category: FormatCategory            # Format category
    processor: ProcessorType            # Primary processor that handles this
    handler_class: Optional[str] = None # Handler class name in formats.py
    description: Optional[str] = None   # Human-readable description
    mime_type: Optional[str] = None     # MIME type if applicable
    binary: bool = False                # Whether format is binary
    
    @property
    def primary_extension(self) -> str:
        """Get the primary (preferred) extension."""
        return self.extensions[0] if self.extensions else ""


# Define all supported formats
FORMATS = {
    # Structure formats
    "PDB": FileFormat(
        name="PDB",
        extensions=[".pdb"],
        category=FormatCategory.STRUCTURE,
        processor=ProcessorType.STRUCTURE,
        handler_class="PDBHandler",
        description="Protein Data Bank format",
        mime_type="chemical/x-pdb"
    ),
    "MMCIF": FileFormat(
        name="mmCIF",
        extensions=[".cif", ".mmcif"],
        category=FormatCategory.STRUCTURE,
        processor=ProcessorType.STRUCTURE,
        handler_class="CIFHandler",
        description="Macromolecular Crystallographic Information File",
        mime_type="chemical/x-mmcif"
    ),
    "STRUCTURE_PKL": FileFormat(
        name="Structure Pickle",
        extensions=[".pkl"],
        category=FormatCategory.STRUCTURE,
        processor=ProcessorType.STRUCTURE,
        handler_class="PickleHandler",
        description="Cached structure DataFrame",
        binary=True
    ),
    
    # Sequence formats
    "FASTA": FileFormat(
        name="FASTA",
        extensions=[".fasta", ".fa", ".faa"],
        category=FormatCategory.SEQUENCE,
        processor=ProcessorType.SEQUENCE,
        handler_class="FASTAHandler",
        description="FASTA sequence format",
        mime_type="text/x-fasta"
    ),
    "CLUSTAL": FileFormat(
        name="Clustal",
        extensions=[".aln", ".clustal"],
        category=FormatCategory.ALIGNMENT,
        processor=ProcessorType.SEQUENCE,
        handler_class="ClustalHandler",
        description="Clustal alignment format"
    ),
    
    # Table formats
    "CSV": FileFormat(
        name="CSV",
        extensions=[".csv"],
        category=FormatCategory.TABLE,
        processor=ProcessorType.PROPERTY,  # Default, but used by multiple
        handler_class="CSVHandler",
        description="Comma-separated values",
        mime_type="text/csv"
    ),
    "TSV": FileFormat(
        name="TSV",
        extensions=[".tsv", ".tab"],
        category=FormatCategory.TABLE,
        processor=ProcessorType.PROPERTY,
        handler_class="TSVHandler",
        description="Tab-separated values",
        mime_type="text/tab-separated-values"
    ),
    "EXCEL": FileFormat(
        name="Excel",
        extensions=[".xlsx", ".xls"],
        category=FormatCategory.TABLE,
        processor=ProcessorType.PROPERTY,
        handler_class="ExcelHandler",
        description="Microsoft Excel spreadsheet",
        mime_type="application/vnd.ms-excel",
        binary=True
    ),
    
    # Embedding formats
    "NUMPY": FileFormat(
        name="NumPy",
        extensions=[".npy", ".npz"],
        category=FormatCategory.EMBEDDING,
        processor=ProcessorType.EMBEDDING,
        handler_class="NumpyHandler",
        description="NumPy array format",
        binary=True
    ),
    "PYTORCH": FileFormat(
        name="PyTorch",
        extensions=[".pt", ".pth"],
        category=FormatCategory.EMBEDDING,
        processor=ProcessorType.EMBEDDING,
        handler_class="PyTorchHandler",
        description="PyTorch tensor format",
        binary=True
    ),
    "EMBEDDING_PKL": FileFormat(
        name="Embedding Pickle",
        extensions=[".pkl"],
        category=FormatCategory.EMBEDDING,
        processor=ProcessorType.EMBEDDING,
        handler_class="PickleHandler",
        description="Pickled embeddings",
        binary=True
    ),
    
    # Dataset and config formats
    "JSON": FileFormat(
        name="JSON",
        extensions=[".json"],
        category=FormatCategory.CONFIG,
        processor=ProcessorType.STRUCTURE,  # Default
        handler_class="JSONHandler",
        description="JavaScript Object Notation",
        mime_type="application/json"
    ),
    "YAML": FileFormat(
        name="YAML",
        extensions=[".yaml", ".yml"],
        category=FormatCategory.CONFIG,
        processor=ProcessorType.STRUCTURE,
        handler_class="YAMLHandler",
        description="YAML configuration",
        mime_type="text/yaml"
    ),
    
    # Ligand formats
    "SDF": FileFormat(
        name="SDF",
        extensions=[".sdf", ".mol"],
        category=FormatCategory.STRUCTURE,
        processor=ProcessorType.LIGAND,
        handler_class="SDFHandler",
        description="Structure Data File for small molecules"
    ),
    "MOL2": FileFormat(
        name="MOL2",
        extensions=[".mol2"],
        category=FormatCategory.STRUCTURE,
        processor=ProcessorType.LIGAND,
        handler_class="Mol2Handler",
        description="Tripos MOL2 format"
    ),
    
    # Network formats
    "GRAPHML": FileFormat(
        name="GraphML",
        extensions=[".graphml"],
        category=FormatCategory.NETWORK,
        processor=ProcessorType.GRAPH,
        handler_class="GraphMLHandler",
        description="Graph Markup Language",
        mime_type="application/graphml+xml"
    ),
}


class FormatRegistry:
    """Registry for file format definitions and lookups."""
    
    def __init__(self):
        """Initialize the format registry."""
        self.formats = FORMATS
        self._build_indexes()
    
    def _build_indexes(self):
        """Build lookup indexes for fast access."""
        # Extension to format mapping
        self.extension_map: Dict[str, FileFormat] = {}
        for fmt in self.formats.values():
            for ext in fmt.extensions:
                # Handle case where multiple formats share extension
                if ext not in self.extension_map:
                    self.extension_map[ext] = fmt
        
        # Processor to formats mapping
        self.processor_map: Dict[ProcessorType, List[FileFormat]] = {}
        for fmt in self.formats.values():
            if fmt.processor not in self.processor_map:
                self.processor_map[fmt.processor] = []
            self.processor_map[fmt.processor].append(fmt)
        
        # Category to formats mapping
        self.category_map: Dict[FormatCategory, List[FileFormat]] = {}
        for fmt in self.formats.values():
            if fmt.category not in self.category_map:
                self.category_map[fmt.category] = []
            self.category_map[fmt.category].append(fmt)
    
    def get_format_by_extension(self, extension: str) -> Optional[FileFormat]:
        """
        Get format definition by file extension.
        
        Args:
            extension: File extension (with or without dot)
            
        Returns:
            FileFormat if found, None otherwise
        """
        if not extension.startswith('.'):
            extension = f'.{extension}'
        return self.extension_map.get(extension.lower())
    
    def get_format_by_path(self, path: Path) -> Optional[FileFormat]:
        """
        Get format definition by file path.
        
        Args:
            path: File path
            
        Returns:
            FileFormat if found, None otherwise
        """
        return self.get_format_by_extension(path.suffix)
    
    def get_formats_for_processor(self, processor: ProcessorType) -> List[FileFormat]:
        """
        Get all formats handled by a processor.
        
        Args:
            processor: Processor type
            
        Returns:
            List of FileFormat objects
        """
        return self.processor_map.get(processor, [])
    
    def get_extensions_for_processor(self, processor: ProcessorType) -> Set[str]:
        """
        Get all file extensions handled by a processor.
        
        Args:
            processor: Processor type
            
        Returns:
            Set of extensions
        """
        extensions = set()
        for fmt in self.get_formats_for_processor(processor):
            extensions.update(fmt.extensions)
        return extensions
    
    def get_formats_by_category(self, category: FormatCategory) -> List[FileFormat]:
        """
        Get all formats in a category.
        
        Args:
            category: Format category
            
        Returns:
            List of FileFormat objects
        """
        return self.category_map.get(category, [])
    
    def detect_format(self, path: Path, content: Optional[bytes] = None) -> Optional[FileFormat]:
        """
        Detect format by extension and optionally content.
        
        Args:
            path: File path
            content: Optional file content for validation
            
        Returns:
            FileFormat if detected, None otherwise
        """
        # First try by extension
        fmt = self.get_format_by_extension(path.suffix)
        
        # If ambiguous (e.g., .pkl used by multiple), try content detection
        if fmt and content is not None and path.suffix == '.pkl':
            # TODO: Implement content-based detection for ambiguous formats
            pass
        
        return fmt
    
    def get_processor_for_file(self, path: Path) -> Optional[ProcessorType]:
        """
        Get the appropriate processor type for a file.
        
        Args:
            path: File path
            
        Returns:
            ProcessorType if found, None otherwise
        """
        fmt = self.get_format_by_path(path)
        return fmt.processor if fmt else None


# Global registry instance
format_registry = FormatRegistry()


# Convenience functions
def get_file_format(path: Path) -> Optional[FileFormat]:
    """Get format definition for a file path."""
    return format_registry.get_format_by_path(path)


def get_processor_type(path: Path) -> Optional[ProcessorType]:
    """Get processor type for a file path."""
    return format_registry.get_processor_for_file(path)


def get_valid_extensions(processor: ProcessorType) -> Set[str]:
    """Get valid file extensions for a processor."""
    return format_registry.get_extensions_for_processor(processor)