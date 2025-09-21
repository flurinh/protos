#!/usr/bin/env python3
"""Project data bootstrapper for Protos."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root / "src"))

# Project data directory
PROJECT_DATA_DIR = project_root / "data"


def reinstall_reference_data(paths, reason: str) -> None:
    marker_file = Path(paths.data_root) / ".protos_initialized"
    if marker_file.exists():
        marker_file.unlink()
    print(f"  ↻ Reinstalling GRN reference data ({reason})")
    paths._install_reference_data()  # pylint: disable=protected-access


def setup_project_data(refresh: bool = False) -> bool:
    """Set up the project data directory."""
    print("=== Protos Project Data Setup ===")
    print(f"Project root: {project_root}")
    print(f"Data directory: {PROJECT_DATA_DIR}\n")

    import protos
    from protos.io.paths import get_protos_paths

    print("Setting up data directory...")
    protos.set_data_path(str(PROJECT_DATA_DIR))
    paths = get_protos_paths()
    print(f"ProtosPaths configured with: {paths.data_root}")

    if refresh:
        reinstall_reference_data(paths, "--refresh requested")

    print("\nInitializing directory structure...")
    paths.get_processor_path("structure")

    if not PROJECT_DATA_DIR.exists():
        print("❌ Failed to create data directory!")
        return False

    print("✅ Data directory created successfully!")
    expected_dirs = [
        "structure",
        "grn",
        "sequence",
        "property",
        "embedding",
        "ligand",
        "graph",
        "input",
        "temp",
    ]

    all_exist = True
    for dir_name in expected_dirs:
        if (PROJECT_DATA_DIR / dir_name).exists():
            print(f"  ✓ {dir_name}/")
        else:
            print(f"  ✗ {dir_name}/ (missing)")
            all_exist = False

    registry_path = PROJECT_DATA_DIR / "global_registry.json"
    if registry_path.exists():
        print("  ✓ global_registry.json")
    else:
        print("  ✗ global_registry.json (missing)")
        all_exist = False

    grn_config = PROJECT_DATA_DIR / "grn" / "configs" / "config.json"
    grn_reference_dir = PROJECT_DATA_DIR / "grn" / "reference"
    grn_reference_files = list(grn_reference_dir.glob("*") if grn_reference_dir.exists() else [])

    if not grn_reference_files:
        reinstall_reference_data(paths, "missing grn/reference contents")
        grn_reference_files = list(grn_reference_dir.glob("*"))

    if grn_config.exists() and grn_reference_files:
        ref_names = ", ".join(f.name for f in grn_reference_files[:3])
        if len(grn_reference_files) > 3:
            ref_names += ", ..."
        print("\n✅ Reference data installed")
        print(f"  • Config: {grn_config.name}")
        print(f"  • Reference files: {ref_names}")
    else:
        print("\n⚠️  Reference data incomplete; please reinstall the package data files.")
        all_exist = False

    if all_exist:
        print("\n🎉 Project data setup complete!")
        print(f"📁 Location: {PROJECT_DATA_DIR}")
        print("\n📝 Example usage:")
        print("```python")
        print("from protos.processing.structure import StructureProcessor")
        print()
        print("# Create a processor - it will use the project data directory")
        print('processor = StructureProcessor("my_analysis")')
        print()
        print("# Load a structure")
        print('processor.load_structure("1ubq")')
        print("```")
    else:
        print("\n⚠️  Some directories or reference files were not created. Please check the output above.")

    return all_exist


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Initialize Protos project data directory.")
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Force reinstallation of bundled reference data (GRN configs/ref).",
    )
    args = parser.parse_args()

    try:
        success = setup_project_data(refresh=args.refresh)
        sys.exit(0 if success else 1)
    except Exception as exc:  # noqa: BLE001
        print(f"\n❌ Error: {exc}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
