"""
Test if Enamine credentials are properly loaded from .env file.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

# First, let's check if .env exists
env_file = Path(__file__).parent / '.env'
print(f".env file exists: {env_file.exists()}")
print(f".env file path: {env_file.absolute()}")

# Import and check dotenv
try:
    from dotenv import load_dotenv
    print("\n✓ python-dotenv is installed")
    
    # Load .env file
    load_dotenv()
    print("✓ .env file loaded")
except ImportError:
    print("\n✗ python-dotenv is NOT installed")
    sys.exit(1)

# Check environment variables
print("\n=== Environment Variables ===")
username1 = os.environ.get('ENAMINE_USERNAME')
username2 = os.environ.get('enamine_username')
password1 = os.environ.get('ENAMINE_PASSWORD')
password2 = os.environ.get('enamine_password')

print(f"ENAMINE_USERNAME: {'SET' if username1 else 'NOT SET'}")
print(f"enamine_username: {'SET' if username2 else 'NOT SET'}")
print(f"ENAMINE_PASSWORD: {'SET' if password1 else 'NOT SET'}")
print(f"enamine_password: {'SET' if password2 else 'NOT SET'}")

# Now test through the loader
print("\n=== Testing through Enamine loader ===")
from protos.io.ingest import enamine_loader

username, password = enamine_loader.get_enamine_credentials()
print(f"Credentials found: {bool(username and password)}")
if username:
    print(f"Username: {username}")
else:
    print("No credentials found by loader")

# List available datasets
print("\n=== Available Datasets ===")
datasets = enamine_loader.list_available_datasets()
print(f"Found {len(datasets)} datasets")
for name in list(datasets.keys())[:3]:
    print(f"  - {name}")

print("\nNote: The Enamine download URLs in the loader are placeholders.")
print("Real Enamine URLs would need to be configured based on your subscription.")