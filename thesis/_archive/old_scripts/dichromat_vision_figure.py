#!/usr/bin/env python3
"""Generate side-by-side dichromatic vision comparison figure.

Creates a figure with:
- Left: Original tiger image (human trichromatic vision)
- Right: Dichromatic vision simulation (bovid prey vision)

This demonstrates why tigers evolved orange coloration - their prey
cannot distinguish orange from green due to having only two cone types.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Paths
THESIS_DIR = Path(__file__).resolve().parent.parent
TIGER_IMAGE = THESIS_DIR / "tiger" / "tiger_green_orange.png"
OUTPUT_DIR = THESIS_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def simulate_dichromat_vision(rgb_img: np.ndarray, sws_peak: float = 451, lws_peak: float = 555) -> np.ndarray:
    """Simulate dichromatic vision based on cone sensitivity peaks.

    Bovids have SWS (blue, ~451nm) and LWS (yellow-green, ~555nm) cones only.
    They cannot distinguish colors that differ only in the M-cone (green) region.
    Orange and green appear similar because both primarily stimulate LWS.
    """
    # Approximate cone response weights based on peak wavelengths
    # RGB roughly maps to: R~600nm, G~550nm, B~450nm

    # SWS cone response (primarily blue channel, some green)
    sws_weight_b = 1.0
    sws_weight_g = 0.1 * max(0, 1 - abs(sws_peak - 500) / 100)
    sws_weight_r = 0.0

    # LWS cone response (yellow-green: R and G channels)
    lws_weight_r = 0.7
    lws_weight_g = 1.0
    lws_weight_b = 0.0

    # Calculate cone responses
    r, g, b = rgb_img[:, :, 0], rgb_img[:, :, 1], rgb_img[:, :, 2]

    sws_response = sws_weight_r * r + sws_weight_g * g + sws_weight_b * b
    lws_response = lws_weight_r * r + lws_weight_g * g + lws_weight_b * b

    # Normalize responses
    sws_response = np.clip(sws_response, 0, 1)
    lws_response = np.clip(lws_response, 0, 1)

    # Map back to RGB for visualization
    # SWS -> Blue channel, LWS -> Yellow (R+G)
    out_r = lws_response * 0.8
    out_g = lws_response * 0.7 + sws_response * 0.1
    out_b = sws_response * 0.6

    return np.stack([out_r, out_g, out_b], axis=-1)


def main():
    """Generate the dichromatic vision comparison figure."""
    print("Loading tiger image...")
    img = Image.open(TIGER_IMAGE)
    img_array = np.array(img).astype(np.float32) / 255.0

    # Handle RGBA images
    if img_array.shape[-1] == 4:
        img_array = img_array[:, :, :3]

    print("Simulating dichromatic vision...")
    dichromat_view = simulate_dichromat_vision(img_array)
    dichromat_view = np.clip(dichromat_view, 0, 1)

    # Create side-by-side figure
    print("Creating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Left: Original (human vision)
    axes[0].imshow(img_array)
    axes[0].axis('off')

    # Right: Dichromatic vision
    axes[1].imshow(dichromat_view)
    axes[1].axis('off')

    plt.tight_layout()

    # Save figure
    output_path = OUTPUT_DIR / "dichromat_vision_tiger.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")

    # Also save as PDF for thesis
    output_pdf = OUTPUT_DIR / "dichromat_vision_tiger.pdf"
    fig.savefig(output_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_pdf}")

    plt.close()

    print("\nDone! The figure shows:")
    print("  Left:  Tiger as humans see it (orange stands out from green)")
    print("  Right: Tiger as prey animals see it (orange blends with green)")


if __name__ == "__main__":
    main()
