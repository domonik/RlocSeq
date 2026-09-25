#!/usr/bin/env python3

"""Calculate mean 16S rRNA TAMRA intensity for each segmented cell."""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


TAMRA_CHANNEL = 1
OUTPUT_COLUMNS = [
    "replicate",
    "timepoint",
    "cellidx",
    "imagename",
    "meanintensity",
]

IMAGE_NUMBER_PATTERN = re.compile(r"^(?:Composite|TM)_Image (?P<number>\d+)")
REPLICATE_PATTERN = re.compile(r"^Replicate\s*(?P<number>\d+)$", re.IGNORECASE)
TIMEPOINT_PATTERN = re.compile(
    r"^(?P<hours>\d+)h_(?P<condition>light|dark)$", re.IGNORECASE
)


def _classify_file(path):
    """Return the role of a relevant TIFF, or None for unrelated TIFFs."""
    name = path.name
    if name.startswith("Composite_Image") and "corrected" not in name.lower():
        return "intensity"
    if name.startswith("TM_Image") and (
        "{TotalCell}" in name or "_corrected_total_cell" in name
    ):
        return "mask"
    return None


def _image_number(path):
    match = IMAGE_NUMBER_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Cannot extract an image number from {path}")
    return int(match.group("number"))


def _replicate_label(part):
    match = REPLICATE_PATTERN.match(part)
    if match is None:
        raise ValueError(f"Cannot extract a replicate number from {part!r}")
    return f"Replicate {int(match.group('number'))}"


def _sample_metadata(path):
    """Extract normalized replicate and timepoint labels from an input path."""
    parts = path.parts

    if "16sRNAFISHLightDark" in parts:
        root_index = parts.index("16sRNAFISHLightDark")
        try:
            replicate_part = parts[root_index + 1]
            timepoint_part = parts[root_index + 2]
        except IndexError as error:
            raise ValueError(f"Unexpected light/dark input path: {path}") from error

        timepoint_match = TIMEPOINT_PATTERN.match(timepoint_part)
        if timepoint_match is None:
            raise ValueError(f"Cannot parse light/dark timepoint from {path}")
        timepoint = (
            f"{int(timepoint_match.group('hours'))}h_"
            f"{timepoint_match.group('condition').lower()}"
        )
        return _replicate_label(replicate_part), timepoint

    if "AllFISHData" in parts:
        root_index = parts.index("AllFISHData")
        try:
            probe_part = parts[root_index + 1]
            replicate_part = parts[root_index + 2]
            treatment_part = parts[root_index + 3]
        except IndexError as error:
            raise ValueError(f"Unexpected AllFISHData input path: {path}") from error

        if probe_part != "16S rRNA" or treatment_part.lower() != "segmentation_untreated":
            raise ValueError(f"Input is not untreated 16S rRNA data: {path}")
        return _replicate_label(replicate_part), "continuous_light"

    raise ValueError(f"Input is outside the supported 16S rRNA datasets: {path}")


def _collect_image_pairs(input_paths):
    samples = {}
    for input_path in input_paths:
        path = Path(input_path)
        paths = sorted(path.glob("*.tif*")) if path.is_dir() else [path]
        for path in paths:
            role = _classify_file(path)
            if role is None:
                continue

            replicate, timepoint = _sample_metadata(path)
            image_number = _image_number(path)
            key = replicate, timepoint, image_number
            sample = samples.setdefault(key, {})
            if role in sample:
                raise ValueError(
                    f"Multiple {role} files found for {replicate}, {timepoint}, "
                    f"Image {image_number}: {sample[role]} and {path}"
                )
            sample[role] = path

    pairs = []
    for key, sample in sorted(samples.items()):
        if "intensity" not in sample:
            print(
                f"Skipping TotalCell mask without a composite image: {sample['mask']}",
                file=sys.stderr,
            )
            continue
        if "mask" not in sample:
            print(
                f"Skipping composite image without a TotalCell mask: "
                f"{sample['intensity']}",
                file=sys.stderr,
            )
            continue
        pairs.append((key, sample["intensity"], sample["mask"]))

    if not pairs:
        raise ValueError("No complete composite image and TotalCell mask pairs were found")
    return pairs


def _read_tamra_channel(path):
    with tifffile.TiffFile(path) as tif:
        image = tif.asarray()
        axes = tif.series[0].axes

    if "C" in axes:
        channel_axis = axes.index("C")
    elif image.ndim == 3 and image.shape[0] > TAMRA_CHANNEL:
        channel_axis = 0
    else:
        raise ValueError(
            f"Expected a multichannel composite image with a TAMRA channel in {path}; "
            f"found shape {image.shape} and axes {axes!r}"
        )

    if image.shape[channel_axis] <= TAMRA_CHANNEL:
        raise ValueError(
            f"Composite image {path} has no TAMRA channel at index {TAMRA_CHANNEL}"
        )

    tamra = np.take(image, TAMRA_CHANNEL, axis=channel_axis)
    tamra = np.squeeze(tamra)
    if tamra.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional TAMRA channel in {path}, found {tamra.shape}"
        )
    return tamra


def _mean_intensity_per_cell(intensity_path, mask_path):
    tamra = _read_tamra_channel(intensity_path)
    labels = np.squeeze(tifffile.imread(mask_path))

    if labels.ndim != 2:
        raise ValueError(f"Expected a two-dimensional TotalCell mask in {mask_path}")
    if labels.shape != tamra.shape:
        raise ValueError(
            f"Image/mask shape mismatch for {intensity_path} and {mask_path}: "
            f"{tamra.shape} versus {labels.shape}"
        )
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError(f"TotalCell mask labels must be integers in {mask_path}")

    foreground = labels > 0
    cell_labels, inverse = np.unique(labels[foreground], return_inverse=True)
    if cell_labels.size == 0:
        print(f"Skipping mask without labeled cells: {mask_path}", file=sys.stderr)
        return cell_labels, np.array([], dtype=float)

    pixel_counts = np.bincount(inverse)
    intensity_sums = np.bincount(
        inverse, weights=tamra[foreground].astype(float, copy=False)
    )
    return cell_labels, intensity_sums / pixel_counts


def calculate_mean_intensities(input_paths):
    rows = []
    for (replicate, timepoint, image_number), intensity_path, mask_path in (
        _collect_image_pairs(input_paths)
    ):
        cell_labels, means = _mean_intensity_per_cell(intensity_path, mask_path)
        for cell_label, mean_intensity in zip(cell_labels, means):
            rows.append(
                {
                    "replicate": replicate,
                    "timepoint": timepoint,
                    "cellidx": int(cell_label),
                    "imagename": f"Image {image_number}",
                    "meanintensity": float(mean_intensity),
                }
            )

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def main(input_paths, output_path):
    result = calculate_mean_intensities(input_paths)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake.input.data_dirs, snakemake.output.tsv)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("output", help="Output TSV path")
        parser.add_argument("inputs", nargs="+", help="Composite and TotalCell TIFFs")
        arguments = parser.parse_args()
        main(arguments.inputs, arguments.output)
