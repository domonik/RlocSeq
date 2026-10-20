import base64
from io import BytesIO
import re
import xml.etree.ElementTree as ET

from PIL import Image


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)


def is_square_image_element(elem):
    width = elem.get("width")
    height = elem.get("height")

    if width is None or height is None:
        return False

    try:
        w = float(re.sub(r"[^\d.]", "", width))
        h = float(re.sub(r"[^\d.]", "", height))
        return abs(w - h) < 1e-6
    except ValueError:
        return False


def encode_tif_to_data_uri(tif_path):
    """Embed a TIFF as lossless PNG at its original pixel dimensions."""
    with Image.open(tif_path) as image:
        if image.n_frames != 1:
            raise ValueError(f"Expected a single-frame FISH panel: {tif_path}")
        if image.mode not in {"1", "L", "LA", "P", "RGB", "RGBA", "I;16"}:
            raise ValueError(f"Cannot encode {image.mode} as PNG without conversion: {tif_path}")
        buffer = BytesIO()
        image.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def replace_images(svg_path, tif_list, output_svg, square_only=False):
    tree = ET.parse(svg_path)
    root = tree.getroot()

    namespaces = {"svg": SVG_NS}
    images = root.findall(".//svg:image", namespaces)

    tif_iter = iter(tif_list)
    replacements = 0

    for img in images:
        href = (
            img.get(f"{{{XLINK_NS}}}href")
            or img.get("href")
        )

        if not href:
            continue

        if square_only and not is_square_image_element(img):
            continue

        try:
            tif_path = next(tif_iter)
        except StopIteration:
            break

        new_data_uri = encode_tif_to_data_uri(tif_path)

        # Replace href while preserving attribute style
        if img.get(f"{{{XLINK_NS}}}href") is not None:
            img.set(f"{{{XLINK_NS}}}href", new_data_uri)
        else:
            img.set("href", new_data_uri)

        replacements += 1

    tree.write(output_svg, encoding="utf-8", xml_declaration=True)
    return replacements


def replace_square_images(svg_path, tif_list, output_svg):
    return replace_images(svg_path, tif_list, output_svg, square_only=True)


if __name__ == "__main__":
    svg_input = "Pipeline/Paper/Figure6/F6FISH.svg"
    svg_output = "Pipeline/Paper/Figure6/F6FISHfixed.svg"

    tif_images = [
        "initialData/AllFISHData/psbA2/Replicate 1/Composite_Image 36_untreated_green.tif",
        "initialData/AllFISHData/psbA2/Replicate 1/Composite_Image 54_Puromycin_green.tif",
        "initialData/AllFISHData/rpl10_12/Replicate 1/Kopie von Composite_Image 83_untreated_green.tif",
        "initialData/AllFISHData/rpl10_12/Replicate 1/Kopie von Composite_Image 95_Puromycin_green.tif"
    ]
    replace_square_images(svg_input, tif_images, svg_output)

    
