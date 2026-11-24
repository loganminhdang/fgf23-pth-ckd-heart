"""Make a composite SVG fully vector by replacing its pasted panel pictures.

`submission_figures.py` builds Figure 1 and Supplementary Figure S1 by pasting
the panel PNGs written by `figure1.py` and `figureS1.py` (the TIF and PNG for
submission are rasters in any case). In the SVG copy, each pasted picture is
replaced here by the panel's own SVG from the same script, placed and scaled
to the same box, so that every line and label can be edited. Element ids of
each panel are prefixed with the panel name so that they stay unique.
"""
import re
from pathlib import Path

from lxml import etree

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
PT_PER_MM = 72 / 25.4


def _panel(path, prefix):
    """Root of the panel SVG with every id prefixed, and its viewBox."""
    root = etree.parse(str(path)).getroot()
    ids = {el.get("id") for el in root.iter() if isinstance(el.tag, str) and el.get("id")}

    def url(match):
        return f"url(#{prefix}{match.group(1)})" if match.group(1) in ids else match.group(0)

    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        for attr, value in list(el.attrib.items()):
            new = value
            if attr == "id":
                new = prefix + value
            elif attr == f"{{{XLINK}}}href" and value.startswith("#") and value[1:] in ids:
                new = "#" + prefix + value[1:]
            elif "url(#" in value:
                new = re.sub(r"url\(#([^)]+)\)", url, value)
            if new != value:
                el.set(attr, new)
    view_box = [float(v) for v in root.get("viewBox").split()]
    return root, view_box


def replace_pasted_panels(svg_path, placements):
    """Replace the pictures in svg_path, in order, by the panel SVGs.

    placements: [(panel_svg_path, x0_mm, top_mm, width_mm)], as given to paste().
    """
    svg_path = Path(svg_path)
    tree = etree.parse(str(svg_path))
    images = list(tree.getroot().iter(f"{{{SVG}}}image"))
    if len(images) != len(placements):
        raise ValueError(f"{svg_path.name}: {len(images)} pictures, {len(placements)} panels")
    for image, (panel_path, x0, top, width) in zip(images, placements):
        panel_path = Path(panel_path)
        if abs(float(image.get("x")) - x0 * PT_PER_MM) > 1.0:
            raise ValueError(f"{svg_path.name}: picture of {panel_path.stem} is not at x = {x0} mm")
        root, (vx, vy, vw, vh) = _panel(panel_path, f"{panel_path.stem}_")
        scale = width * PT_PER_MM / vw
        group = etree.Element(f"{{{SVG}}}g", nsmap={None: SVG, "xlink": XLINK})
        group.set("id", f"panel_{panel_path.stem}")
        group.set("transform", f"translate({x0 * PT_PER_MM:.6f} {top * PT_PER_MM:.6f}) "
                               f"scale({scale:.6f}) translate({-vx:g} {-vy:g})")
        for child in root:
            if isinstance(child.tag, str) and child.tag != f"{{{SVG}}}metadata":
                group.append(child)
        # The picture sits alone in a group clipped to its axes; replace that group.
        target = image
        parent = image.getparent()
        if parent.get("clip-path") and len(parent) == 1:
            target = parent
        target.getparent().replace(target, group)
    tree.write(str(svg_path), xml_declaration=True, encoding="utf-8", standalone=False)
