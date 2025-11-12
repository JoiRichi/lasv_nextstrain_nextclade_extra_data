
import argparse
import json
import re
from collections import Counter
from pathlib import Path
from typing import List, Tuple

# Default palette (pair-list, Auspice-compatible)
DEFAULT_SCALE = [
    ("Lassa", "#2812b8"),
    ("Lassa_original", "#7f7f7f"),
    ("Gairo", "#2ca02c"),
    ("Luna", "#17becf"),
    ("Mopeia/Morogoro", "#ff7f0e"),
    ("Choriomeningitidis (LCMV)", "#9467bd"),
    ("Rat mammarenavirus", "#8c564b"),
    ("Solweziense", "#e377c2"),
    ("Dhati-welense", "#bcbd22"),
    ("Mammarenavirus sp.", "#1f9e89"),
]

def build_keyword_map(map_args: List[str]) -> List[Tuple[str, List[str]]]:
    mapping: List[Tuple[str, List[str]]] = []
    for item in map_args or []:
        if ":" not in item:
            raise ValueError(f"--map '{item}' must look like Label:kw1,kw2")
        label, kws = item.split(":", 1)
        keywords = [k.strip().lower() for k in kws.split(",") if k.strip()]
        if not label.strip() or not keywords:
            raise ValueError(f"--map '{item}' has empty label/keywords")
        mapping.append((label.strip(), keywords))
    return mapping

def compile_patterns(mapping: List[Tuple[str, List[str]]]) -> List[Tuple[str, List[re.Pattern]]]:
    compiled = []
    for label, kws in mapping:
        pats = [re.compile(rf"(?<![A-Za-z]){re.escape(kw)}(?![A-Za-z])", re.IGNORECASE) for kw in kws]
        compiled.append((label, pats))
    return compiled

def infer_from_text(text: str, patterns: List[Tuple[str, List[re.Pattern]]], unknown_label: str) -> str:
    t = (text or "")
    for label, pats in patterns:
        for p in pats:
            if p.search(t):
                return label
    return unknown_label

def extract_text_from_node(node: dict) -> str:
    na = node.get("node_attrs", {})
    parts = []
    name = na.get("name", {})
    if isinstance(name, dict):
        parts.append(name.get("value", ""))
    elif isinstance(name, str):
        parts.append(name)
    if isinstance(node.get("name"), str):
        parts.append(node["name"])
    for k in ("strain", "virus", "organism", "clade_membership"):
        v = na.get(k, {})
        if isinstance(v, dict):
            parts.append(v.get("value", ""))
        elif isinstance(v, str):
            parts.append(v)
    return " | ".join([p for p in parts if p])

def traverse_annotate(node: dict, patterns: List[Tuple[str, List[re.Pattern]]], unknown_label: str, species_counts: Counter) -> str:
    na = node.setdefault("node_attrs", {})
    text = extract_text_from_node(node)
    species = infer_from_text(text, patterns, unknown_label)

    children = node.get("children") or []
    child_species = []
    for ch in children:
        child_species.append(traverse_annotate(ch, patterns, unknown_label, species_counts))

    if species == unknown_label and child_species:
        counts = Counter(child_species)
        species = counts.most_common(1)[0][0]

    na["pathogen_species"] = {"value": species}

    if not children:
        species_counts[species] += 1
    return species

def ensure_coloring(meta: dict, title: str, scale_pairs: List[Tuple[str, str]]):
    colorings = meta.setdefault("colorings", [])
    colorings[:] = [c for c in colorings if c.get("key") != "pathogen_species"]
    colorings.append({
        "type": "categorical",
        "key": "pathogen_species",
        "title": title,
        "scale": [[lab, col] for (lab, col) in scale_pairs]
    })

def ensure_filter(meta: dict):
    filters = meta.setdefault("filters", [])
    if "pathogen_species" not in filters and not any(isinstance(f, dict) and f.get("key") == "pathogen_species" for f in filters):
        filters.append("pathogen_species")

def mirror_to_clade(node: dict):
    na = node.setdefault("node_attrs", {})
    ps = na.get("pathogen_species", {}).get("value")
    if ps:
        na["clade_membership"] = {"value": ps}
    for ch in node.get("children") or []:
        mirror_to_clade(ch)

def run(in_path: Path, out_path: Path, title: str, unknown_label: str,
        mirror_to_clade_flag: bool, add_filter_flag: bool,
        map_args: List[str], scale_args: List[str]):
    with in_path.open("r") as f:
        data = json.load(f)

    default_map = [
        ("Lassa", ["lassa", "lassaense"]),
        ("Gairo", ["gairo", "gairoense"]),
        ("Luna", ["luna", "lunaense"]),
        ("Mopeia/Morogoro", ["mopeia", "mopeiaense", "morogoro"]),
        ("Choriomeningitidis (LCMV)", ["choriomeningitidis", "lcmv", "lymphocytic choriomeningitis"]),
        ("Rat mammarenavirus", ["rat mammarenavirus", "rat-mammarenavirus", "ratmammarenavirus"]),
        ("Solweziense", ["solweziense", "solwezi"]),
        ("Dhati-welense", ["dhati-welense", "dhat-welense", "dhatiwelense"]),
        ("Mammarenavirus sp.", ["mammarenavirus sp", "mammarenavirus sp.", "sp. mammarenavirus", "mammarenavirus sp "])
    ]
    user_map = build_keyword_map(map_args) if map_args else []
    combined_map = user_map + default_map
    patterns = compile_patterns(combined_map)

    species_counts = Counter()

    if "tree" in data:
        traverse_annotate(data["tree"], patterns, unknown_label, species_counts)
    elif "trees" in data and isinstance(data["trees"], list):
        for t in data["trees"]:
            traverse_annotate(t, patterns, unknown_label, species_counts)
    else:
        raise ValueError("Unrecognized Auspice JSON: expected 'tree' or 'trees'.")

    # Color scale
    if scale_args:
        pairs: List[Tuple[str, str]] = []
        for s in scale_args:
            if ":" not in s:
                raise ValueError(f"--scale '{s}' must look like Label:#RRGGBB")
            lab, col = s.split(":", 1)
            lab, col = lab.strip(), col.strip()
            if not lab or not re.match(r"^#([0-9a-fA-F]{6})$", col):
                raise ValueError(f"--scale '{s}' has invalid label or color")
            pairs.append((lab, col))
        scale_pairs = pairs
    else:
        scale_dict = dict(DEFAULT_SCALE)
        if unknown_label not in scale_dict:
            scale_dict[unknown_label] = "#7f7f7f"
        order = ["Lassa", unknown_label, "Gairo", "Luna", "Mopeia/Morogoro",
                 "Choriomeningitidis (LCMV)", "Rat mammarenavirus", "Solweziense", "Dhati-welense"]
        for lab in species_counts.keys():
            if lab not in order:
                order.append(lab)
        scale_pairs = [(lab, scale_dict.get(lab, "#999999")) for lab in order]

    meta = data.setdefault("meta", {})
    ensure_coloring(meta, title, scale_pairs)

    if add_filter_flag:
        ensure_filter(meta)

    dd = meta.setdefault("display_defaults", {})
    dd["color_by"] = "pathogen_species"

    if mirror_to_clade_flag:
        if "tree" in data:
            mirror_to_clade(data["tree"])
        elif "trees" in data and data["trees"]:
            for t in data["trees"]:
                mirror_to_clade(t)

    with out_path.open("w") as f:
        json.dump(data, f, indent=2)

    return species_counts

def main():
    ap = argparse.ArgumentParser(description="Annotate Auspice JSON with pathogen species & color scale.")
    ap.add_argument("--in", dest="in_path", required=True, help="Input Auspice JSON")
    ap.add_argument("--out", dest="out_path", required=True, help="Output Auspice JSON")
    ap.add_argument("--title", default="Pathogen species", help="Legend/title for the coloring")
    ap.add_argument("--unknown-label", default="Lassa_original", help="Label to use when no keywords match")
    ap.add_argument("--mirror-to-clade", action="store_true", help="Also mirror into clade_membership to color by Clade")
    ap.add_argument("--add-filter", action="store_true", help="Add a sidebar filter for pathogen_species")
    ap.add_argument("--map", action="append", help="Keyword map like 'Label:kw1,kw2'. Repeatable.")
    ap.add_argument("--scale", action="append", help="Explicit color scale 'Label:#RRGGBB'. Repeatable.")
    args = ap.parse_args()

    in_path = Path(args.in_path)
    out_path = Path(args.out_path)

    species_counts = run(
        in_path=in_path,
        out_path=out_path,
        title=args.title,
        unknown_label=args.unknown_label,
        mirror_to_clade_flag=args.mirror_to_clade,
        add_filter_flag=args.add_filter,
        map_args=args.map,
        scale_args=args.scale
    )

    print(f"Annotated pathogen species written to: {out_path}")
    print("Counts (tips):")
    for lab, n in species_counts.most_common():
        print(f"  {lab}: {n}")

if __name__ == "__main__":
    main()
