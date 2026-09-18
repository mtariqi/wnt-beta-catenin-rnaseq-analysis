#!/usr/bin/env python3
"""Keep the TF in the center when it also occurs in its own overlap list."""

from pathlib import Path


TARGET = Path("scripts/17_tf_target_deg_spider_maps.py")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Patch '{label}' expected one match; found {count}.")
    return text.replace(old, new, 1)


def main() -> None:
    if not TARGET.exists():
        raise FileNotFoundError(f"Run from the repository root; missing {TARGET}")
    text = TARGET.read_text(encoding="utf-8")
    if "peripheral_genes = [gene for gene in all_genes if gene != tf]" in text:
        print("Central-TF fix is already installed; no changes made.")
        return
    if "def draw_spider_web" not in text:
        raise RuntimeError("Install the spider-web upgrade before applying this fix.")

    text = replace_once(
        text,
        '    all_genes = row["overlap_genes"]\n'
        '    genes = all_genes[:MAX_GENES_PER_MAP]\n',
        '    all_genes = row["overlap_genes"]\n'
        '    # A TF can occur in its own enrichment overlap. Keep that evidence\n'
        '    # in the table, but do not let it overwrite the central TF position.\n'
        '    peripheral_genes = [gene for gene in all_genes if gene != tf]\n'
        '    genes = peripheral_genes[:MAX_GENES_PER_MAP]\n',
        "drawing selection",
    )

    text = replace_once(
        text,
        '    if len(all_genes) > len(genes):\n'
        '        display_note = f"Showing {len(genes)} of {len(all_genes)} exact overlap genes; complete list is in the evidence table."\n',
        '    if len(peripheral_genes) > len(genes):\n'
        '        self_note = " (self-overlap retained in table)" if tf in all_genes else ""\n'
        '        display_note = (f"Showing {len(genes)} of {len(peripheral_genes)} peripheral overlap genes"\n'
        '                        f"{self_note}; complete list is in the evidence table.")\n',
        "caption",
    )

    text = replace_once(
        text,
        '            shown = row["overlap_genes"][:MAX_GENES_PER_MAP]\n',
        '            shown = [gene for gene in row["overlap_genes"] if gene != tf][:MAX_GENES_PER_MAP]\n',
        "manifest selection",
    )

    TARGET.write_text(text, encoding="utf-8")
    print(f"Installed central-TF fix in {TARGET}.")


if __name__ == "__main__":
    main()
