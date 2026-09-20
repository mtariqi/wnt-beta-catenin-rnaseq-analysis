#!/usr/bin/env python3
"""Upgrade scripts/17_tf_target_deg_spider_maps.py from star to spider-web layout."""

from pathlib import Path


TARGET = Path("scripts/17_tf_target_deg_spider_maps.py")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Patch '{label}' expected exactly one match; found {count}.")
    return text.replace(old, new, 1)


def main() -> None:
    if not TARGET.exists():
        raise FileNotFoundError(f"Run from the repository root; missing {TARGET}")
    text = TARGET.read_text(encoding="utf-8")
    if "def draw_spider_web" in text:
        print(f"{TARGET} already contains the spider-web upgrade; no changes made.")
        return

    text = replace_once(
        text,
        'GREY = "#737373"\n',
        'GREY = "#737373"\nWEB_GREY = "#D5D9DE"\n\n'
        '# Canonical display label for the RNF2/RING1B synonym pair.\n'
        'REGULATOR_ALIASES = {"RNF2": "RING1B"}\n',
        "constants",
    )

    text = replace_once(
        text,
        '    full["regulator"] = full.get(\n'
        '        "regulator", full["Term"].map(regulator_from_term)\n'
        '    ).astype(str).str.upper()\n',
        '    full["regulator"] = full.get(\n'
        '        "regulator", full["Term"].map(regulator_from_term)\n'
        '    ).astype(str).str.upper().replace(REGULATOR_ALIASES)\n',
        "alias normalization",
    )

    marker = '\n\ndef draw_map(row: pd.Series) -> plt.Figure:\n'
    web_function = '''

def draw_spider_web(ax: plt.Axes, genes: list[str]) -> None:
    """Draw polygonal guide rings; these are not biological gene-gene edges."""
    angles = np.linspace(0, 2 * math.pi, len(genes), endpoint=False) + math.pi / 2
    for radius, width, alpha in [(0.34, 1.0, 0.60), (0.67, 1.2, 0.75), (1.0, 1.5, 0.90)]:
        xs = [radius * math.cos(angle) for angle in angles]
        ys = [radius * math.sin(angle) for angle in angles]
        ax.plot(xs + [xs[0]], ys + [ys[0]], color=WEB_GREY,
                linewidth=width, alpha=alpha, zorder=0)
    for angle in angles:
        ax.plot([0, math.cos(angle)], [0, math.sin(angle)],
                color=WEB_GREY, linewidth=0.9, alpha=0.60, zorder=0)


def draw_map(row: pd.Series) -> plt.Figure:
'''
    text = replace_once(text, marker, web_function, "web function")

    text = replace_once(
        text,
        '    pos = radial_positions(tf, genes)\n\n'
        '    fig, ax = plt.subplots(figsize=(13.33, 7.5), constrained_layout=False)\n',
        '    pos = radial_positions(tf, genes)\n'
        '    hub_size = 4200 + max(0, len(tf) - 4) * 1100\n\n'
        '    fig, ax = plt.subplots(figsize=(13.33, 7.5), constrained_layout=False)\n',
        "hub size",
    )

    text = replace_once(
        text,
        '    ax.set_facecolor("white")\n'
        '    nx.draw_networkx_edges(graph, pos, ax=ax, width=2.0, edge_color=GREY, alpha=0.85)\n',
        '    ax.set_facecolor("white")\n'
        '    draw_spider_web(ax, genes)\n'
        '    nx.draw_networkx_edges(graph, pos, ax=ax, width=2.4, edge_color=GREY, alpha=0.95)\n',
        "draw web",
    )

    text = replace_once(
        text,
        'graph, pos, nodelist=[tf], ax=ax, node_size=4200,',
        'graph, pos, nodelist=[tf], ax=ax, node_size=hub_size,',
        "dynamic hub",
    )

    text = replace_once(
        text,
        '        font_size=20, font_weight="bold", font_color="white",\n',
        '        font_size=18, font_weight="bold", font_color="white",\n',
        "hub font",
    )

    text = replace_once(
        text,
        '    direction_label = direction.replace("regulated", "regulated candidate-gene")\n'
        '    fig.suptitle(\n'
        '        f"{tf} target-set overlap: Top-2,000 {direction_label} input",\n'
        '        fontsize=23, fontweight="bold", y=0.965,\n'
        '    )\n',
        '    direction_label = direction.replace("regulated", "regulated candidate genes")\n'
        '    fig.suptitle(\n'
        '        f"{tf}-centered target-gene spider map - Top-2,000 {direction_label}",\n'
        '        fontsize=20, fontweight="bold", y=0.965,\n'
        '    )\n',
        "title",
    )

    text = replace_once(
        text,
        '        "Edges indicate Enrichr target-set overlap only; they do not establish TF activity, direct binding, direction, or causality.",\n'
        '        ha="center", fontsize=12.5, color="#4B5563",\n',
        '        "Spokes indicate Enrichr overlap; pale web rings are layout guides only. "\n'
        '        "No TF activity, direct binding, gene-to-gene link, or causality is established.",\n'
        '        ha="center", fontsize=11.5, color="#4B5563",\n',
        "caption",
    )

    TARGET.write_text(text, encoding="utf-8")
    print(f"Upgraded {TARGET} to the spider-web layout.")
    print("Next: compile, verify markers, and regenerate the PDF.")


if __name__ == "__main__":
    main()
