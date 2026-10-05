#!/usr/bin/env python3
"""Build the manuscript-ready dataset table and simple selection-flow exports."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor, white
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parent
REGISTRY = ROOT / "dataset_candidate_registry.csv"


def load_rows() -> list[dict[str, str]]:
    with REGISTRY.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def accession(row: dict[str, str]) -> str:
    values = [row[key] for key in ("GEO", "BioProject", "SRA") if row[key]]
    return "; ".join(values) if values else (f"PMID {row['PMID']}" if row["PMID"] else "none")


def groups(row: dict[str, str]) -> str:
    counts = f"control={row['sample_count_control']}; SD={row['sample_count_SD']}"
    return counts if row["other_groups"] in ("", "none") else f"{counts}; {row['other_groups']}"


def availability(row: dict[str, str]) -> str:
    return (
        f"raw={row['raw_data_available']}; processed={row['processed_data_available']}; "
        f"public={row['publicly_accessible']}"
    )


def exclusion(row: dict[str, str]) -> str:
    if row["eligibility_status"] == "EXCLUDED":
        return f"{row['exclusion_code']}: {row['exclusion_reason']}"
    if row["eligibility_status"] == "POTENTIALLY_ELIGIBLE":
        return row["notes"]
    return "Not excluded"


def write_tables(rows: list[dict[str, str]]) -> None:
    columns = [
        "Study",
        "Accession",
        "Species",
        "Sleep paradigm",
        "Tissue",
        "Assay",
        "Groups",
        "Replication",
        "Availability",
        "Eligibility",
        "Reason for exclusion",
    ]
    output_rows = []
    for row in rows:
        replication = (
            f"control={row['sample_count_control']}; SD={row['sample_count_SD']} biological samples"
        )
        output_rows.append(
            {
                "Study": row["title"],
                "Accession": accession(row),
                "Species": row["species"],
                "Sleep paradigm": f"{row['sleep_paradigm']}; duration={row['duration']}",
                "Tissue": f"{row['tissue']}; segment={row['tissue_segment']}",
                "Assay": row["assay"],
                "Groups": groups(row),
                "Replication": replication,
                "Availability": availability(row),
                "Eligibility": row["eligibility_status"],
                "Reason for exclusion": exclusion(row),
            }
        )

    csv_path = ROOT / "Supplementary_Table_dataset_search.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(output_rows)

    md_path = ROOT / "Supplementary_Table_dataset_search.md"
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("# Supplementary Table: public-dataset search and eligibility audit\n\n")
        handle.write("Search cutoff: **2026-10-04**. Unknown values were not inferred.\n\n")
        handle.write("| " + " | ".join(columns) + " |\n")
        handle.write("| " + " | ".join(["---"] * len(columns)) + " |\n")
        for row in output_rows:
            cells = [row[column].replace("|", "\\|").replace("\n", " ") for column in columns]
            handle.write("| " + " | ".join(cells) + " |\n")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    path = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf")
    if path.exists():
        return ImageFont.truetype(str(path), size)
    return ImageFont.truetype(name, size)


def wrapped(draw: ImageDraw.ImageDraw, text: str, box: tuple[int, int, int, int], size: int, bold: bool = False) -> None:
    x1, y1, x2, y2 = box
    fnt = font(size, bold)
    words = text.split()
    lines: list[str] = []
    line = ""
    for word in words:
        candidate = word if not line else f"{line} {word}"
        if draw.textbbox((0, 0), candidate, font=fnt)[2] <= x2 - x1 - 48:
            line = candidate
        else:
            lines.append(line)
            line = word
    if line:
        lines.append(line)
    line_height = size + 10
    total_height = len(lines) * line_height
    y = y1 + (y2 - y1 - total_height) / 2
    for line in lines:
        width = draw.textbbox((0, 0), line, font=fnt)[2]
        draw.text((x1 + (x2 - x1 - width) / 2, y), line, fill="#162337", font=fnt)
        y += line_height


def box(draw: ImageDraw.ImageDraw, coords: tuple[int, int, int, int], fill: str, outline: str, text: str, size: int = 34) -> None:
    draw.rounded_rectangle(coords, radius=24, fill=fill, outline=outline, width=4)
    wrapped(draw, text, coords, size=size)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int]) -> None:
    draw.line([start, end], fill="#5B6675", width=8)
    x, y = end
    draw.polygon([(x, y), (x - 16, y - 28), (x + 16, y - 28)], fill="#5B6675")


def write_png(rows: list[dict[str, str]]) -> None:
    status = Counter(row["eligibility_status"] for row in rows)
    reasons = Counter(row["exclusion_code"] for row in rows if row["eligibility_status"] == "EXCLUDED")
    image = Image.new("RGB", (2000, 1500), "#F7F9FC")
    draw = ImageDraw.Draw(image)
    title = "Public-dataset identification and selection audit"
    subtitle = "Search cutoff: 2026-10-04  |  Not a PRISMA systematic review"
    draw.text((100, 55), title, fill="#10243E", font=font(52, True))
    draw.text((100, 125), subtitle, fill="#53657A", font=font(29))

    box(draw, (310, 210, 1690, 420), "#E8F0FA", "#55789F",
        "Database searches and accession/citation checks: GEO, SRA, BioProject, PubMed, Europe PMC, ENA, and BioStudies/ArrayExpress", 35)
    arrow(draw, (1000, 420), (1000, 505))
    box(draw, (480, 505, 1520, 675), "#EDF4F1", "#4C7B6D",
        f"{len(rows)} unique candidate studies after study-level merging", 40)
    arrow(draw, (1000, 675), (1000, 760))

    draw.line([(1000, 760), (560, 760), (560, 835)], fill="#5B6675", width=8)
    draw.line([(1000, 760), (1440, 760), (1440, 835)], fill="#5B6675", width=8)
    draw.polygon([(560, 835), (544, 807), (576, 807)], fill="#5B6675")
    draw.polygon([(1440, 835), (1424, 807), (1456, 807)], fill="#5B6675")

    excluded_text = (
        f"{status['EXCLUDED']} EXCLUDED\n"
        f"E1={reasons['E1']}  E2={reasons['E2']}  E4={reasons['E4']}\n"
        f"E6={reasons['E6']}  E7={reasons['E7']}  E10={reasons['E10']}"
    )
    retained_text = (
        f"{status['ELIGIBLE']} ELIGIBLE\nGSE289089 / PRJNA1221288\nPRJNA1167170\nPRJNA862187\n\n"
        f"{status['POTENTIALLY_ELIGIBLE']} POTENTIALLY ELIGIBLE\nPRJNA1395688"
    )
    box(draw, (120, 835, 990, 1190), "#FBEDEE", "#B9686B", excluded_text, 35)
    box(draw, (1010, 835, 1880, 1190), "#EAF5ED", "#4A8660", retained_text, 34)

    note = "Duplicate accessions and publication records were merged at study level. No FASTQ files were downloaded and no candidate dataset was analysed."
    wrapped(draw, note, (170, 1260, 1830, 1415), size=29)
    image.save(ROOT / "dataset_selection_flow.png", dpi=(300, 300))


def pdf_text(c: canvas.Canvas, text: str, x: float, y: float, width: float, size: float, bold: bool = False) -> None:
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    words = text.split()
    lines: list[str] = []
    line = ""
    for word in words:
        candidate = word if not line else f"{line} {word}"
        if c.stringWidth(candidate, "Helvetica-Bold" if bold else "Helvetica", size) <= width:
            line = candidate
        else:
            lines.append(line)
            line = word
    if line:
        lines.append(line)
    for item in lines:
        c.drawCentredString(x + width / 2, y, item)
        y -= size * 1.25


def pdf_box(c: canvas.Canvas, x: float, y: float, width: float, height: float, fill: str, stroke: str, text: str, size: float) -> None:
    c.setFillColor(HexColor(fill))
    c.setStrokeColor(HexColor(stroke))
    c.setLineWidth(2)
    c.roundRect(x, y, width, height, 10, fill=1, stroke=1)
    c.setFillColor(HexColor("#162337"))
    pdf_text(c, text, x + 15, y + height / 2 + size, width - 30, size)


def write_pdf(rows: list[dict[str, str]]) -> None:
    status = Counter(row["eligibility_status"] for row in rows)
    reasons = Counter(row["exclusion_code"] for row in rows if row["eligibility_status"] == "EXCLUDED")
    path = ROOT / "dataset_selection_flow.pdf"
    width, height = landscape(letter)
    c = canvas.Canvas(str(path), pagesize=(width, height))
    c.setFillColor(HexColor("#F7F9FC"))
    c.rect(0, 0, width, height, fill=1, stroke=0)
    c.setFillColor(HexColor("#10243E"))
    c.setFont("Helvetica-Bold", 22)
    c.drawString(40, height - 42, "Public-dataset identification and selection audit")
    c.setFillColor(HexColor("#53657A"))
    c.setFont("Helvetica", 11)
    c.drawString(40, height - 62, "Search cutoff: 2026-10-04 | Not a PRISMA systematic review")

    pdf_box(c, 150, 395, 492, 90, "#E8F0FA", "#55789F", "Database searches and accession/citation checks: GEO, SRA, BioProject, PubMed, Europe PMC, ENA, and BioStudies/ArrayExpress", 12)
    c.setStrokeColor(HexColor("#5B6675")); c.setLineWidth(3); c.line(396, 395, 396, 355)
    pdf_box(c, 220, 285, 352, 70, "#EDF4F1", "#4C7B6D", f"{len(rows)} unique candidate studies after study-level merging", 14)
    c.line(396, 285, 396, 255); c.line(396, 255, 200, 255); c.line(396, 255, 592, 255); c.line(200, 255, 200, 230); c.line(592, 255, 592, 230)
    excluded = f"{status['EXCLUDED']} EXCLUDED | E1={reasons['E1']}, E2={reasons['E2']}, E4={reasons['E4']}, E6={reasons['E6']}, E7={reasons['E7']}, E10={reasons['E10']}"
    retained = f"{status['ELIGIBLE']} ELIGIBLE: GSE289089, PRJNA1167170, PRJNA862187 | {status['POTENTIALLY_ELIGIBLE']} POTENTIALLY ELIGIBLE: PRJNA1395688"
    pdf_box(c, 45, 110, 310, 120, "#FBEDEE", "#B9686B", excluded, 12)
    pdf_box(c, 437, 110, 310, 120, "#EAF5ED", "#4A8660", retained, 12)
    c.setFillColor(HexColor("#53657A"))
    pdf_text(c, "Duplicate accessions and publication records were merged at study level. No FASTQ files were downloaded and no candidate dataset was analysed.", 90, 75, 612, 10)
    c.save()


def main() -> None:
    rows = load_rows()
    assert len(rows) == 23
    assert Counter(row["eligibility_status"] for row in rows) == {
        "ELIGIBLE": 3,
        "POTENTIALLY_ELIGIBLE": 1,
        "EXCLUDED": 19,
    }
    write_tables(rows)
    write_png(rows)
    write_pdf(rows)


if __name__ == "__main__":
    main()
