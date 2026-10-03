#!/usr/bin/env python3
"""Summarize matched somatic VCFs and build the Somvar report using Python's standard library."""

import argparse
import csv
import gzip
import hashlib
import html
import json
import re
import shutil
from collections import Counter, defaultdict, deque
from datetime import datetime
from pathlib import Path


CORE = ["chromosome", "position", "ref", "alt", "filter", "variant_type",
        "annotation_status", "genes", "consequences", "impacts"]
SAMPLE_FIELDS = [f"{sample}_{field}" for sample in ("normal", "tumor")
                 for field in ("GT", "DP", "AD", "AF")]
IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def open_vcf(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_header(path):
    annotations, samples, version = [], [], None
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r"Format: ([^\"]+)", line)
                if match:
                    annotations = match.group(1).rstrip(">").split("|")
            elif line.startswith("##VEP="):
                version = line.strip()
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    return annotations, samples, version


def records(path):
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record in {path}")
            yield fields


def key(fields):
    return fields[0], fields[1], fields[3], fields[4]


def parse_csq(fields, names):
    csq = next((item[4:] for item in fields[7].split(";") if item.startswith("CSQ=")), "")
    if csq and not names:
        raise ValueError("CSQ values are present but the annotation header is missing")
    return [dict(zip(names, entry.split("|"))) for entry in csq.split(",") if entry and entry != "."]


def write_table(path, columns, rows, delimiter="\t"):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def variant_type(ref, alt):
    types = set()
    for allele in alt.split(","):
        if allele == "*" or allele.startswith("<") or "[" in allele or "]" in allele:
            types.add("other")
        elif len(ref) == len(allele) == 1:
            types.add("SNV")
        elif len(ref) != len(allele):
            types.add("indel")
        else:
            types.add("MNV")
    return next(iter(types)) if len(types) == 1 else "mixed"


def summarize(args):
    output = Path(args.outdir)
    output.mkdir(parents=True, exist_ok=True)
    names, _, vep_version = read_header(args.vep_vcf)
    _, samples, _ = read_header(args.pass_vcf)
    missing = {"normal", "tumor"} - set(samples)
    if missing:
        raise ValueError("PASS VCF is missing sample columns: " + ", ".join(sorted(missing)))
    annotations = defaultdict(deque)
    vep_records = 0
    for fields in records(args.vep_vcf):
        if not names:
            raise ValueError("Nonempty VEP VCF is missing its CSQ header")
        annotations[key(fields)].append(parse_csq(fields, names))
        vep_records += 1

    counts = Counter(pass_records=0, annotated_records=0, unannotated_output_records=0,
                     not_in_vep_output_records=0)
    genes, consequences, impacts, types = Counter(), Counter(), Counter(), Counter()
    transcript_columns = ["chromosome", "position", "ref", "alt"] + names
    with (output / "variants.tsv").open("w", newline="") as tsv, \
         (output / "variants.csv").open("w", newline="") as csv_handle, \
         (output / "consequences.tsv").open("w", newline="") as transcripts:
        writers = [csv.DictWriter(tsv, fieldnames=CORE + SAMPLE_FIELDS, delimiter="\t",
                                  lineterminator="\n"),
                   csv.DictWriter(csv_handle, fieldnames=CORE + SAMPLE_FIELDS, lineterminator="\n")]
        transcript_writer = csv.DictWriter(transcripts, fieldnames=transcript_columns,
                                            delimiter="\t", lineterminator="\n")
        for writer in [*writers, transcript_writer]:
            writer.writeheader()
        for fields in records(args.pass_vcf):
            if fields[6] != "PASS":
                raise ValueError("PASS input contains a record whose FILTER is not PASS")
            counts["pass_records"] += 1
            group = annotations.get(key(fields))
            if group:
                csq = group.popleft()
                status = "annotated" if csq else "unannotated_output"
            else:
                csq, status = [], "not_in_vep_output"
            counts[status + "_records"] += 1
            gene_set = {a.get("SYMBOL") or a.get("Gene") for a in csq} - {None, ""}
            consequence_set = {c for a in csq for c in a.get("Consequence", "").split("&") if c}
            impact_set = {a.get("IMPACT") for a in csq} - {None, ""}
            genes.update(gene_set)
            consequences.update(consequence_set)
            if impact_set:
                impacts[min(impact_set, key=lambda value: IMPACT_ORDER.get(value, 4))] += 1
            kind = variant_type(fields[3], fields[4])
            types[kind] += 1
            row = dict(zip(CORE[:5], [fields[0], fields[1], fields[3], fields[4], fields[6]]))
            row.update(variant_type=kind, annotation_status=status, genes=";".join(sorted(gene_set)),
                       consequences=";".join(sorted(consequence_set)), impacts=";".join(sorted(impact_set)))
            format_names = fields[8].split(":") if len(fields) > 8 else []
            for sample in ("normal", "tumor"):
                index = 9 + samples.index(sample)
                values = dict(zip(format_names, fields[index].split(":"))) if index < len(fields) else {}
                row.update({f"{sample}_{name}": values.get(name, ".") for name in ("GT", "DP", "AD", "AF")})
            for writer in writers:
                writer.writerow(row)
            for annotation in csq:
                transcript_writer.writerow({**dict(zip(transcript_columns[:4], key(fields))), **annotation})
    if any(group for group in annotations.values()):
        raise ValueError("VEP contains records that could not be matched to the PASS input")

    warning_lines = [line for line in Path(args.warnings).read_text(errors="replace").splitlines() if line.strip()]
    skipped_count = sum(line.startswith("Line ") for line in Path(args.skipped).read_text(errors="replace").splitlines())
    contamination = list(csv.DictReader(Path(args.contamination).read_text().splitlines(), delimiter="\t"))
    summary = {"schema_version": 1, "counts": {**counts, "vep_records": vep_records},
               "samples": samples, "variant_types": dict(types), "highest_impact": dict(impacts),
               "genes": dict(sorted(genes.items())), "consequences": dict(sorted(consequences.items())),
               "vep_header": vep_version, "contamination": contamination,
               "warnings": {"line_count": len(warning_lines), "skipped_log_records": skipped_count,
                            "examples": warning_lines[:20]},
               "count_definition": "VCF records; multi-allelic records count once. Gene and consequence counts count each term once per record."}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    for name, counter in (("gene", genes), ("consequence", consequences)):
        write_table(output / f"{name}_counts.tsv", [name, "records"],
                    [{name: label, "records": count} for label, count in sorted(counter.items(), key=lambda item: (-item[1], item[0]))])
    custom = {"id": "somvar_variant_counts", "section_name": "Somvar variant counts",
              "description": summary["count_definition"], "plot_type": "table",
              "pconfig": {"id": "somvar_variant_counts_table", "title": "Somvar variant records"},
              "data": {"Matched tumor / normal": summary["counts"]}}
    (output / "somvar_counts_mqc.json").write_text(json.dumps(custom, indent=2) + "\n")
    print(json.dumps(summary["counts"]))


def escape(value):
    return html.escape(str(value), quote=True)


def table(columns, rows, table_id=None):
    identifier = f' id="{escape(table_id)}"' if table_id else ""
    head = "".join(f"<th>{escape(column)}</th>" for column in columns)
    body = "".join("<tr>" + "".join(f"<td>{escape(row.get(column, ''))}</td>" for column in columns) + "</tr>" for row in rows)
    return f'<div class="scroll"><table{identifier}><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def duplicate_metrics(path):
    lines = path.read_text(errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.startswith("LIBRARY\t") and index + 1 < len(lines):
            return dict(zip(line.split("\t"), lines[index + 1].split("\t")))
    return {}


def percent(value):
    try:
        return f"{float(value) * 100:.2f}%"
    except (ValueError, TypeError):
        return "Unavailable"


def build_report(args):
    output = Path(args.outdir)
    if output.exists():
        raise ValueError(f"Report output already exists: {output}")
    output.mkdir(parents=True)
    summary_dir = Path(args.summary_dir)
    shutil.copytree(summary_dir, output / "summary")
    for name, source in (("multiqc_report.html", args.multiqc_report),
                         ("vep_summary.html", args.vep_summary),
                         ("vep.warnings.txt", args.warnings),
                         ("vep.skipped_variants.txt", args.skipped)):
        (output / "qc").mkdir(exist_ok=True)
        shutil.copyfile(source, output / "qc" / name)
    shutil.copytree(args.multiqc_data, output / "qc" / "multiqc_data")
    (output / "variants").mkdir()
    for source in [args.pass_vcf, args.pass_index, args.vep_vcf, args.vep_index]:
        shutil.copyfile(source, output / "variants" / Path(source).name)
    (output / "provenance").mkdir()
    manifest = json.loads(Path(args.manifest).read_text())
    manifest["generated_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    manifest["reporting_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest["multiqc_version"] = Path(args.multiqc_version).read_text().strip()
    (output / "provenance" / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    summary = json.loads((summary_dir / "summary.json").read_text())
    counts = summary["counts"]
    qc_rows = []
    for source in sorted(Path(args.fastp_dir).glob("*_fastp.json")):
        sample = source.name.removesuffix("_fastp.json")
        metrics = json.loads(source.read_text())["summary"]
        before, after = metrics["before_filtering"], metrics["after_filtering"]
        dup = duplicate_metrics(Path(args.duplicate_dir) / f"{sample}.duplicate_metrics.txt")
        qc_rows.append({"Sample": sample, "Raw reads": f"{before['total_reads']:,}",
                        "Clean reads": f"{after['total_reads']:,}",
                        "Clean Q30 bases": percent(after.get("q30_rate")),
                        "Duplicate reads": percent(dup.get("PERCENT_DUPLICATION"))})
    with (summary_dir / "variants.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        variant_rows = []
        for index, row in enumerate(reader):
            if index == 1000:
                break
            variant_rows.append(row)
    count_rows = lambda key_name, label: [{label: name, "Records": count} for name, count in
        sorted(summary[key_name].items(), key=lambda item: (-item[1], item[0]))[:20]]
    cards = "".join(f'<article><strong>{counts[key_name]:,}</strong><span>{label}</span></article>' for key_name, label in
                    (("pass_records", "PASS records"), ("annotated_records", "Annotated"),
                     ("not_in_vep_output_records", "Absent from VEP output"),
                     ("unannotated_output_records", "Output without consequences")))
    warning_text = (f"{counts['not_in_vep_output_records']} PASS records are absent from the VEP output. "
                    "They remain in the PASS VCF and variant table. Review the skipped variants and warnings below."
                    if counts["not_in_vep_output_records"] else "All PASS records are represented in the VEP output.")
    download_names = ["variants.tsv", "variants.csv", "consequences.tsv", "summary.json", "gene_counts.tsv", "consequence_counts.tsv"]
    downloads = "".join(f'<a href="summary/{name}">{name}</a>' for name in download_names)
    displayed = f"Showing {len(variant_rows):,} of {counts['pass_records']:,} records. Download the full TSV or CSV for all records."
    variant_columns = ["chromosome", "position", "ref", "alt", "annotation_status", "genes",
                       "consequences", "tumor_AF", "tumor_DP", "normal_AF", "normal_DP"]
    settings_rows = [{"Setting": name, "Value": json.dumps(value) if isinstance(value, (list, dict)) else value}
                     for name, value in manifest.items()]
    warning_examples = "".join(f"<li>{escape(line)}</li>" for line in summary["warnings"]["examples"])
    contamination_rows = [{"Sample": row.get("sample", ""), "Estimated contamination": percent(row.get("contamination")),
                           "Error (fraction)": row.get("error", "")} for row in summary["contamination"]]
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Somvar analysis report</title><style>
:root{{color-scheme:light;--ink:#172b3a;--muted:#526675;--line:#d9e3e9;--accent:#176b78}}
*{{box-sizing:border-box}}body{{margin:0;background:#f3f6f8;color:var(--ink);font:15px/1.55 system-ui,sans-serif}}
main{{max-width:1280px;margin:auto;padding:32px 24px}}header{{padding:28px;background:#123a4b;color:white;border-radius:12px}}
h1{{margin:0;font-size:32px}}header p{{margin:8px 0 0}}h2{{margin-top:0;font-size:21px}}a{{color:var(--accent)}}
section{{margin-top:22px;padding:24px;background:white;border:1px solid var(--line);border-radius:10px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin-top:22px}}
article{{background:white;padding:22px;border:1px solid var(--line);border-radius:10px}}article strong{{display:block;font-size:32px}}article span,p.note{{color:var(--muted)}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:22px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}th{{background:#eef4f6;white-space:nowrap}}td{{overflow-wrap:anywhere;min-width:85px}}
.downloads{{display:flex;gap:12px;flex-wrap:wrap}}.downloads a{{padding:8px 12px;background:#eef4f6;border-radius:6px}}
input{{padding:10px;width:100%;max-width:450px;border:1px solid #9cb1bf;border-radius:6px;margin:8px 0 16px}}
.notice{{border-left:4px solid #b7791f;padding:10px 16px;background:#fff7e8}}li{{overflow-wrap:anywhere}}
@media(max-width:800px){{.cards{{grid-template-columns:1fr 1fr}}.grid{{grid-template-columns:1fr}}main{{padding:16px}}}}
@media print{{body{{background:white}}main{{padding:0}}header{{color:black;background:white}}input{{display:none}}section,article{{break-inside:avoid}}}}
</style></head><body><main>
<header><h1>Somvar analysis report</h1><p>Matched tumor and normal · GRCh38 · {escape(manifest.get('run_name', ''))}</p>
<p>Generated {escape(manifest['generated_at'])}</p></header>
<div class="cards">{cards}</div>
<section><h2>Annotation completeness</h2><p class="notice">{escape(warning_text)}</p>
<p class="note">{escape(summary['count_definition'])} Allele frequencies, depths and genotypes are copied from the VCF; comma-separated values retain the ALT allele order.</p>
<div class="downloads">{downloads}</div></section>
<section><h2>Read quality and duplicates</h2>{table(['Sample', 'Raw reads', 'Clean reads', 'Clean Q30 bases', 'Duplicate reads'], qc_rows)}
<p><a href="qc/multiqc_report.html">Open MultiQC report</a> · <a href="qc/vep_summary.html">Open VEP summary</a></p>
<p class="note">Raw and clean FastQC are separate sections in MultiQC. Read counts include both mates.</p></section>
<div class="grid"><section><h2>Consequences</h2>{table(['Consequence', 'Records'], count_rows('consequences', 'Consequence'))}</section>
<section><h2>Genes</h2>{table(['Gene', 'Records'], count_rows('genes', 'Gene'))}<p class="note">Up to 20 terms shown per table; full counts are available as TSV.</p></section></div>
<section><h2>Tumor contamination</h2>{table(['Sample', 'Estimated contamination', 'Error (fraction)'], contamination_rows)}</section>
<section><h2>PASS variant table</h2><p class="note">{displayed}</p>
<label for="search">Search visible records by gene, consequence, or position</label><br><input id="search" type="search" placeholder="Search variants">
{table(variant_columns, variant_rows, 'variants')}</section>
<section><h2>Warnings and skipped variants</h2><p>{summary['warnings']['line_count']} warning lines; {summary['warnings']['skipped_log_records']} entries in the VEP skipped log.</p>
<p><a href="qc/vep.warnings.txt">All VEP warnings</a> · <a href="qc/vep.skipped_variants.txt">Skipped variants</a></p>
<details><summary>First {len(summary['warnings']['examples'])} warning lines</summary><ul>{warning_examples}</ul></details></section>
<section><h2>VCFs and provenance</h2><p><a href="variants/{escape(Path(args.pass_vcf).name)}">All PASS variants</a> ·
<a href="variants/{escape(Path(args.vep_vcf).name)}">VEP annotated variants</a> · <a href="provenance/run_manifest.json">Run manifest</a></p>
<details><summary>Run settings</summary>{table(['Setting', 'Value'], settings_rows)}</details>
<details><summary>VEP version and annotation sources</summary><p>{escape(summary['vep_header'] or 'No variants available for annotation')}</p></details></section>
</main><script>
document.getElementById('search').addEventListener('input', function() {{
 const query = this.value.toLowerCase();
 document.querySelectorAll('#variants tbody tr').forEach(row => {{row.hidden = !row.textContent.toLowerCase().includes(query);}});
}});
</script></body></html>'''
    (output / "index.html").write_text(document)
    print(f"Report saved to {output / 'index.html'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    summary = commands.add_parser("summary")
    report = commands.add_parser("report")
    for command in (summary, report):
        for name in ("pass-vcf", "vep-vcf", "warnings", "skipped", "outdir"):
            command.add_argument("--" + name, required=True)
    summary.add_argument("--contamination", required=True)
    for name in ("pass-index", "vep-index", "summary-dir", "multiqc-report", "multiqc-data", "multiqc-version",
                 "vep-summary", "manifest", "fastp-dir", "duplicate-dir"):
        report.add_argument("--" + name, required=True)
    args = parser.parse_args()
    try:
        (summarize if args.command == "summary" else build_report)(args)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"Somvar reporting failed: {error}\n")


if __name__ == "__main__":
    main()
