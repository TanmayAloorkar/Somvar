#!/usr/bin/env python3
"""Run Somvar with a compact live progress display.

Provide paired normal and tumor FASTQs for each run. Defaults: the project
reference, the project results folder, four threads per task, and -resume.
Other Nextflow options can follow the wrapper options.
"""

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


PROJECT = Path(__file__).resolve().parent
DEFAULT_REFERENCE = PROJECT / "references" / "hg38_broad_v0" / "Homo_sapiens_assembly38.fasta"
DEFAULT_OUTDIR = PROJECT / "results"
DEFAULT_THREADS = 4
DEFAULT_COMMON_SITES = PROJECT / "references" / "hg38_broad_v0" / "mutect2" / "small_exac_common_3.hg38.vcf.gz"
DEFAULT_VEP_CACHE = Path.home() / "references" / "vep"
DEFAULT_VEP_CACHE_VERSION = 116
SINGLE_TASK_STAGES = {"MUTECT2", "CALCULATE_CONTAMINATION", "FILTER_MUTECT_CALLS", "PASS_VARIANTS", "VEP", "VARIANT_SUMMARY", "MULTIQC", "FINAL_REPORT"}
MIN_SYSTEM_RAM_BYTES = 16_000_000_000
STAGES = ("FASTQC", "FASTP", "FASTQC_CLEAN", "BWA_MEM", "BAM_PROCESSING")
BQSR_NOTICE = (
    "Note: BQSR runs two GATK passes per sample. Its bar advances when a sample "
    "finishes, so it may appear paused while GATK is processing reads. "
    "Live read counts are in the BQSR task's .command.err log."
)
MUTECT2_NOTICE = (
    "Note: Mutect2's latest GATK ProgressMeter entry appears below the bar while it runs."
)
SAMPLES = ("normal", "tumor")
TASK_EVENT = re.compile(
    r"\[([0-9a-f]{2}/[0-9a-f]+)\] (Submitted|Cached) process > "
    r"(GET_PILEUP_SUMMARIES|CALCULATE_CONTAMINATION|FILTER_MUTECT_CALLS|PASS_VARIANTS|VARIANT_SUMMARY|FINAL_REPORT|MULTIQC|FASTQC_CLEAN|BAM_PROCESSING|BWA_MEM|FASTQC|FASTP|BQSR|MUTECT2|VEP)"
    r"(?: \(([^)]+)\))?"
)
BWA_BATCH = re.compile(r"Processed (\d+) reads in [\d.]+ CPU sec, [\d.]+ real sec")
MUTECT2_PROGRESS = re.compile(
    r"^\d{2}:\d{2}:\d{2}\.\d+\s+INFO\s+ProgressMeter\s+-\s+"
    r"\S+:\d+\s+[\d.]+\s+\d+\s+[\d.]+$"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalR1", required=True, help="normal read 1 FASTQ")
    parser.add_argument("--normalR2", required=True, help="normal read 2 FASTQ")
    parser.add_argument("--tumorR1", required=True, help="tumor read 1 FASTQ")
    parser.add_argument("--tumorR2", required=True, help="tumor read 2 FASTQ")
    parser.add_argument("--ref", default=str(DEFAULT_REFERENCE),
                        help="reference FASTA (default: references/hg38_broad_v0/Homo_sapiens_assembly38.fasta)")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR),
                        help="result directory (default: project results/)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                        help="threads per sample task (default: 4)")
    parser.add_argument("--known_sites", help="comma-separated indexed VCFs to enable optional BQSR")
    parser.add_argument("--germline_resource", help="reference-matched population AF VCF for Mutect2")
    parser.add_argument("--pon", help="reference-matched panel-of-normals VCF for Mutect2")
    parser.add_argument("--intervals", help="optional exome target BED or interval list for Mutect2")
    parser.add_argument("--common_sites", default=str(DEFAULT_COMMON_SITES),
                        help="indexed common biallelic SNP VCF with INFO/AF for contamination estimation (default: bundled hg38 ExAC sites)")
    parser.add_argument("--no-resume", action="store_true", help="start without Nextflow -resume")
    parser.add_argument("--vep_cache", default=str(DEFAULT_VEP_CACHE),
                        help="VEP cache root containing homo_sapiens/ (default: ~/references/vep)")
    parser.add_argument("--vep_cache_version", type=int, default=DEFAULT_VEP_CACHE_VERSION,
                        help="Ensembl cache release matching your VEP software (default: 116)")
    parser.add_argument("--skip_vep", "--skip-vep", action="store_true",
                        help="stop after PASS variants while VEP is being installed")
    args, extra = parser.parse_known_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.vep_cache_version < 1:
        parser.error("--vep_cache_version must be at least 1")
    removed_resource_options = ("--ram-gb", "--ram", "--max-cpus")
    if any(option == token or token.startswith(option + "=")
           for token in extra for option in removed_resource_options):
        parser.error("RAM and total CPU budgets are detected by Nextflow; use --threads for each task")
    return args, extra


def read_text(path):
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


def work_directory(short_hash):
    prefix, suffix = short_hash.split("/", 1)
    matches = list((PROJECT / "work" / prefix).glob(suffix + "*"))
    return matches[0] if matches else None


def bwa_fraction(sample, work, fastp_workdirs):
    report = next(
        (directory / f"{sample}_fastp.json" for directory in reversed(fastp_workdirs)
         if (directory / f"{sample}_fastp.json").is_file()),
        None,
    )
    if report is None:
        return None
    try:
        total = json.loads(report.read_text())["summary"]["after_filtering"]["total_reads"]
    except (OSError, ValueError, KeyError, TypeError):
        return None
    if not isinstance(total, int) or total <= 0 or work is None:
        return None
    processed = sum(int(n) for n in BWA_BATCH.findall(read_text(work / ".command.err")))
    return min(1.0, processed / total)


def task_states(log_path):
    tasks = {}
    events = list(TASK_EVENT.finditer(read_text(log_path)))
    fastp_workdirs = [directory for match in events if match.group(3) == "FASTP"
                      if (directory := work_directory(match.group(1))) is not None]
    for match in events:
        short_hash, event, stage, tag = match.groups()
        key = (stage, tag or short_hash)
        work = work_directory(short_hash)
        if event == "Cached":
            status = "done"
        else:
            exit_code = read_text(work / ".exitcode").strip() if work else ""
            status = "done" if exit_code == "0" else ("failed" if exit_code else "running")
        fraction = 0.0
        if stage == "BWA_MEM" and status == "running" and tag in SAMPLES:
            fraction = bwa_fraction(tag, work, fastp_workdirs) or 0.0
        tasks[key] = (status, fraction)
    return tasks


def latest_mutect2_progress(tasks):
    for (stage, identifier), (status, _) in tasks.items():
        if stage != "MUTECT2" or status != "running":
            continue
        work = work_directory(identifier)
        if work is None:
            continue
        for line in reversed(read_text(work / ".command.err").splitlines()):
            if MUTECT2_PROGRESS.fullmatch(line.strip()):
                return line.strip()
    return None


def progress(tasks, stages=STAGES, finished=False):
    if finished:
        return "Complete", 1.0, 1.0
    overall_units = 0.0
    active_stage = None
    active_fraction = 0.0
    for stage in stages:
        expected = 1 if stage in SINGLE_TASK_STAGES else 2
        states = [value for (name, _), value in tasks.items() if name == stage]
        units = min(expected, sum(1.0 if status == "done" else fraction
                                  for status, fraction in states))
        overall_units += units
        if active_stage is None and units < expected:
            active_stage = stage
            active_fraction = units / expected
    if active_stage is None:
        active_stage = stages[-1]
        active_fraction = 1.0
    total_units = sum(1 if stage in SINGLE_TASK_STAGES else 2 for stage in stages)
    return active_stage, active_fraction, min(1.0, overall_units / total_units)


def bar(fraction, width=24):
    filled = min(width, max(0, round(fraction * width)))
    return "#" * filled + "-" * (width - filled)


def stamped(message):
    return "\n".join(f"[{datetime.now():%H:%M:%S}] {line}"
                     for line in str(message).splitlines())


def print_message(message):
    print(stamped(message), flush=True)


def close_progress_line():
    if sys.stdout.isatty() and show.open_line:
        if show.detail_open:
            sys.stdout.write("\r\033[2K\033[1A\r")
        sys.stdout.write("\n")
        sys.stdout.flush()
        show.open_line = False
        show.detail_open = False


def fit_detail(line):
    columns = shutil.get_terminal_size((80, 24)).columns
    if len(line) < columns:
        return line
    compact = " ".join(line.split())
    return compact if len(compact) < columns else compact[:max(0, columns - 4)] + "..."


def show(line, stage, detail=None, final=False):
    changed = stage != show.stage
    if changed:
        close_progress_line()
        show.stage = stage
        show.stamp = f"[{datetime.now():%H:%M:%S}] "
    display = show.stamp + line
    if sys.stdout.isatty():
        if show.detail_open:
            prefix = "\r\033[2K\033[1A\r\033[2K"
        else:
            prefix = "\r\033[2K" if show.open_line else ""
        output = prefix + display
        if detail and not final:
            output += "\n" + fit_detail(detail)
        if final:
            output += "\n"
        sys.stdout.write(output)
        sys.stdout.flush()
        show.open_line = not final
        show.detail_open = bool(detail) and not final
    elif changed or final:
        print(display, flush=True)


show.stage = None
show.stamp = ""
show.open_line = False
show.detail_open = False


def progress_line(status, stage, stage_fraction, overall):
    head = f"Somvar - {status} {stage} - ["
    tail = f"] - {stage_fraction:.0%} | Overall {overall:.0%}"
    columns = shutil.get_terminal_size((80, 24)).columns
    width = max(4, min(24, columns - len(head) - len(tail) - 12))
    return f"{head}{bar(stage_fraction, width)}{tail}"


def stop_child(child):
    if child.poll() is not None:
        return
    os.killpg(child.pid, signal.SIGINT)
    try:
        child.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()


def main():
    args, extra = parse_args()
    somatic_options = (args.germline_resource, args.pon)
    if (any(somatic_options) or args.intervals) and not all(somatic_options):
        raise SystemExit(stamped(
            "Mutect2 requires --germline_resource and --pon together; --intervals is optional."
        ))
    stages = STAGES + (("BQSR",) if args.known_sites else ())
    if all(somatic_options):
        stages += ("MUTECT2", "GET_PILEUP_SUMMARIES", "CALCULATE_CONTAMINATION", "FILTER_MUTECT_CALLS", "PASS_VARIANTS")
        if not args.skip_vep:
            stages += ("VEP", "VARIANT_SUMMARY", "MULTIQC", "FINAL_REPORT")
    system_ram_bytes = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    if system_ram_bytes < MIN_SYSTEM_RAM_BYTES:
        raise SystemExit(stamped("Somvar requires at least 16 GB of installed RAM."))
    nextflow = shutil.which("nextflow")
    if nextflow is None:
        raise SystemExit(stamped("Nextflow is not on PATH."))
    pipeline = PROJECT / "scripts" / "main.nf"
    if not pipeline.is_file():
        raise SystemExit(stamped(f"Cannot find {pipeline}"))
    reads = {}
    for name in ("normalR1", "normalR2", "tumorR1", "tumorR2"):
        path = Path(getattr(args, name)).expanduser().resolve()
        if not path.is_file():
            raise SystemExit(stamped(f"Missing --{name}: {path}"))
        reads[name] = str(path)

    reference = Path(args.ref).expanduser().resolve()
    if not reference.is_file():
        raise SystemExit(stamped(f"Missing reference FASTA: {reference}. Use --ref to choose one."))
    missing_indexes = [str(reference) + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
                       if not Path(str(reference) + suffix).is_file()]
    if missing_indexes:
        raise SystemExit(stamped("Missing BWA reference indexes:\n" + "\n".join(missing_indexes)))

    somatic_resources = {}
    if all(somatic_options):
        for option, value in (("germline_resource", args.germline_resource),
                              ("pon", args.pon), ("intervals", args.intervals),
                              ("common_sites", args.common_sites)):
            if value is None:
                continue
            path = Path(value).expanduser().resolve()
            if not path.is_file():
                if option == "common_sites":
                    raise SystemExit(stamped(
                        f"Missing --common_sites: {path}. Download the common-SNP VCF and index "
                        f"using {PROJECT / 'README.md'}, or provide --common_sites."
                    ))
                raise SystemExit(stamped(f"Missing --{option}: {path}"))
            if option != "intervals" and not any(
                Path(str(path) + suffix).is_file() for suffix in (".tbi", ".idx")
            ):
                raise SystemExit(stamped(f"Missing .tbi or .idx index for --{option}: {path}"))
            somatic_resources[option] = str(path)

    vep_cache = Path(args.vep_cache).expanduser().resolve()
    if all(somatic_options):
        if shutil.which("bcftools") is None:
            raise SystemExit(stamped("bcftools is required for PASS variant extraction and indexing."))
        if not args.skip_vep:
            for tool in ("vep", "bgzip", "python3", "multiqc"):
                if shutil.which(tool) is None:
                    raise SystemExit(stamped(
                        f"{tool} is not on PATH. See {PROJECT / 'README.md'} "
                        "for setup commands, or use --skip_vep to stop after PASS variants."
                    ))
            if not (PROJECT / "scripts" / "reporting.py").is_file():
                raise SystemExit(stamped("Missing scripts/reporting.py for summaries and final reporting."))
            cache_info = vep_cache / "homo_sapiens" / f"{args.vep_cache_version}_GRCh38" / "info.txt"
            if not cache_info.is_file():
                raise SystemExit(stamped(
                    f"Missing GRCh38 VEP cache: {cache_info}. Use --vep_cache and --vep_cache_version "
                    "to choose an installed cache, or follow README.md (VEP cache)."
                ))

    outdir = Path(args.outdir).expanduser().resolve()

    log_dir = PROJECT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d-%H%M%S") + f"-{os.getpid()}"
    run_dir = log_dir / run_id
    run_dir.mkdir()
    console_log = run_dir / "console.log"
    nextflow_log = run_dir / ".nextflow.log"
    command = [nextflow, "-log", str(nextflow_log), "run", "scripts/main.nf", "-ansi-log", "false"]
    if not args.no_resume and "-resume" not in extra:
        command.append("-resume")
    for name, value in reads.items():
        command.extend((f"--{name}", value))
    command.extend(("--threads", str(args.threads)))
    command.extend(("--ref", str(reference), "--outdir", str(outdir)))
    if args.known_sites:
        command.extend(("--known_sites", args.known_sites))
    for name, value in somatic_resources.items():
        command.extend((f"--{name}", value))
    if all(somatic_options):
        command.extend(("--skip_vep", str(args.skip_vep).lower()))
        if not args.skip_vep:
            command.extend(("--vep_cache", str(vep_cache), "--vep_cache_version", str(args.vep_cache_version)))
    command.extend(extra)

    (run_dir / "settings.json").write_text(json.dumps({
        "reads": reads,
        "reference": str(reference),
        "outdir": str(outdir),
        "threads_per_task": args.threads,
        "known_sites": args.known_sites,
        "somatic_resources": somatic_resources,
        "annotation": {"enabled": all(somatic_options) and not args.skip_vep,
                       "vep_cache": str(vep_cache), "vep_cache_version": args.vep_cache_version},
        "reporting_enabled": all(somatic_options) and not args.skip_vep,
        "resume": not args.no_resume,
        "command": command,
    }, indent=2) + "\n")

    with console_log.open("w") as output:
        child = subprocess.Popen(
            command, cwd=PROJECT, stdout=output, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        bqsr_notice_shown = False
        mutect2_notice_shown = False
        try:
            while child.poll() is None:
                tasks = task_states(nextflow_log)
                stage, stage_fraction, overall = progress(tasks, stages=stages)
                if stage == "BQSR" and not bqsr_notice_shown:
                    close_progress_line()
                    print_message(BQSR_NOTICE)
                    bqsr_notice_shown = True
                if stage == "MUTECT2" and not mutect2_notice_shown:
                    close_progress_line()
                    print_message(MUTECT2_NOTICE)
                    mutect2_notice_shown = True
                status = "Running" if any(stage == key[0] for key in tasks) else "Waiting for"
                detail = latest_mutect2_progress(tasks) if stage == "MUTECT2" else None
                show(progress_line(status, stage, stage_fraction, overall), stage, detail=detail)
                time.sleep(1)
        except KeyboardInterrupt:
            stop_child(child)
            show("Somvar - Stopped", "Stopped", final=True)
            print_message(f"Logs: {run_dir}")
            return 130

    if child.returncode == 0:
        show(f"Somvar - Complete - [{bar(1.0)}] - 100%", "Complete", final=True)
        print_message(f"Logs: {run_dir}")
        return 0

    show("Somvar - Failed", "Failed", final=True)
    print_message(f"Logs: {run_dir}")
    for line in read_text(console_log).splitlines()[-12:]:
        print_message(line)
    return child.returncode or 1


if __name__ == "__main__":
    sys.exit(main())
