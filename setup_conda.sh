#!/usr/bin/env bash
# Create a Somvar conda environment after interactive confirmation.
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    printf 'Run this installer with: bash setup_conda.sh\n' >&2
    return 1
fi

environment_name=somvar
case "$#" in
    0) ;;
    1)
        if [[ "$1" == --help || "$1" == -h ]]; then
            printf 'Usage: bash setup_conda.sh [--name ENVIRONMENT]\n'
            exit 0
        fi
        printf 'Usage: bash setup_conda.sh [--name ENVIRONMENT]\n' >&2
        exit 2
        ;;
    2)
        [[ "$1" == --name ]] || { printf 'Expected --name ENVIRONMENT\n' >&2; exit 2; }
        environment_name=$2
        ;;
    *) printf 'Usage: bash setup_conda.sh [--name ENVIRONMENT]\n' >&2; exit 2 ;;
esac
if [[ ! "$environment_name" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ || "${environment_name,,}" == base ]]; then
    printf 'Choose a name using letters, digits, hyphens or underscores; base is reserved.\n' >&2
    exit 2
fi

printf 'Somvar conda setup: create "%s" and install all workflow software.\n' "$environment_name"
read -r -p 'Proceed with building the conda environment? [y/N]: ' answer || answer=''
case "${answer,,}" in
    y|yes) ;;
    *) printf 'Setup cancelled.\n'; exit 0 ;;
esac

creation_started=false
on_error() {
    local status=$?
    printf '\nSetup failed. Review the error above.\n' >&2
    if [[ "$creation_started" == true ]]; then
        printf 'A partial environment may remain. Use --name to choose another name for a fresh attempt.\n' >&2
    fi
    exit "$status"
}
trap on_error ERR

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    printf 'This environment specification supports Linux x86_64, including compatible WSL2 installations.\n' >&2
    exit 1
fi
conda_exe=${CONDA_EXE:-$(command -v conda || true)}
if [[ -z "$conda_exe" || ! -x "$conda_exe" ]]; then
    printf 'Install Miniforge or Miniconda and make conda available in this terminal. See README.md (Conda setup).\n' >&2
    exit 1
fi
project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
specification="$project_dir/environment.yml"
if [[ ! -f "$specification" ]]; then
    printf 'Missing environment definition: %s\n' "$specification" >&2
    exit 1
fi

# Use conda's own Python to read its JSON; no extra system Python is needed.
conda_base=$("$conda_exe" info --base)
conda_info=$("$conda_exe" info --json)
existing=$("$conda_base/bin/python" -c '
import json, pathlib, sys
info = json.load(sys.stdin)
name = sys.argv[1]
candidates = [pathlib.Path(p) for p in info.get("envs", []) if pathlib.Path(p).name == name]
candidates += [pathlib.Path(p) / name for p in info["envs_dirs"] if (pathlib.Path(p) / name).exists()]
print(str(candidates[0]) if candidates else "")
' "$environment_name" <<< "$conda_info")
if [[ -n "$existing" ]]; then
    printf 'Environment path already exists: %s\nChoose another name with --name or activate your existing environment.\n' "$existing" >&2
    exit 1
fi

printf 'Creating the environment. Conda will display installation progress.\n'
creation_started=true
CONDA_CHANNEL_PRIORITY=strict "$conda_exe" env create \
    --name "$environment_name" --file "$specification" --no-default-packages --yes

printf 'Checking the installed environment...\n'
"$conda_exe" run --no-capture-output --name "$environment_name" bash -c '
set -euo pipefail
for tool in python3 nextflow java fastqc fastp bwa samtools gatk bcftools \
    bgzip tabix vep multiqc aria2c curl unzip bash tar gzip md5sum; do
    executable=$(command -v "$tool") || { printf "Missing executable: %s\n" "$tool" >&2; exit 1; }
    case "$executable" in
        "$CONDA_PREFIX"/*) ;;
        *) printf "Executable is outside this environment: %s\n" "$executable" >&2; exit 1 ;;
    esac
done
check_command() {
    local output
    if output=$("$@" 2>&1); then
        printf "OK: %s\n" "$*"
    else
        printf "Check failed: %s\n%s\n" "$*" "$output" >&2
        exit 1
    fi
}
check_command nextflow -version
check_command gatk --version
check_command vep --help
check_command multiqc --version
printf "All required executables are present in the conda environment.\n"
'

lock_file="$project_dir/conda-${environment_name}-linux-64.lock.txt"
"$conda_exe" list --name "$environment_name" --explicit > "$lock_file"
printf '\nSetup complete. Exact package builds saved to %s\n' "$lock_file"
printf 'Run these commands in your terminal:\n  conda activate %s\n  python3 run_somvar.py --help\n' "$environment_name"
printf 'Then prepare the reference files and VEP cache using README.md.\n'
