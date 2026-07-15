#!/usr/bin/env bash
set -u -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

LOG_DIR="logs"
RUN_ID="$(date '+%Y%m%d-%H%M%S')"
LOG_FILE="$LOG_DIR/leb-regeneration-$RUN_ID.log"
LOCK_DIR=".leb-regeneration.lock"
LOCK_WAIT_SECONDS="${LOCK_WAIT_SECONDS:-60}"
LOCK_ACQUIRED=0

TIER="all"
GROUP_SPEC="all"
LEB2_GROUP_SPEC=""
VERIFY=1
VERIFY_SAMPLES=500
LEB2_VERIFY_SAMPLES=200
KEEP_SPK_CACHE=0
DO_LEB1=1
DO_LEB2=1
DO_MERGE=1
INSTALL=0
INSTALL_DIR="${LIBEPHEMERIS_DATA_DIR:-$HOME/.libephemeris}/leb"
DRY_RUN=0
QUIET=0

LEB1_GROUPS=(planets asteroids exotics analytical)
# Keep this build-orchestration list in exactly the same order as
# ``libephemeris.leb_groups.LEB2_GROUPS``.  It is intentionally repeated here
# because this script must be able to parse its arguments before a project
# Python interpreter has been selected.  The focused regression test compares
# the two definitions so future group changes cannot silently drift apart.
#
# There is deliberately no ``uranians`` companion.  That historical group was
# retired together with its unsupported hypothetical-body inventory.  The
# surviving analytical LEB1 group contributes only to ``core`` (nodes) and
# ``apogee`` (apsides) in LEB2.
LEB2_GROUPS=(core asteroids exotics apogee)
TIERS=(base medium extended)
TNO_SPK_IDS=(
  136199  # Eris
  90377   # Sedna
  136108  # Haumea
  136472  # Makemake
  28978   # Ixion
  90482   # Orcus
  50000   # Quaoar
)

FAILED_COMMANDS=()
FAILED_STATUSES=()
GENERATED_LEB1_MAIN=()
GENERATED_LEB2_FILES=()

usage() {
  cat <<'EOF'
Usage:
  ./regenerate-leb.sh [base|medium|extended|all] [options]
  ./regenerate-leb.sh --tier base|medium|extended|all [options]

Regenerate LibEphemeris LEB1 and LEB2 files. This is intended to be the
single safe entrypoint for full tier regeneration and for partial group
updates such as "medium exotics only".

Default:
  ./regenerate-leb.sh all

Default workflow per selected tier:
  1. if exotics are selected, move cached SPKs for the 7 refitted TNOs to a
     timestamped logs/ backup so Horizons re-downloads fresh kernels
  2. generate selected LEB1 group files:
       planets, asteroids, exotics, analytical
  3. merge:
       - all groups: rebuild data/leb/ephemeris_<tier>.leb from all partials
       - subset: replace only selected group bodies inside the existing main LEB1
  4. verify LEB1 unless --no-verify
  5. convert affected LEB2 groups unless --leb1-only
  6. verify LEB2 unless --no-verify

Tier selection:
  --tier TIER             base, medium, extended, all. Default: all.
  base|medium|extended|all may also be passed positionally.

LEB1 group selection:
  --group GROUP           One LEB1 group: planets, asteroids, exotics,
                          analytical, all.
  --groups LIST           Comma-separated LEB1 groups, e.g.
                          --groups planets,exotics.
  --no-merge              Generate partial LEB1 group files only; do not update
                          data/leb/ephemeris_<tier>.leb.

LEB2 selection:
  --leb1-only             Generate/merge/verify LEB1 only.
  --leb2-only             Convert/verify LEB2 only from existing LEB1 files.
  --leb2-groups LIST      Comma-separated LEB2 groups: core, asteroids,
                          exotics, apogee, all. If omitted, the
                          affected LEB2 groups are inferred from selected LEB1
                          groups.

Verification:
  --no-verify             Skip LEB1 and LEB2 verification.
  --verify-samples N      LEB1 verify samples per body. Default: 500.
  --leb2-verify-samples N LEB2 verify samples per body. Default: 200.

SPK cache:
  --keep-spk-cache        Do not move cached SPKs before exotics generation.
                          Use only after auditing ~/.libephemeris/spk.

Install:
  --install               Copy generated main LEB1 and affected LEB2 files to
                          ${LIBEPHEMERIS_DATA_DIR:-~/.libephemeris}/leb with
                          timestamped backups.
  --install-dir DIR       Override install destination.

Other:
  --dry-run               Print actions without running generators or copying.
  -q, --quiet             Pass quiet mode to generators where supported.
  -h, --help              Show this help.

Examples:
  ./regenerate-leb.sh medium --group exotics --install
  ./regenerate-leb.sh base --groups planets,asteroids
  ./regenerate-leb.sh all --group all --no-verify
  ./regenerate-leb.sh medium --leb2-only --leb2-groups exotics --install
  ./regenerate-leb.sh extended --group exotics --leb1-only

Notes:
  - Extended generation can take hours.
  - "planets" is the LEB1 Skyfield planet group. LEB2 "core" is different and
    also contains node/apogee analytical bodies, so LEB2 groups are inferred
    conservatively from the changed LEB1 groups.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

contains_word() {
  local needle="$1"
  shift
  local item
  for item in "$@"; do
    [[ "$item" == "$needle" ]] && return 0
  done
  return 1
}

normalize_csv() {
  local value="$1"
  value="${value// /}"
  value="${value//;/,}"
  printf '%s\n' "$value"
}

parse_args() {
  while (($#)); do
    case "$1" in
      --tier)
        [[ $# -ge 2 ]] || die "--tier requires a value"
        TIER="$2"
        shift 2
        ;;
      --tier=*)
        TIER="${1#--tier=}"
        shift
        ;;
      --group)
        [[ $# -ge 2 ]] || die "--group requires a value"
        GROUP_SPEC="$2"
        shift 2
        ;;
      --group=*)
        GROUP_SPEC="${1#--group=}"
        shift
        ;;
      --groups)
        [[ $# -ge 2 ]] || die "--groups requires a value"
        GROUP_SPEC="$2"
        shift 2
        ;;
      --groups=*)
        GROUP_SPEC="${1#--groups=}"
        shift
        ;;
      --leb2-groups)
        [[ $# -ge 2 ]] || die "--leb2-groups requires a value"
        LEB2_GROUP_SPEC="$2"
        shift 2
        ;;
      --leb2-groups=*)
        LEB2_GROUP_SPEC="${1#--leb2-groups=}"
        shift
        ;;
      --leb1-only)
        DO_LEB1=1
        DO_LEB2=0
        shift
        ;;
      --leb2-only)
        DO_LEB1=0
        DO_LEB2=1
        shift
        ;;
      --no-merge)
        DO_MERGE=0
        shift
        ;;
      --no-verify)
        VERIFY=0
        shift
        ;;
      --verify-samples)
        [[ $# -ge 2 ]] || die "--verify-samples requires a value"
        VERIFY_SAMPLES="$2"
        shift 2
        ;;
      --verify-samples=*)
        VERIFY_SAMPLES="${1#--verify-samples=}"
        shift
        ;;
      --leb2-verify-samples)
        [[ $# -ge 2 ]] || die "--leb2-verify-samples requires a value"
        LEB2_VERIFY_SAMPLES="$2"
        shift 2
        ;;
      --leb2-verify-samples=*)
        LEB2_VERIFY_SAMPLES="${1#--leb2-verify-samples=}"
        shift
        ;;
      --keep-spk-cache)
        KEEP_SPK_CACHE=1
        shift
        ;;
      --install)
        INSTALL=1
        shift
        ;;
      --install-dir)
        [[ $# -ge 2 ]] || die "--install-dir requires a value"
        INSTALL_DIR="$2"
        shift 2
        ;;
      --install-dir=*)
        INSTALL_DIR="${1#--install-dir=}"
        shift
        ;;
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      -q|--quiet)
        QUIET=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      base|medium|extended|all)
        TIER="$1"
        shift
        ;;
      *)
        die "unknown argument: $1"
        ;;
    esac
  done

  case "$TIER" in
    base|medium|extended|all) ;;
    *) die "invalid tier: $TIER" ;;
  esac

  [[ "$VERIFY_SAMPLES" =~ ^[0-9]+$ ]] || die "--verify-samples must be an integer"
  [[ "$LEB2_VERIFY_SAMPLES" =~ ^[0-9]+$ ]] || die "--leb2-verify-samples must be an integer"

  # Validate group specifications here, in the parent shell.  Most execution
  # loops consume the selectors through process substitution; an ``exit``
  # raised while producing that substitution would terminate only its child
  # shell and could otherwise leave the orchestrator running with no groups.
  # Eager validation therefore makes stale/unknown names fail closed before a
  # lock, backup, output file, or generator process is touched.
  selected_leb1_groups >/dev/null
  selected_leb2_groups >/dev/null
}

cleanup() {
  # Only remove the lock this process actually acquired. Interrupting while
  # still waiting for another holder's lock must not delete that holder's lock.
  if ((LOCK_ACQUIRED)); then
    rm -rf "$LOCK_DIR"
  fi
}

run() {
  local statuses
  local cmd_status
  local tee_status

  echo
  echo "================================================================"
  echo "START: $*"
  echo "TIME : $(date '+%Y-%m-%d %H:%M:%S')"
  echo "================================================================"

  if ((DRY_RUN)); then
    echo "DRY RUN: $*"
    echo "DONE : $*"
    return 0
  fi

  "$@" 2>&1 | tee -a "$LOG_FILE"
  statuses=("${PIPESTATUS[@]}")
  cmd_status="${statuses[0]}"
  tee_status="${statuses[1]:-0}"

  echo
  if ((cmd_status == 0)); then
    echo "DONE : $*"
  else
    echo "FAILED: $*"
    echo "EXIT CODE: $cmd_status"
    FAILED_COMMANDS+=("$*")
    FAILED_STATUSES+=("$cmd_status")
  fi

  if ((tee_status != 0)); then
    echo "WARNING: failed to write command output to log file: $LOG_FILE"
  fi

  echo "TIME : $(date '+%Y-%m-%d %H:%M:%S')"
  return "$cmd_status"
}

select_python() {
  if [[ -x ".venv/bin/python" ]]; then
    PYTHON=".venv/bin/python"
  elif command -v python >/dev/null 2>&1; then
    PYTHON="python"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON="python3"
  else
    die "no Python interpreter found"
  fi
  export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
}

selected_tiers() {
  if [[ "$TIER" == "all" ]]; then
    printf '%s\n' "${TIERS[@]}"
  else
    printf '%s\n' "$TIER"
  fi
}

selected_leb1_groups() {
  local spec
  local group
  local candidate
  # macOS ships Bash 3.2, where expanding an empty array under ``set -u`` is
  # an error even after ``local selected=()``.  A private sentinel keeps the
  # array bound while preserving portable ordered de-duplication.
  local selected=(__leph_no_group__)
  local candidates=()

  spec="$(normalize_csv "$GROUP_SPEC")"
  if [[ -z "$spec" || "$spec" == "all" ]]; then
    printf '%s\n' "${LEB1_GROUPS[@]}"
    return 0
  fi

  IFS=',' read -r -a _groups <<< "$spec"
  for group in "${_groups[@]}"; do
    case "$group" in
      planets|asteroids|exotics|analytical) candidates=("$group") ;;
      all)
        candidates=("${LEB1_GROUPS[@]}")
        ;;
      *) die "invalid LEB1 group: $group" ;;
    esac
    for candidate in "${candidates[@]}"; do
      if ! contains_word "$candidate" "${selected[@]}"; then
        selected+=("$candidate")
      fi
    done
  done
  for candidate in "${selected[@]}"; do
    [[ "$candidate" == "__leph_no_group__" ]] || printf '%s\n' "$candidate"
  done
}

selected_leb2_groups() {
  local spec
  local group
  local candidate
  # See ``selected_leb1_groups`` for the Bash-3.2/nounset sentinel rationale.
  local selected=(__leph_no_group__)
  local candidates=()

  spec="$(normalize_csv "$LEB2_GROUP_SPEC")"
  if [[ -n "$spec" ]]; then
    if [[ "$spec" == "all" ]]; then
      printf '%s\n' "${LEB2_GROUPS[@]}"
      return 0
    fi
    IFS=',' read -r -a _groups <<< "$spec"
    for group in "${_groups[@]}"; do
      case "$group" in
        core|asteroids|exotics|apogee) candidates=("$group") ;;
        all) candidates=("${LEB2_GROUPS[@]}") ;;
        *) die "invalid LEB2 group: $group" ;;
      esac
      for candidate in "${candidates[@]}"; do
        if ! contains_word "$candidate" "${selected[@]}"; then
          selected+=("$candidate")
        fi
      done
    done
    for candidate in "${selected[@]}"; do
      [[ "$candidate" == "__leph_no_group__" ]] || printf '%s\n' "$candidate"
    done
    return 0
  fi

  affected_leb2_groups
}

affected_leb2_groups() {
  local group

  while IFS= read -r group; do
    case "$group" in
      planets)
        printf '%s\n' core
        ;;
      asteroids)
        printf '%s\n' asteroids
        ;;
      exotics)
        printf '%s\n' exotics
        ;;
      analytical)
        # The analytical LEB1 partition contains the public node/apsis models.
        # In LEB2 those bodies are distributed between core and apogee only;
        # the former hypothetical ``uranians`` output no longer exists.
        printf '%s\n' core apogee
        ;;
    esac
  done < <(selected_leb1_groups) | awk '!seen[$0]++'
}

all_leb1_groups_selected() {
  local count=0
  local group
  while IFS= read -r group; do
    count=$((count + 1))
  done < <(selected_leb1_groups)
  [[ "$count" -eq "${#LEB1_GROUPS[@]}" ]]
}

exotics_selected() {
  local group
  while IFS= read -r group; do
    [[ "$group" == "exotics" ]] && return 0
  done < <(selected_leb1_groups)
  return 1
}

main_leb1_path() {
  local tier="$1"
  printf 'data/leb/ephemeris_%s.leb\n' "$tier"
}

partial_leb1_path() {
  local tier="$1"
  local group="$2"
  printf 'data/leb/ephemeris_%s_%s.leb\n' "$tier" "$group"
}

leb2_path() {
  local tier="$1"
  local group="$2"
  printf 'data/leb2/%s_%s.leb2\n' "$tier" "$group"
}

backup_path_if_exists() {
  local path="$1"
  local backup_dir="$LOG_DIR/leb-backup-$RUN_ID"

  [[ -e "$path" ]] || return 0
  mkdir -p "$backup_dir"
  if ((DRY_RUN)); then
    echo "DRY RUN: cp -p $path $backup_dir/$(basename "$path")"
    return 0
  fi
  cp -p "$path" "$backup_dir/$(basename "$path")"
  echo "Backed up: $path -> $backup_dir/$(basename "$path")"
}

unique_spk_dirs() {
  local data_dir="${LIBEPHEMERIS_DATA_DIR:-$HOME/.libephemeris}"
  local dirs=("$data_dir/spk")

  if [[ -n "${LIBEPHEMERIS_SPK_DIR:-}" ]]; then
    dirs+=("$LIBEPHEMERIS_SPK_DIR")
  fi

  printf '%s\n' "${dirs[@]}" | awk '!seen[$0]++'
}

backup_problem_tno_spks() {
  local backup_dir="$LOG_DIR/spk-backup-$RUN_ID"
  local moved=0
  local spk_dir
  local id
  local path

  if ((KEEP_SPK_CACHE)); then
    echo "SPK cache refresh skipped (--keep-spk-cache)."
    return 0
  fi

  if ! ((DO_LEB1)) || ! exotics_selected; then
    echo "SPK cache refresh skipped (exotics LEB1 generation not selected)."
    return 0
  fi

  shopt -s nullglob
  while IFS= read -r spk_dir; do
    [[ -d "$spk_dir" ]] || continue
    for id in "${TNO_SPK_IDS[@]}"; do
      for path in "$spk_dir"/"${id}"_*.bsp; do
        [[ -e "$path" ]] || continue
        mkdir -p "$backup_dir"
        if ((DRY_RUN)); then
          echo "DRY RUN: mv $path $backup_dir/$(basename "$path")"
        else
          mv "$path" "$backup_dir/$(basename "$path")"
        fi
        moved=$((moved + 1))
      done
    done
  done < <(unique_spk_dirs)
  shopt -u nullglob

  if ((moved == 0)); then
    echo "No cached SPKs found for the 7 refitted TNOs."
  else
    echo "Moved $moved SPK file(s) to $backup_dir."
    echo "Generation will re-download fresh JPL Horizons kernels as needed."
  fi
}

generate_leb1_groups_for_tier() {
  local tier="$1"
  local group
  local args

  echo
  echo "### LEB1 group generation: $tier"

  mkdir -p data/leb
  while IFS= read -r group; do
    args=(scripts/generate_leb.py --tier "$tier" --group "$group" --output "$(partial_leb1_path "$tier" "$group")")
    ((QUIET)) && args+=(--quiet)
    run "$PYTHON" "${args[@]}" || return 1
  done < <(selected_leb1_groups)
}

replace_leb1_from_partials() {
  local tier="$1"
  local output
  local partials=()
  local group
  local statuses
  local cmd_status
  local tee_status

  output="$(main_leb1_path "$tier")"
  if ! ((DRY_RUN)); then
    [[ -f "$output" ]] || die "cannot update selected groups: missing existing main LEB1: $output"
  fi

  while IFS= read -r group; do
    partials+=("$(partial_leb1_path "$tier" "$group")")
  done < <(selected_leb1_groups)

  if ! ((DRY_RUN)); then
    for group in "${partials[@]}"; do
      [[ -f "$group" ]] || die "missing generated partial: $group"
    done
  fi

  backup_path_if_exists "$output"

  echo
  echo "================================================================"
  echo "START: replace selected LEB1 bodies in $output"
  echo "TIME : $(date '+%Y-%m-%d %H:%M:%S')"
  echo "================================================================"

  if ((DRY_RUN)); then
    echo "DRY RUN: replace bodies from: ${partials[*]}"
    echo "DONE : replace selected LEB1 bodies in $output"
    return 0
  fi

  "$PYTHON" - "$output" "$output" "${partials[@]}" <<'PY' 2>&1 | tee -a "$LOG_FILE"
from __future__ import annotations

import mmap
import os
import sys
import tempfile
import time
from pathlib import Path

from libephemeris.leb_format import (
    BODY_ENTRY_SIZE,
    HEADER_SIZE,
    MAGIC,
    SECTION_BODY_INDEX,
    SECTION_CHEBYSHEV,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    VERSION,
    BodyEntry,
    FileHeader,
    SectionEntry,
    read_body_entry,
    read_header,
    read_section_dir,
    segment_byte_size,
    write_body_entry,
    write_header,
    write_section_dir,
)

NUM_SECTIONS = 5
J2000 = 2451545.0
CHUNK = 16 * 1024 * 1024


class LebFile:
    def __init__(self, path: str) -> None:
        self.path = Path(path)
        self.file = self.path.open("rb")
        self.mm = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)
        self.header = read_header(self.mm, 0)
        if self.header.magic != MAGIC or self.header.version != VERSION:
            raise ValueError(f"invalid LEB file: {self.path}")
        self.sections = {}
        for idx in range(self.header.section_count):
            sec = read_section_dir(self.mm, HEADER_SIZE + idx * SECTION_DIR_SIZE)
            self.sections[sec.section_id] = sec
        self.bodies = {}
        sec = self.sections[SECTION_BODY_INDEX]
        for idx in range(self.header.body_count):
            entry = read_body_entry(self.mm, sec.offset + idx * BODY_ENTRY_SIZE)
            self.bodies[entry.body_id] = entry

    def close(self) -> None:
        self.mm.close()
        self.file.close()


def copy_range(out, mm: mmap.mmap, offset: int, size: int) -> None:
    mm.seek(offset)
    remaining = size
    while remaining:
        n = min(CHUNK, remaining)
        out.write(mm.read(n))
        remaining -= n


main_path = sys.argv[1]
output_path = sys.argv[2]
partial_paths = sys.argv[3:]

main = LebFile(main_path)
partials = [LebFile(path) for path in partial_paths]
try:
    for partial in partials:
        if (
            abs(partial.header.jd_start - main.header.jd_start) > 0.5
            or abs(partial.header.jd_end - main.header.jd_end) > 0.5
        ):
            raise ValueError(
                f"JD range mismatch: {partial.path} has "
                f"{partial.header.jd_start:.1f}..{partial.header.jd_end:.1f}, "
                f"main has {main.header.jd_start:.1f}..{main.header.jd_end:.1f}"
            )

    source: dict[int, tuple[LebFile, BodyEntry]] = {
        body_id: (main, entry) for body_id, entry in main.bodies.items()
    }
    replaced = 0
    for partial in partials:
        for body_id, entry in partial.bodies.items():
            source[body_id] = (partial, entry)
            replaced += 1

    body_ids = sorted(source)
    body_index_size = len(body_ids) * BODY_ENTRY_SIZE
    cheb_size = sum(
        entry.segment_count * segment_byte_size(entry.degree, entry.components)
        for _, entry in source.values()
    )

    def section_size(section_id: int) -> int:
        sec = main.sections.get(section_id)
        return 0 if sec is None else sec.size

    nut_size = section_size(SECTION_NUTATION)
    dt_size = section_size(SECTION_DELTA_T)
    star_size = section_size(SECTION_STARS)

    body_index_offset = HEADER_SIZE + NUM_SECTIONS * SECTION_DIR_SIZE
    cheb_offset = body_index_offset + body_index_size
    nut_offset = cheb_offset + cheb_size
    dt_offset = nut_offset + nut_size
    star_offset = dt_offset + dt_size

    fd, tmp_name = tempfile.mkstemp(
        prefix=Path(output_path).name + ".",
        suffix=".tmp",
        dir=str(Path(output_path).parent),
    )
    os.close(fd)

    try:
        with open(tmp_name, "wb") as out:
            header_buf = bytearray(HEADER_SIZE + NUM_SECTIONS * SECTION_DIR_SIZE)
            header = FileHeader(
                magic=MAGIC,
                version=VERSION,
                section_count=NUM_SECTIONS,
                body_count=len(body_ids),
                jd_start=main.header.jd_start,
                jd_end=main.header.jd_end,
                generation_epoch=J2000 + (time.time() / 86400.0 - 10957.5),
                flags=main.header.flags,
            )
            write_header(header_buf, header)
            sections = [
                SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
                SectionEntry(SECTION_CHEBYSHEV, cheb_offset, cheb_size),
                SectionEntry(SECTION_NUTATION, nut_offset, nut_size),
                SectionEntry(SECTION_DELTA_T, dt_offset, dt_size),
                SectionEntry(SECTION_STARS, star_offset, star_size),
            ]
            for idx, sec in enumerate(sections):
                write_section_dir(header_buf, HEADER_SIZE + idx * SECTION_DIR_SIZE, sec)
            out.write(header_buf)

            body_index = bytearray(body_index_size)
            coeff_offset = cheb_offset
            coeff_sources: list[tuple[LebFile, int, int]] = []
            for idx, body_id in enumerate(body_ids):
                src, entry = source[body_id]
                size = entry.segment_count * segment_byte_size(
                    entry.degree, entry.components
                )
                new_entry = BodyEntry(
                    body_id=entry.body_id,
                    coord_type=entry.coord_type,
                    segment_count=entry.segment_count,
                    jd_start=entry.jd_start,
                    jd_end=entry.jd_end,
                    interval_days=entry.interval_days,
                    degree=entry.degree,
                    components=entry.components,
                    data_offset=coeff_offset,
                )
                write_body_entry(body_index, idx * BODY_ENTRY_SIZE, new_entry)
                coeff_sources.append((src, entry.data_offset, size))
                coeff_offset += size
            out.write(body_index)

            for src, offset, size in coeff_sources:
                copy_range(out, src.mm, offset, size)

            for section_id in (SECTION_NUTATION, SECTION_DELTA_T, SECTION_STARS):
                sec = main.sections.get(section_id)
                if sec and sec.size:
                    copy_range(out, main.mm, sec.offset, sec.size)

        os.replace(tmp_name, output_path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise

    print(
        f"updated {output_path}: bodies={len(body_ids)} "
        f"replaced_entries={replaced} partials={len(partials)}"
    )
finally:
    for partial in partials:
        partial.close()
    main.close()
PY
  statuses=("${PIPESTATUS[@]}")
  cmd_status="${statuses[0]}"
  tee_status="${statuses[1]:-0}"

  echo
  if ((cmd_status == 0)); then
    echo "DONE : replace selected LEB1 bodies in $output"
  else
    echo "FAILED: replace selected LEB1 bodies in $output"
    echo "EXIT CODE: $cmd_status"
    FAILED_COMMANDS+=("replace selected LEB1 bodies in $output")
    FAILED_STATUSES+=("$cmd_status")
  fi

  if ((tee_status != 0)); then
    echo "WARNING: failed to write command output to log file: $LOG_FILE"
  fi

  echo "TIME : $(date '+%Y-%m-%d %H:%M:%S')"
  return "$cmd_status"
}

merge_leb1_tier() {
  local tier="$1"
  local output
  local partials=()
  local group
  local args

  ((DO_MERGE)) || {
    echo "LEB1 merge skipped (--no-merge)."
    return 0
  }

  output="$(main_leb1_path "$tier")"
  if all_leb1_groups_selected; then
    while IFS= read -r group; do
      partials+=("$(partial_leb1_path "$tier" "$group")")
    done < <(selected_leb1_groups)
    backup_path_if_exists "$output"
    args=(scripts/generate_leb.py --tier "$tier" --merge "${partials[@]}" --output "$output")
    ((QUIET)) && args+=(--quiet)
    run "$PYTHON" "${args[@]}" || return 1
  else
    replace_leb1_from_partials "$tier" || return 1
  fi

  GENERATED_LEB1_MAIN+=("$output")
}

verify_leb1_tier() {
  local tier="$1"
  local group
  local target
  local args

  ((VERIFY)) || {
    echo "LEB1 verification skipped (--no-verify)."
    return 0
  }

  if ((DO_MERGE)); then
    target="$(main_leb1_path "$tier")"
    args=(scripts/generate_leb.py --tier "$tier" --verify-only --output "$target" --verify-samples "$VERIFY_SAMPLES")
    ((QUIET)) && args+=(--quiet)
    run "$PYTHON" "${args[@]}" || return 1
    return 0
  fi

  while IFS= read -r group; do
    target="$(partial_leb1_path "$tier" "$group")"
    args=(scripts/generate_leb.py --tier "$tier" --verify-only --output "$target" --verify-samples "$VERIFY_SAMPLES")
    ((QUIET)) && args+=(--quiet)
    run "$PYTHON" "${args[@]}" || return 1
  done < <(selected_leb1_groups)
}

generate_leb1_tier() {
  local tier="$1"

  echo
  echo "### LEB1: $tier"
  generate_leb1_groups_for_tier "$tier" || return 1
  merge_leb1_tier "$tier" || return 1
  verify_leb1_tier "$tier" || return 1
}

convert_leb2_tier() {
  local tier="$1"
  local group
  local source
  local output
  local args

  echo
  echo "### LEB2 conversion: $tier"

  mkdir -p data/leb2
  source="$(main_leb1_path "$tier")"
  if ! ((DRY_RUN)); then
    [[ -f "$source" ]] || die "cannot convert LEB2: missing LEB1 source: $source"

    # Freshness guard: LEB2 is converted from the merged main LEB1. If any
    # partial LEB1 group file for this tier is newer than the main (e.g. a
    # --no-merge run regenerated groups without updating the main), the LEB2
    # output would be built from stale data. Fail loudly rather than ship it.
    local newest_partial="" partial
    for partial in data/leb/ephemeris_"$tier"_*.leb; do
      [[ -e "$partial" ]] || continue
      if [[ -z "$newest_partial" || "$partial" -nt "$newest_partial" ]]; then
        newest_partial="$partial"
      fi
    done
    if [[ -n "$newest_partial" && "$newest_partial" -nt "$source" ]]; then
      die "cannot convert LEB2: main LEB1 '$source' is older than partial" \
        "group file '$newest_partial'; merge the LEB1 groups first (drop" \
        "--no-merge) so LEB2 is not built from a stale main."
    fi
  fi

  while IFS= read -r group; do
    output="$(leb2_path "$tier" "$group")"
    backup_path_if_exists "$output"
    # ``--tier`` is not merely output naming metadata.  The converter uses it
    # to authenticate the exact body inventory permitted for that tier (most
    # notably, the extended exotics set excludes the near-Earth asteroids).
    # Omitting it makes a valid extended file look like a malformed generic
    # exotics file, so always propagate the tier selected by this orchestrator.
    args=(scripts/generate_leb2.py convert "$source" -o "$output" --group "$group" --tier "$tier")
    ((QUIET)) && args+=(-q)
    run "$PYTHON" "${args[@]}" || return 1
    GENERATED_LEB2_FILES+=("$output")
  done < <(selected_leb2_groups)
}

verify_leb2_tier() {
  local tier="$1"
  local group
  local input
  local reference
  local args

  ((VERIFY)) || {
    echo "LEB2 verification skipped (--no-verify)."
    return 0
  }

  reference="$(main_leb1_path "$tier")"
  while IFS= read -r group; do
    input="$(leb2_path "$tier" "$group")"
    # Verification must apply the same authenticated group/tier contract as
    # conversion.  Supplying both values also prevents a correctly readable
    # file with the wrong companion inventory from passing unnoticed.
    args=(scripts/generate_leb2.py verify "$input" --reference "$reference" --samples "$LEB2_VERIFY_SAMPLES" --group "$group" --tier "$tier")
    ((QUIET)) && args+=(-q)
    run "$PYTHON" "${args[@]}" || return 1
  done < <(selected_leb2_groups)
}

convert_and_verify_leb2_tier() {
  local tier="$1"

  convert_leb2_tier "$tier" || return 1
  verify_leb2_tier "$tier" || return 1
}

install_file() {
  local src="$1"
  local dst="$2"
  local backup_dir="$INSTALL_DIR/backup-$RUN_ID"

  if ! ((DRY_RUN)); then
    [[ -f "$src" ]] || die "install source not found: $src"
  fi
  mkdir -p "$INSTALL_DIR"
  if [[ -e "$dst" ]]; then
    mkdir -p "$backup_dir"
    if ((DRY_RUN)); then
      echo "DRY RUN: cp -p $dst $backup_dir/$(basename "$dst")"
    else
      cp -p "$dst" "$backup_dir/$(basename "$dst")"
      echo "Install backup: $dst -> $backup_dir/$(basename "$dst")"
    fi
  fi

  if ((DRY_RUN)); then
    echo "DRY RUN: cp -p $src $dst"
  else
    cp -p "$src" "$dst"
    echo "Installed: $src -> $dst"
  fi
}

install_outputs_for_tier() {
  local tier="$1"
  local group
  local src

  ((INSTALL)) || return 0

  echo
  echo "### Install: $tier -> $INSTALL_DIR"

  if ((DO_LEB1)) && ((DO_MERGE)); then
    src="$(main_leb1_path "$tier")"
    install_file "$src" "$INSTALL_DIR/$(basename "$src")" || return 1
  fi

  if ((DO_LEB2)); then
    while IFS= read -r group; do
      src="$(leb2_path "$tier" "$group")"
      install_file "$src" "$INSTALL_DIR/$(basename "$src")" || return 1
    done < <(selected_leb2_groups)
  fi
}

print_configuration() {
  local leb1_groups
  local leb2_groups

  leb1_groups="$(selected_leb1_groups | paste -sd, -)"
  leb2_groups="$(selected_leb2_groups | paste -sd, -)"

  echo "Repo: $SCRIPT_DIR"
  echo "Python: $PYTHON"
  echo "Tier: $TIER"
  echo "LEB1: $([[ "$DO_LEB1" == 1 ]] && echo yes || echo no)"
  echo "LEB1 groups: $leb1_groups"
  echo "Merge: $([[ "$DO_MERGE" == 1 ]] && echo yes || echo no)"
  echo "LEB2: $([[ "$DO_LEB2" == 1 ]] && echo yes || echo no)"
  echo "LEB2 groups: $leb2_groups"
  echo "Verify: $([[ "$VERIFY" == 1 ]] && echo yes || echo no)"
  echo "Install: $([[ "$INSTALL" == 1 ]] && echo "$INSTALL_DIR" || echo no)"
  echo "Dry run: $([[ "$DRY_RUN" == 1 ]] && echo yes || echo no)"
  echo "Log: $LOG_FILE"
}

main() {
  local tier

  parse_args "$@"
  mkdir -p "$LOG_DIR"
  trap cleanup EXIT

  while ! mkdir "$LOCK_DIR" 2>/dev/null; do
    echo "WARNING: lock exists: $LOCK_DIR"
    echo "Waiting ${LOCK_WAIT_SECONDS}s before retrying."
    sleep "$LOCK_WAIT_SECONDS"
  done
  # Lock is now held by this process; the EXIT trap may remove it.
  LOCK_ACQUIRED=1

  select_python
  print_configuration

  backup_problem_tno_spks | tee -a "$LOG_FILE"

  while IFS= read -r tier; do
    if ((DO_LEB1)); then
      generate_leb1_tier "$tier"
      if ((${#FAILED_COMMANDS[@]})); then
        break
      fi
    fi

    if ((DO_LEB2)); then
      convert_and_verify_leb2_tier "$tier"
      if ((${#FAILED_COMMANDS[@]})); then
        break
      fi
    fi

    install_outputs_for_tier "$tier"
    if ((${#FAILED_COMMANDS[@]})); then
      break
    fi
  done < <(selected_tiers)

  echo
  if ((${#FAILED_COMMANDS[@]} == 0)); then
    echo "SUCCESS: LEB regeneration completed."
    echo "Log: $LOG_FILE"
    exit 0
  fi

  echo "FAILED:"
  for i in "${!FAILED_COMMANDS[@]}"; do
    echo "  - ${FAILED_COMMANDS[$i]} (exit ${FAILED_STATUSES[$i]})"
  done
  echo "Log: $LOG_FILE"
  exit 1
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
