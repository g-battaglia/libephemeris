#!/usr/bin/env bash
#
# Regenerate, verify, convert, and optionally install LibEphemeris LEB data.
#
# Compatibility:
#   - Bash 3.2+ (including the version shipped with macOS)
#   - No associative arrays, namerefs, or mapfile
#
# Error handling note:
#   `errexit` (`set -e`) is intentionally NOT enabled. Every external command
#   is executed through `run_logged`, which records failures and lets the main
#   workflow print one clear summary before exiting.

set -u -o pipefail

# ==============================================================================
# Paths and immutable project metadata
# ==============================================================================

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

readonly LOG_DIR="logs"
readonly RUN_ID="$(date '+%Y%m%d-%H%M%S')"
readonly LOG_FILE="$LOG_DIR/leb-regeneration-$RUN_ID.log"
readonly LOCK_DIR=".leb-regeneration.lock"
readonly LOCK_WAIT_SECONDS="${LOCK_WAIT_SECONDS:-60}"

# Supported values. Keep LEB2 groups in the same order as
# `libephemeris.leb_groups.LEB2_GROUPS`.
TIERS=(base medium extended)
LEB1_GROUPS=(planets asteroids exotics analytical)
LEB2_GROUPS=(core asteroids exotics apogee)

# These seven TNO kernels must be refreshed before exotics regeneration unless
# the operator explicitly opts into reusing the SPK cache.
TNO_SPK_IDS=(
  136199  # Eris
  90377   # Sedna
  136108  # Haumea
  136472  # Makemake
  28978   # Ixion
  90482   # Orcus
  50000   # Quaoar
)
readonly TIERS LEB1_GROUPS LEB2_GROUPS TNO_SPK_IDS

# ==============================================================================
# Terminal formatting and UI helpers
# Colors are emitted only on an interactive terminal and are disabled when the
# NO_COLOR convention (https://no-color.org) is set, so redirected output and
# log files stay clean. Bash 3.2 compatible: plain ANSI escapes, no tput.
# ==============================================================================

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  C_RESET=$'\033[0m'
  C_BOLD=$'\033[1m'
  C_DIM=$'\033[2m'
  C_CYAN=$'\033[36m'
  C_GREEN=$'\033[32m'
  C_YELLOW=$'\033[33m'
  C_RED=$'\033[31m'
else
  C_RESET=""
  C_BOLD=""
  C_DIM=""
  C_CYAN=""
  C_GREEN=""
  C_YELLOW=""
  C_RED=""
fi
readonly C_RESET C_BOLD C_DIM C_CYAN C_GREEN C_YELLOW C_RED

RULE="────────────────────────────────────────────────────────────────────────"

# Dim horizontal divider.
print_rule() {
  printf '%s%s%s\n' "$C_DIM" "$RULE" "$C_RESET"
}

# Major stage header: blank line, bold title with a marker, an optional dim
# subtitle, then a divider. Use for top-level workflow stages.
print_header() {
  local title="$1"
  local subtitle="${2:-}"

  echo
  printf '%s● %s%s\n' "$C_BOLD$C_CYAN" "$title" "$C_RESET"
  if [[ -n "$subtitle" ]]; then
    printf '%s  %s%s\n' "$C_DIM" "$subtitle" "$C_RESET"
  fi
  print_rule
}

# Lighter phase header for a group of commands inside a major stage (no divider).
print_subheader() {
  local title="$1"

  echo
  printf '%s▸ %s%s\n' "$C_BOLD$C_CYAN" "$title" "$C_RESET"
}

# Dim informational note: backups, cache moves, skips.
print_note() {
  printf '%s· %s%s\n' "$C_DIM" "$*" "$C_RESET"
}

# Dim, tagged line for an action that --dry-run would have taken.
print_dry_run() {
  printf '%s· dry-run  %s%s\n' "$C_DIM" "$*" "$C_RESET"
}

# ==============================================================================
# Mutable configuration (populated by parse_arguments)
# ==============================================================================

TIER_SELECTION="all"
LEB1_GROUP_SELECTION="all"
LEB2_GROUP_SELECTION=""

ENABLE_LEB1=1
ENABLE_LEB2=1
ENABLE_MERGE=1
ENABLE_VERIFICATION=1
ENABLE_INSTALL=0
KEEP_SPK_CACHE=0
DRY_RUN=0
QUIET=0
DOCTOR_ONLY=0
SKIP_PREFLIGHT=0

LEB1_VERIFY_SAMPLES=500
LEB2_VERIFY_SAMPLES=200
INSTALL_DIR="${LIBEPHEMERIS_DATA_DIR:-$HOME/.libephemeris}/leb"

# Set after argument parsing. Resolving them once avoids repeatedly spawning
# subshells and guarantees that every workflow phase uses the same selection.
SELECTED_TIERS=()
SELECTED_LEB1_GROUPS=()
SELECTED_LEB2_GROUPS=()

# Runtime state.
PYTHON=""
LOCK_ACQUIRED=0
FAILED_COMMANDS=()
FAILED_STATUSES=()

# ==============================================================================
# User interface and argument parsing
# ==============================================================================

print_usage() {
  cat <<'EOF'
Usage:
  ./regenerate-leb.sh [base|medium|extended|all] [options]
  ./regenerate-leb.sh --tier base|medium|extended|all [options]

Regenerate LibEphemeris LEB1 and LEB2 files. This is the single safe entrypoint
for complete tier regeneration and partial group updates such as “medium
exotics only”.

Default:
  ./regenerate-leb.sh all

Default workflow for each selected tier:
  1. Refresh cached SPKs for the seven refitted TNOs when exotics are selected.
  2. Generate the selected LEB1 partial group files.
  3. Rebuild the complete LEB1 file or replace only the selected groups.
  4. Verify LEB1 unless --no-verify is used.
  5. Convert the affected LEB2 groups unless --leb1-only is used.
  6. Verify LEB2 unless --no-verify is used.

Tier selection:
  --tier TIER             base, medium, extended, or all. Default: all.
  base|medium|extended|all may also be passed positionally.

LEB1 group selection:
  --group GROUP           planets, asteroids, exotics, analytical, or all.
  --groups LIST           Comma-separated groups, e.g. planets,exotics.
  --no-merge              Generate partial LEB1 files without updating the main
                          data/leb/ephemeris_<tier>.leb file.

LEB2 selection:
  --leb1-only             Generate, merge, and verify LEB1 only.
  --leb2-only             Convert and verify LEB2 from existing LEB1 files only.
  --leb2-groups LIST      Comma-separated groups: core, asteroids, exotics,
                          apogee, or all. When omitted, groups are inferred from
                          the selected LEB1 groups.

Verification:
  --no-verify             Skip LEB1 and LEB2 verification.
  --verify-samples N      LEB1 samples per body. Default: 500.
  --leb2-verify-samples N LEB2 samples per body. Default: 200.

SPK cache:
  --keep-spk-cache        Reuse cached SPKs during exotics generation. Use only
                          after auditing ~/.libephemeris/spk.

Install:
  --install               Install generated main LEB1 and affected LEB2 files.
  --install-dir DIR       Override the installation directory.

Preflight:
  --doctor                Run only dependency/data checks and exit. No lock and
                          no generation. Exit 0 when ready, 2 when items are
                          missing.
  --skip-deps-check       Skip the normal preflight. Use only when the execution
                          environment is already known to be complete.

Other:
  --dry-run               Print planned actions without generating or copying.
  -q, --quiet             Pass quiet mode to generators where supported.
  -h, --help              Show this help.

Examples:
  ./regenerate-leb.sh medium --group exotics --install
  ./regenerate-leb.sh base --groups planets,asteroids
  ./regenerate-leb.sh all --group all --no-verify
  ./regenerate-leb.sh medium --leb2-only --leb2-groups exotics --install
  ./regenerate-leb.sh extended --group exotics --leb1-only
  ./regenerate-leb.sh --doctor
  ./regenerate-leb.sh base --doctor

Notes:
  - Extended generation can take hours.
  - Extended exotics require rebound, assist, and the ASSIST data files.
  - LEB1 “planets” and LEB2 “core” are not equivalent groups. LEB2 group
    inference is deliberately conservative.
EOF
}

fatal() {
  printf '%sERROR: %s%s\n' "$C_BOLD$C_RED" "$*" "$C_RESET" >&2
  exit 2
}

require_option_value() {
  local option_name="$1"
  local remaining_argument_count="$2"

  ((remaining_argument_count >= 2)) || fatal "$option_name requires a value"
}

normalize_list() {
  local raw_value="$1"

  raw_value="${raw_value// /}"
  raw_value="${raw_value//;/,}"
  printf '%s\n' "$raw_value"
}

parse_arguments() {
  while (($#)); do
    case "$1" in
      --tier)
        require_option_value "$1" "$#"
        TIER_SELECTION="$2"
        shift 2
        ;;
      --tier=*)
        TIER_SELECTION="${1#--tier=}"
        shift
        ;;

      --group|--groups)
        require_option_value "$1" "$#"
        LEB1_GROUP_SELECTION="$2"
        shift 2
        ;;
      --group=*|--groups=*)
        LEB1_GROUP_SELECTION="${1#*=}"
        shift
        ;;

      --leb2-groups)
        require_option_value "$1" "$#"
        LEB2_GROUP_SELECTION="$2"
        shift 2
        ;;
      --leb2-groups=*)
        LEB2_GROUP_SELECTION="${1#--leb2-groups=}"
        shift
        ;;

      --leb1-only)
        ENABLE_LEB1=1
        ENABLE_LEB2=0
        shift
        ;;
      --leb2-only)
        ENABLE_LEB1=0
        ENABLE_LEB2=1
        shift
        ;;
      --no-merge)
        ENABLE_MERGE=0
        shift
        ;;
      --no-verify)
        ENABLE_VERIFICATION=0
        shift
        ;;

      --verify-samples)
        require_option_value "$1" "$#"
        LEB1_VERIFY_SAMPLES="$2"
        shift 2
        ;;
      --verify-samples=*)
        LEB1_VERIFY_SAMPLES="${1#--verify-samples=}"
        shift
        ;;
      --leb2-verify-samples)
        require_option_value "$1" "$#"
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
        ENABLE_INSTALL=1
        shift
        ;;
      --install-dir)
        require_option_value "$1" "$#"
        INSTALL_DIR="$2"
        shift 2
        ;;
      --install-dir=*)
        INSTALL_DIR="${1#--install-dir=}"
        shift
        ;;

      --doctor)
        DOCTOR_ONLY=1
        shift
        ;;
      --skip-deps-check)
        SKIP_PREFLIGHT=1
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
        print_usage
        exit 0
        ;;

      base|medium|extended|all)
        TIER_SELECTION="$1"
        shift
        ;;
      *)
        fatal "unknown argument: $1"
        ;;
    esac
  done
}

# ==============================================================================
# Configuration resolution
# ==============================================================================

append_unique_leb1_group() {
  local candidate="$1"
  local existing

  if ((${#SELECTED_LEB1_GROUPS[@]} > 0)); then
    for existing in "${SELECTED_LEB1_GROUPS[@]}"; do
      [[ "$existing" == "$candidate" ]] && return 0
    done
  fi
  SELECTED_LEB1_GROUPS+=("$candidate")
}

append_unique_leb2_group() {
  local candidate="$1"
  local existing

  if ((${#SELECTED_LEB2_GROUPS[@]} > 0)); then
    for existing in "${SELECTED_LEB2_GROUPS[@]}"; do
      [[ "$existing" == "$candidate" ]] && return 0
    done
  fi
  SELECTED_LEB2_GROUPS+=("$candidate")
}

resolve_tier_selection() {
  case "$TIER_SELECTION" in
    all)
      SELECTED_TIERS=("${TIERS[@]}")
      ;;
    base|medium|extended)
      SELECTED_TIERS=("$TIER_SELECTION")
      ;;
    *)
      fatal "invalid tier: $TIER_SELECTION"
      ;;
  esac
}

resolve_leb1_group_selection() {
  local normalized_selection
  local requested_group
  local candidate
  local requested_groups=()
  local candidates=()

  normalized_selection="$(normalize_list "$LEB1_GROUP_SELECTION")"

  if [[ -z "$normalized_selection" || "$normalized_selection" == "all" ]]; then
    SELECTED_LEB1_GROUPS=("${LEB1_GROUPS[@]}")
    return 0
  fi

  IFS=',' read -r -a requested_groups <<< "$normalized_selection"
  for requested_group in "${requested_groups[@]}"; do
    case "$requested_group" in
      planets|asteroids|exotics|analytical)
        candidates=("$requested_group")
        ;;
      all)
        candidates=("${LEB1_GROUPS[@]}")
        ;;
      *)
        fatal "invalid LEB1 group: $requested_group"
        ;;
    esac

    for candidate in "${candidates[@]}"; do
      append_unique_leb1_group "$candidate"
    done
  done
}

infer_leb2_groups_from_leb1() {
  local leb1_group

  for leb1_group in "${SELECTED_LEB1_GROUPS[@]}"; do
    case "$leb1_group" in
      planets)
        append_unique_leb2_group core
        ;;
      asteroids)
        append_unique_leb2_group asteroids
        ;;
      exotics)
        append_unique_leb2_group exotics
        ;;
      analytical)
        # Public analytical bodies are split across LEB2 core and apogee.
        # The historical “uranians” output was intentionally retired.
        append_unique_leb2_group core
        append_unique_leb2_group apogee
        ;;
    esac
  done
}

resolve_leb2_group_selection() {
  local normalized_selection
  local requested_group
  local candidate
  local requested_groups=()
  local candidates=()

  normalized_selection="$(normalize_list "$LEB2_GROUP_SELECTION")"

  if [[ -z "$normalized_selection" ]]; then
    infer_leb2_groups_from_leb1
    return 0
  fi

  if [[ "$normalized_selection" == "all" ]]; then
    SELECTED_LEB2_GROUPS=("${LEB2_GROUPS[@]}")
    return 0
  fi

  IFS=',' read -r -a requested_groups <<< "$normalized_selection"
  for requested_group in "${requested_groups[@]}"; do
    case "$requested_group" in
      core|asteroids|exotics|apogee)
        candidates=("$requested_group")
        ;;
      all)
        candidates=("${LEB2_GROUPS[@]}")
        ;;
      *)
        fatal "invalid LEB2 group: $requested_group"
        ;;
    esac

    for candidate in "${candidates[@]}"; do
      append_unique_leb2_group "$candidate"
    done
  done
}

resolve_configuration() {
  [[ "$LEB1_VERIFY_SAMPLES" =~ ^[0-9]+$ ]] \
    || fatal "--verify-samples must be an integer"
  [[ "$LEB2_VERIFY_SAMPLES" =~ ^[0-9]+$ ]] \
    || fatal "--leb2-verify-samples must be an integer"

  resolve_tier_selection
  resolve_leb1_group_selection
  resolve_leb2_group_selection
}

array_contains() {
  local needle="$1"
  shift

  local item
  for item in "$@"; do
    [[ "$item" == "$needle" ]] && return 0
  done
  return 1
}

all_leb1_groups_selected() {
  [[ "${#SELECTED_LEB1_GROUPS[@]}" -eq "${#LEB1_GROUPS[@]}" ]]
}

exotics_selected() {
  array_contains exotics "${SELECTED_LEB1_GROUPS[@]}"
}

nbody_required() {
  exotics_selected || return 1
  array_contains extended "${SELECTED_TIERS[@]}"
}

spk_network_needed() {
  array_contains asteroids "${SELECTED_LEB1_GROUPS[@]}" \
    || array_contains exotics "${SELECTED_LEB1_GROUPS[@]}"
}

join_with_commas() {
  local IFS=','
  printf '%s' "$*"
}

yes_or_no() {
  local enabled="$1"
  ((enabled)) && printf 'yes' || printf 'no'
}

# ==============================================================================
# Process execution, logging, and failure tracking
# ==============================================================================

format_command() {
  local formatted=""

  printf -v formatted '%q ' "$@"
  printf '%s' "${formatted% }"
}

record_failure() {
  local description="$1"
  local status="$2"

  FAILED_COMMANDS+=("$description")
  FAILED_STATUSES+=("$status")
}

has_failures() {
  ((${#FAILED_COMMANDS[@]} > 0))
}

print_action_header() {
  local description="$1"

  echo
  printf '%s›%s %s%s%s\n' "$C_DIM" "$C_RESET" "$C_BOLD" "$description" "$C_RESET"
  if ((DRY_RUN)); then
    printf '%s  dry-run — no action taken%s\n' "$C_DIM$C_CYAN" "$C_RESET"
  else
    printf '%s  start %s%s\n' "$C_DIM" "$(date '+%H:%M:%S')" "$C_RESET"
  fi
}

# Execute a command, mirror stdout/stderr to the terminal and log file, preserve
# the command's real exit status through `tee`, and add failures to the final
# summary. The optional explicit description is useful for heredoc commands.
run_logged_as() {
  local description="$1"
  shift

  local pipeline_statuses=()
  local command_status
  local tee_status

  print_action_header "$description"

  if ((DRY_RUN)); then
    return 0
  fi

  "$@" 2>&1 | tee -a "$LOG_FILE"
  pipeline_statuses=("${PIPESTATUS[@]}")
  command_status="${pipeline_statuses[0]}"
  tee_status="${pipeline_statuses[1]:-0}"

  echo
  if ((command_status == 0)); then
    printf '%s✓ done%s   %s%s%s\n' \
      "$C_GREEN$C_BOLD" "$C_RESET" "$C_DIM" "$(date '+%H:%M:%S')" "$C_RESET"
  else
    printf '%s✗ failed%s (exit %d)   %s%s%s\n' \
      "$C_RED$C_BOLD" "$C_RESET" "$command_status" "$C_DIM" "$(date '+%H:%M:%S')" "$C_RESET"
    record_failure "$description" "$command_status"
  fi

  if ((tee_status != 0)); then
    printf '%s⚠ could not append output to %s%s\n' "$C_YELLOW" "$LOG_FILE" "$C_RESET"
  fi

  return "$command_status"
}

run_logged() {
  local description

  description="$(format_command "$@")"
  run_logged_as "$description" "$@"
}

select_python_interpreter() {
  if [[ -x ".venv/bin/python" ]]; then
    PYTHON=".venv/bin/python"
  elif command -v python >/dev/null 2>&1; then
    PYTHON="python"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON="python3"
  else
    fatal "no Python interpreter found"
  fi

  # Ensure local imports resolve to this repository before any installed copy.
  export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
}

# ==============================================================================
# Lock lifecycle
# ==============================================================================

cleanup_lock() {
  # Never remove another process's lock. The flag is set only after our mkdir
  # succeeds, so an interrupted waiter cannot delete the active holder's lock.
  if ((LOCK_ACQUIRED)); then
    rm -rf "$LOCK_DIR"
  fi
}

acquire_lock() {
  while ! mkdir "$LOCK_DIR" 2>/dev/null; do
    printf '%s⚠ lock exists: %s — waiting %ss before retrying.%s\n' \
      "$C_YELLOW" "$LOCK_DIR" "$LOCK_WAIT_SECONDS" "$C_RESET"
    sleep "$LOCK_WAIT_SECONDS"
  done

  LOCK_ACQUIRED=1
}

# ==============================================================================
# File names, backups, and cache management
# ==============================================================================

leb1_main_path() {
  local tier="$1"
  printf 'data/leb/ephemeris_%s.leb\n' "$tier"
}

leb1_partial_path() {
  local tier="$1"
  local group="$2"
  printf 'data/leb/ephemeris_%s_%s.leb\n' "$tier" "$group"
}

leb2_output_path() {
  local tier="$1"
  local group="$2"
  printf 'data/leb2/%s_%s.leb2\n' "$tier" "$group"
}

ensure_directory() {
  local directory="$1"

  if ((DRY_RUN)); then
    print_dry_run "mkdir -p $(printf '%q' "$directory")"
  else
    mkdir -p "$directory"
  fi
}

backup_file_if_present() {
  local source_path="$1"
  local backup_directory="$LOG_DIR/leb-backup-$RUN_ID"
  local backup_path="$backup_directory/$(basename "$source_path")"

  [[ -e "$source_path" ]] || return 0

  if ((DRY_RUN)); then
    print_dry_run "cp -p $(printf '%q' "$source_path") $(printf '%q' "$backup_path")"
    return 0
  fi

  mkdir -p "$backup_directory"
  cp -p "$source_path" "$backup_path"
  print_note "backed up $source_path → $backup_path"
}

list_spk_cache_directories() {
  local data_directory="${LIBEPHEMERIS_DATA_DIR:-$HOME/.libephemeris}"
  local directories=("$data_directory/spk")

  if [[ -n "${LIBEPHEMERIS_SPK_DIR:-}" ]]; then
    directories+=("$LIBEPHEMERIS_SPK_DIR")
  fi

  printf '%s\n' "${directories[@]}" | awk '!seen[$0]++'
}

refresh_refitted_tno_spk_cache() {
  local backup_directory="$LOG_DIR/spk-backup-$RUN_ID"
  local moved_file_count=0
  local spk_directory
  local object_id
  local spk_path

  if ((KEEP_SPK_CACHE)); then
    print_note "SPK cache refresh skipped (--keep-spk-cache)"
    return 0
  fi

  if ! ((ENABLE_LEB1)) || ! exotics_selected; then
    print_note "SPK cache refresh skipped (exotics LEB1 not selected)"
    return 0
  fi

  # nullglob prevents an unmatched pattern from being treated as a literal path.
  shopt -s nullglob
  while IFS= read -r spk_directory; do
    [[ -d "$spk_directory" ]] || continue

    for object_id in "${TNO_SPK_IDS[@]}"; do
      for spk_path in "$spk_directory"/"${object_id}"_*.bsp; do
        [[ -e "$spk_path" ]] || continue

        if ((DRY_RUN)); then
          print_dry_run "mv $(printf '%q' "$spk_path") $(printf '%q' "$backup_directory/$(basename "$spk_path")")"
        else
          mkdir -p "$backup_directory"
          mv "$spk_path" "$backup_directory/$(basename "$spk_path")"
        fi
        moved_file_count=$((moved_file_count + 1))
      done
    done
  done < <(list_spk_cache_directories)
  shopt -u nullglob

  if ((moved_file_count == 0)); then
    print_note "no cached SPKs found for the seven refitted TNOs"
  else
    print_note "moved $moved_file_count SPK file(s) to $backup_directory"
    print_note "generation will download fresh JPL Horizons kernels as needed"
  fi
}

# ==============================================================================
# LEB1 generation, merging, and verification
# ==============================================================================

generate_leb1_partials_for_tier() {
  local tier="$1"
  local group
  local output_path
  local command_args=()

  print_subheader "LEB1 partial generation — $tier"
  ensure_directory data/leb

  for group in "${SELECTED_LEB1_GROUPS[@]}"; do
    output_path="$(leb1_partial_path "$tier" "$group")"
    command_args=(
      scripts/generate_leb.py
      --tier "$tier"
      --group "$group"
      --output "$output_path"
    )
    ((QUIET)) && command_args+=(--quiet)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
  done
}

# Replace bodies from selected partial files inside an existing main LEB1.
# The embedded Python implementation writes a complete temporary file and uses
# os.replace(), so an interrupted update never leaves a half-written main file.
replace_leb1_groups_from_partials() {
  local tier="$1"
  local main_path
  local group
  local partial_path
  local partial_paths=()
  local description

  main_path="$(leb1_main_path "$tier")"

  if ! ((DRY_RUN)); then
    [[ -f "$main_path" ]] \
      || fatal "cannot update selected groups: missing main LEB1: $main_path"
  fi

  for group in "${SELECTED_LEB1_GROUPS[@]}"; do
    partial_path="$(leb1_partial_path "$tier" "$group")"
    partial_paths+=("$partial_path")

    if ! ((DRY_RUN)); then
      [[ -f "$partial_path" ]] || fatal "missing generated partial: $partial_path"
    fi
  done

  backup_file_if_present "$main_path"
  description="replace selected LEB1 bodies in $main_path"

  run_logged_as "$description" \
    "$PYTHON" - "$main_path" "$main_path" "${partial_paths[@]}" <<'PY'
from __future__ import annotations

"""Atomically replace selected body entries in a main LEB1 file.

The shell orchestrator passes the current main file as both input and output,
followed by one or more freshly generated partial LEB1 files. Metadata sections
(nutation, Delta-T, and stars) are preserved from the main file, while body
entries found in partials replace entries with the same body ID.
"""

import mmap
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import BinaryIO, Iterable

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

SECTION_COUNT = 5
J2000_JULIAN_DAY = 2451545.0
COPY_CHUNK_SIZE = 16 * 1024 * 1024
PRESERVED_SECTION_IDS = (SECTION_NUTATION, SECTION_DELTA_T, SECTION_STARS)


class LebFileReader:
    """Memory-mapped, read-only view of the LEB structures needed here."""

    def __init__(self, path: str) -> None:
        self.path = Path(path)
        self._file = self.path.open("rb")
        self.memory_map = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

        self.header = read_header(self.memory_map, 0)
        if self.header.magic != MAGIC or self.header.version != VERSION:
            self.close()
            raise ValueError(f"invalid LEB file: {self.path}")

        self.sections: dict[int, SectionEntry] = {}
        for index in range(self.header.section_count):
            directory_offset = HEADER_SIZE + index * SECTION_DIR_SIZE
            section = read_section_dir(self.memory_map, directory_offset)
            self.sections[section.section_id] = section

        body_index_section = self.sections[SECTION_BODY_INDEX]
        self.bodies: dict[int, BodyEntry] = {}
        for index in range(self.header.body_count):
            entry_offset = body_index_section.offset + index * BODY_ENTRY_SIZE
            body_entry = read_body_entry(self.memory_map, entry_offset)
            self.bodies[body_entry.body_id] = body_entry

    def close(self) -> None:
        self.memory_map.close()
        self._file.close()

    def __enter__(self) -> "LebFileReader":
        return self

    def __exit__(self, _exc_type, _exc, _traceback) -> None:
        self.close()


def copy_mapped_range(
    destination: BinaryIO,
    source_map: mmap.mmap,
    offset: int,
    size: int,
) -> None:
    """Copy a large mmap range without loading the whole section into RAM."""

    source_map.seek(offset)
    remaining_bytes = size

    while remaining_bytes:
        chunk_size = min(COPY_CHUNK_SIZE, remaining_bytes)
        destination.write(source_map.read(chunk_size))
        remaining_bytes -= chunk_size


def body_data_size(entry: BodyEntry) -> int:
    """Return the serialized Chebyshev coefficient size for one body."""

    bytes_per_segment = segment_byte_size(entry.degree, entry.components)
    return entry.segment_count * bytes_per_segment


def validate_partial_ranges(
    main_file: LebFileReader,
    partial_files: Iterable[LebFileReader],
) -> None:
    """Reject partial files generated for a materially different JD range."""

    for partial_file in partial_files:
        start_mismatch = abs(partial_file.header.jd_start - main_file.header.jd_start)
        end_mismatch = abs(partial_file.header.jd_end - main_file.header.jd_end)

        if start_mismatch > 0.5 or end_mismatch > 0.5:
            raise ValueError(
                f"JD range mismatch: {partial_file.path} has "
                f"{partial_file.header.jd_start:.1f}.."
                f"{partial_file.header.jd_end:.1f}, main has "
                f"{main_file.header.jd_start:.1f}.."
                f"{main_file.header.jd_end:.1f}"
            )


def section_size(main_file: LebFileReader, section_id: int) -> int:
    section = main_file.sections.get(section_id)
    return 0 if section is None else section.size


def build_body_source_map(
    main_file: LebFileReader,
    partial_files: Iterable[LebFileReader],
) -> tuple[dict[int, tuple[LebFileReader, BodyEntry]], int]:
    """Overlay partial body entries on top of the main file's body entries."""

    body_sources = {
        body_id: (main_file, entry)
        for body_id, entry in main_file.bodies.items()
    }

    replaced_entry_count = 0
    for partial_file in partial_files:
        for body_id, entry in partial_file.bodies.items():
            body_sources[body_id] = (partial_file, entry)
            replaced_entry_count += 1

    return body_sources, replaced_entry_count


def write_rebuilt_leb(
    main_file: LebFileReader,
    body_sources: dict[int, tuple[LebFileReader, BodyEntry]],
    output_path: Path,
) -> None:
    """Serialize the rebuilt file to a sibling temp file, then replace output."""

    ordered_body_ids = sorted(body_sources)
    body_index_size = len(ordered_body_ids) * BODY_ENTRY_SIZE
    chebyshev_size = sum(
        body_data_size(entry) for _, entry in body_sources.values()
    )

    nutation_size = section_size(main_file, SECTION_NUTATION)
    delta_t_size = section_size(main_file, SECTION_DELTA_T)
    stars_size = section_size(main_file, SECTION_STARS)

    body_index_offset = HEADER_SIZE + SECTION_COUNT * SECTION_DIR_SIZE
    chebyshev_offset = body_index_offset + body_index_size
    nutation_offset = chebyshev_offset + chebyshev_size
    delta_t_offset = nutation_offset + nutation_size
    stars_offset = delta_t_offset + delta_t_size

    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=output_path.name + ".",
        suffix=".tmp",
        dir=str(output_path.parent),
    )
    os.close(file_descriptor)
    temporary_path = Path(temporary_name)

    try:
        with temporary_path.open("wb") as destination:
            header_and_directory = bytearray(
                HEADER_SIZE + SECTION_COUNT * SECTION_DIR_SIZE
            )

            generation_epoch = J2000_JULIAN_DAY + (
                time.time() / 86400.0 - 10957.5
            )
            output_header = FileHeader(
                magic=MAGIC,
                version=VERSION,
                section_count=SECTION_COUNT,
                body_count=len(ordered_body_ids),
                jd_start=main_file.header.jd_start,
                jd_end=main_file.header.jd_end,
                generation_epoch=generation_epoch,
                flags=main_file.header.flags,
            )
            write_header(header_and_directory, output_header)

            output_sections = (
                SectionEntry(
                    SECTION_BODY_INDEX,
                    body_index_offset,
                    body_index_size,
                ),
                SectionEntry(
                    SECTION_CHEBYSHEV,
                    chebyshev_offset,
                    chebyshev_size,
                ),
                SectionEntry(SECTION_NUTATION, nutation_offset, nutation_size),
                SectionEntry(SECTION_DELTA_T, delta_t_offset, delta_t_size),
                SectionEntry(SECTION_STARS, stars_offset, stars_size),
            )
            for index, section in enumerate(output_sections):
                directory_offset = HEADER_SIZE + index * SECTION_DIR_SIZE
                write_section_dir(
                    header_and_directory,
                    directory_offset,
                    section,
                )

            destination.write(header_and_directory)

            # Rebuild the body index with new contiguous coefficient offsets.
            body_index_buffer = bytearray(body_index_size)
            next_coefficient_offset = chebyshev_offset
            coefficient_sources: list[tuple[LebFileReader, int, int]] = []

            for index, body_id in enumerate(ordered_body_ids):
                source_file, source_entry = body_sources[body_id]
                coefficient_size = body_data_size(source_entry)

                output_entry = BodyEntry(
                    body_id=source_entry.body_id,
                    coord_type=source_entry.coord_type,
                    segment_count=source_entry.segment_count,
                    jd_start=source_entry.jd_start,
                    jd_end=source_entry.jd_end,
                    interval_days=source_entry.interval_days,
                    degree=source_entry.degree,
                    components=source_entry.components,
                    data_offset=next_coefficient_offset,
                )
                write_body_entry(
                    body_index_buffer,
                    index * BODY_ENTRY_SIZE,
                    output_entry,
                )

                coefficient_sources.append(
                    (
                        source_file,
                        source_entry.data_offset,
                        coefficient_size,
                    )
                )
                next_coefficient_offset += coefficient_size

            destination.write(body_index_buffer)

            # Copy each body's coefficients from either the old main file or a
            # selected partial, according to the overlay built above.
            for source_file, offset, size in coefficient_sources:
                copy_mapped_range(
                    destination,
                    source_file.memory_map,
                    offset,
                    size,
                )

            # Non-body metadata always comes from the existing main file.
            for section_id in PRESERVED_SECTION_IDS:
                section = main_file.sections.get(section_id)
                if section is not None and section.size:
                    copy_mapped_range(
                        destination,
                        main_file.memory_map,
                        section.offset,
                        section.size,
                    )

        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def replace_selected_bodies(
    main_path: str,
    output_path: str,
    partial_paths: list[str],
) -> None:
    main_file = LebFileReader(main_path)
    partial_files = [LebFileReader(path) for path in partial_paths]

    try:
        validate_partial_ranges(main_file, partial_files)
        body_sources, replaced_entry_count = build_body_source_map(
            main_file,
            partial_files,
        )
        write_rebuilt_leb(main_file, body_sources, Path(output_path))

        print(
            f"updated {output_path}: bodies={len(body_sources)} "
            f"replaced_entries={replaced_entry_count} "
            f"partials={len(partial_files)}"
        )
    finally:
        for partial_file in partial_files:
            partial_file.close()
        main_file.close()


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(
            "usage: replace_leb1.py MAIN_PATH OUTPUT_PATH PARTIAL_PATH [...]"
        )

    replace_selected_bodies(
        main_path=sys.argv[1],
        output_path=sys.argv[2],
        partial_paths=sys.argv[3:],
    )


if __name__ == "__main__":
    main()
PY
}

merge_leb1_for_tier() {
  local tier="$1"
  local main_path
  local group
  local partial_paths=()
  local command_args=()

  if ! ((ENABLE_MERGE)); then
    echo "LEB1 merge skipped (--no-merge)."
    return 0
  fi

  main_path="$(leb1_main_path "$tier")"

  if all_leb1_groups_selected; then
    for group in "${SELECTED_LEB1_GROUPS[@]}"; do
      partial_paths+=("$(leb1_partial_path "$tier" "$group")")
    done

    backup_file_if_present "$main_path"
    command_args=(
      scripts/generate_leb.py
      --tier "$tier"
      --merge "${partial_paths[@]}"
      --output "$main_path"
    )
    ((QUIET)) && command_args+=(--quiet)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
  else
    replace_leb1_groups_from_partials "$tier" || return 1
  fi
}

verify_leb1_for_tier() {
  local tier="$1"
  local group
  local target_path
  local command_args=()

  if ! ((ENABLE_VERIFICATION)); then
    echo "LEB1 verification skipped (--no-verify)."
    return 0
  fi

  if ((ENABLE_MERGE)); then
    target_path="$(leb1_main_path "$tier")"
    command_args=(
      scripts/generate_leb.py
      --tier "$tier"
      --verify-only
      --output "$target_path"
      --verify-samples "$LEB1_VERIFY_SAMPLES"
    )
    ((QUIET)) && command_args+=(--quiet)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
    return 0
  fi

  # Without a merge, each freshly generated partial is verified independently.
  for group in "${SELECTED_LEB1_GROUPS[@]}"; do
    target_path="$(leb1_partial_path "$tier" "$group")"
    command_args=(
      scripts/generate_leb.py
      --tier "$tier"
      --verify-only
      --output "$target_path"
      --verify-samples "$LEB1_VERIFY_SAMPLES"
    )
    ((QUIET)) && command_args+=(--quiet)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
  done
}

process_leb1_for_tier() {
  local tier="$1"

  print_header "LEB1 workflow" "tier: $tier"
  generate_leb1_partials_for_tier "$tier" || return 1
  merge_leb1_for_tier "$tier" || return 1
  verify_leb1_for_tier "$tier" || return 1
}

# ==============================================================================
# LEB2 conversion and verification
# ==============================================================================

assert_leb1_main_is_fresh() {
  local tier="$1"
  local main_path="$2"
  local newest_partial=""
  local partial_path

  [[ -f "$main_path" ]] \
    || fatal "cannot convert LEB2: missing LEB1 source: $main_path"

  # LEB2 is always converted from the merged main LEB1. A newer partial means
  # somebody regenerated groups with --no-merge, so converting now would use
  # stale data. Fail loudly instead of silently shipping inconsistent output.
  for partial_path in data/leb/ephemeris_"$tier"_*.leb; do
    [[ -e "$partial_path" ]] || continue

    if [[ -z "$newest_partial" || "$partial_path" -nt "$newest_partial" ]]; then
      newest_partial="$partial_path"
    fi
  done

  if [[ -n "$newest_partial" && "$newest_partial" -nt "$main_path" ]]; then
    fatal "cannot convert LEB2: main LEB1 '$main_path' is older than partial" \
      "'$newest_partial'; merge LEB1 first (remove --no-merge)."
  fi
}

convert_leb2_for_tier() {
  local tier="$1"
  local source_path
  local group
  local output_path
  local command_args=()

  print_header "LEB2 conversion" "tier: $tier"
  ensure_directory data/leb2

  source_path="$(leb1_main_path "$tier")"
  if ! ((DRY_RUN)); then
    assert_leb1_main_is_fresh "$tier" "$source_path"
  fi

  for group in "${SELECTED_LEB2_GROUPS[@]}"; do
    output_path="$(leb2_output_path "$tier" "$group")"
    backup_file_if_present "$output_path"

    # The converter uses both group and tier to authenticate the permitted body
    # inventory. Tier is semantic input, not merely part of the output name.
    command_args=(
      scripts/generate_leb2.py
      convert "$source_path"
      -o "$output_path"
      --group "$group"
      --tier "$tier"
    )
    ((QUIET)) && command_args+=(-q)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
  done
}

verify_leb2_for_tier() {
  local tier="$1"
  local reference_path
  local group
  local input_path
  local command_args=()

  if ! ((ENABLE_VERIFICATION)); then
    echo "LEB2 verification skipped (--no-verify)."
    return 0
  fi

  reference_path="$(leb1_main_path "$tier")"

  for group in "${SELECTED_LEB2_GROUPS[@]}"; do
    input_path="$(leb2_output_path "$tier" "$group")"
    command_args=(
      scripts/generate_leb2.py
      verify "$input_path"
      --reference "$reference_path"
      --samples "$LEB2_VERIFY_SAMPLES"
      --group "$group"
      --tier "$tier"
    )
    ((QUIET)) && command_args+=(-q)

    run_logged "$PYTHON" "${command_args[@]}" || return 1
  done
}

process_leb2_for_tier() {
  local tier="$1"

  convert_leb2_for_tier "$tier" || return 1
  verify_leb2_for_tier "$tier" || return 1
}

# ==============================================================================
# Installation
# ==============================================================================

install_file() {
  local source_path="$1"
  local destination_path="$2"
  local backup_directory="$INSTALL_DIR/backup-$RUN_ID"
  local backup_path="$backup_directory/$(basename "$destination_path")"

  if ! ((DRY_RUN)); then
    [[ -f "$source_path" ]] || fatal "install source not found: $source_path"
    mkdir -p "$INSTALL_DIR"
  fi

  if [[ -e "$destination_path" ]]; then
    if ((DRY_RUN)); then
      print_dry_run "cp -p $(printf '%q' "$destination_path") $(printf '%q' "$backup_path")"
    else
      mkdir -p "$backup_directory"
      cp -p "$destination_path" "$backup_path"
      print_note "backup $destination_path → $backup_path"
    fi
  fi

  if ((DRY_RUN)); then
    print_dry_run "cp -p $(printf '%q' "$source_path") $(printf '%q' "$destination_path")"
  else
    cp -p "$source_path" "$destination_path"
    printf '%s✓ installed%s  %s → %s\n' \
      "$C_GREEN" "$C_RESET" "$source_path" "$destination_path"
  fi
}

install_outputs_for_tier() {
  local tier="$1"
  local group
  local source_path

  ((ENABLE_INSTALL)) || return 0

  print_header "Install" "$tier → $INSTALL_DIR"

  if ((ENABLE_LEB1)) && ((ENABLE_MERGE)); then
    source_path="$(leb1_main_path "$tier")"
    install_file \
      "$source_path" \
      "$INSTALL_DIR/$(basename "$source_path")" \
      || return 1
  fi

  if ((ENABLE_LEB2)); then
    for group in "${SELECTED_LEB2_GROUPS[@]}"; do
      source_path="$(leb2_output_path "$tier" "$group")"
      install_file \
        "$source_path" \
        "$INSTALL_DIR/$(basename "$source_path")" \
        || return 1
    done
  fi
}

# ==============================================================================
# Configuration report
# ==============================================================================

print_configuration() {
  print_header "Configuration" "$(date '+%Y-%m-%d %H:%M:%S') · run $RUN_ID"
  printf '  %s%-13s%s %s\n' "$C_DIM" "Repo:" "$C_RESET" "$SCRIPT_DIR"
  printf '  %s%-13s%s %s\n' "$C_DIM" "Python:" "$C_RESET" "$PYTHON"
  printf '  %s%-13s%s %s\n' "$C_DIM" "Tiers:" "$C_RESET" "$(join_with_commas "${SELECTED_TIERS[@]}")"
  printf '  %s%-13s%s %s (%s)\n' "$C_DIM" "LEB1:" "$C_RESET" \
    "$(yes_or_no "$ENABLE_LEB1")" "$(join_with_commas "${SELECTED_LEB1_GROUPS[@]}")"
  printf '  %s%-13s%s %s (%s)\n' "$C_DIM" "LEB2:" "$C_RESET" \
    "$(yes_or_no "$ENABLE_LEB2")" "$(join_with_commas "${SELECTED_LEB2_GROUPS[@]}")"
  printf '  %s%-13s%s %s\n' "$C_DIM" "Merge:" "$C_RESET" "$(yes_or_no "$ENABLE_MERGE")"
  printf '  %s%-13s%s %s\n' "$C_DIM" "Verify:" "$C_RESET" "$(yes_or_no "$ENABLE_VERIFICATION")"

  if ((ENABLE_INSTALL)); then
    printf '  %s%-13s%s %s\n' "$C_DIM" "Install:" "$C_RESET" "$INSTALL_DIR"
  else
    printf '  %s%-13s%s %s\n' "$C_DIM" "Install:" "$C_RESET" "no"
  fi

  printf '  %s%-13s%s %s\n' "$C_DIM" "Dry run:" "$C_RESET" "$(yes_or_no "$DRY_RUN")"
  if ((SKIP_PREFLIGHT)); then
    printf '  %s%-13s%s skipped (--skip-deps-check)\n' "$C_DIM" "Preflight:" "$C_RESET"
  else
    printf '  %s%-13s%s %senabled%s\n' "$C_DIM" "Preflight:" "$C_RESET" "$C_GREEN" "$C_RESET"
  fi

  if nbody_required; then
    printf '  %s%-13s%s %syes%s (extended exotics: rebound/assist)\n' \
      "$C_DIM" "N-body req.:" "$C_RESET" "$C_YELLOW" "$C_RESET"
  else
    printf '  %s%-13s%s no\n' "$C_DIM" "N-body req.:" "$C_RESET"
  fi

  printf '  %s%-13s%s %s\n' "$C_DIM" "Log:" "$C_RESET" "$LOG_FILE"
}

# ==============================================================================
# Dependency and data preflight
# ==============================================================================

# Run dependency/data checks before acquiring the lock or touching outputs.
# In report mode (`--doctor`) failures are returned to the caller. During a real
# run they abort immediately. Dry-run also reports failures without aborting so
# an operator can still inspect the complete execution plan.
run_preflight() {
  local report_only=0
  local module_specification
  local extra_module
  local tiers_csv
  local flags_csv=""
  local exit_status

  [[ "${1:-}" == "--report" ]] && report_only=1

  module_specification="numpy:core,skyfield:core,erfa:core,click:core,zstandard:core,certifi:core,libephemeris:core"
  if nbody_required; then
    module_specification="$module_specification,rebound:nbody,assist:nbody"
  fi

  # Test/diagnostic hook: allow callers to probe additional import names.
  if [[ -n "${LEPH_DOCTOR_EXTRA_MODULES:-}" ]]; then
    for extra_module in ${LEPH_DOCTOR_EXTRA_MODULES//,/ }; do
      module_specification="$module_specification,${extra_module}:core"
    done
  fi

  tiers_csv="$(join_with_commas "${SELECTED_TIERS[@]}")"
  if nbody_required; then
    flags_csv="nbody"
  fi
  if spk_network_needed; then
    flags_csv="${flags_csv:+$flags_csv,}net"
  fi

  print_header "Preflight dependency check" \
    "tiers: $tiers_csv   groups: $(join_with_commas "${SELECTED_LEB1_GROUPS[@]}")"

  "$PYTHON" - "$module_specification" "$tiers_csv" "$flags_csv" <<'PY'
from __future__ import annotations

"""Validate runtime modules, kernels, ASSIST data, and network reachability."""

import importlib
import os
import sys
import urllib.request
from collections.abc import Iterable
from typing import Optional

# Colors mirror the Bash side: only on a TTY and disabled under NO_COLOR, so the
# report stays plain when redirected or piped.
_USE_COLOR = sys.stdout.isatty() and not os.environ.get("NO_COLOR")
if _USE_COLOR:
    _R = "\033[0m"
    _B = "\033[1m"
    _D = "\033[2m"
    _CY = "\033[36m"
    _G = "\033[32m"
    _Y = "\033[33m"
    _RED = "\033[31m"
else:
    _R = _B = _D = _CY = _G = _Y = _RED = ""

MODULE_SPECIFICATION, TIERS_CSV, FLAGS_CSV = sys.argv[1:4]
SELECTED_TIERS = [tier for tier in TIERS_CSV.split(",") if tier]
ENABLED_FLAGS = {flag for flag in FLAGS_CSV.split(",") if flag}

hard_failure_count = 0
warning_count = 0
required_extras: list[str] = []
assist_data_missing = False
printed_fix_commands: set[str] = set()


def require_extra(extra_name: Optional[str]) -> None:
    """Remember an installation extra without printing duplicates."""

    if extra_name and extra_name not in required_extras:
        required_extras.append(extra_name)


def report_ok(label: str) -> None:
    print(f"  {_G}[OK ]{_R} {label}")


def report_missing(label: str, extra_name: Optional[str] = None) -> None:
    global hard_failure_count

    print(f"  {_B}{_RED}[MISS]{_R} {label}")
    hard_failure_count += 1
    require_extra(extra_name)


def report_warning(label: str, note: str = "") -> None:
    global warning_count

    suffix = f"  -- {note}" if note else ""
    print(f"  {_Y}[WARN]{_R} {label}{suffix}")
    warning_count += 1


def check_python_modules(module_specification: str) -> None:
    print("Python modules:")

    for token in (item for item in module_specification.split(",") if item):
        module_name, _, extra_name = token.partition(":")
        extra_name = extra_name or "core"

        try:
            importlib.import_module(module_name)
        except Exception:
            report_missing(f"{module_name}  ({extra_name} extra)", extra_name)
        else:
            report_ok(module_name)


def check_vendored_spk_reader() -> None:
    try:
        importlib.import_module("libephemeris.vendor.spktype21")
    except Exception:
        report_missing("libephemeris.vendor.spktype21 (vendored)", "core")
    else:
        report_ok("libephemeris.vendor.spktype21 (vendored)")


def check_assist_data() -> None:
    """Check files that cannot be downloaded automatically during generation."""

    global assist_data_missing

    if "nbody" not in ENABLED_FLAGS:
        return

    print("\nASSIST data files (required for extended exotics N-body):")

    try:
        from libephemeris.rebound_integration import (
            _ASSIST_DEFAULT_DIR,
            AssistEphemConfig,
        )

        config = AssistEphemConfig()
        required_files = (
            ("planet ephemeris", config.planets_file),
            ("asteroid perturbers", config.asteroids_file),
            # Extended (-5000..+5000) exceeds DE440's 1550..2650 window, so
            # generation resolves to this exact deep-time DE441 file.
            (
                "deep-time DE441",
                str(_ASSIST_DEFAULT_DIR / "linux_m13000p17000.441"),
            ),
        )

        for label, path in required_files:
            if path and os.path.exists(path):
                report_ok(f"{label}: {os.path.basename(path)}")
                continue

            display_name = os.path.basename(path) if path else label
            report_missing(f"{label}: {display_name}", "nbody")
            assist_data_missing = True

    except Exception as error:
        report_missing(
            "ASSIST config resolution "
            f"({type(error).__name__}: {error})",
            "nbody",
        )
        assist_data_missing = True


def check_planetary_kernels(tiers: Iterable[str]) -> None:
    if not tiers:
        return

    print("\nJPL planetary kernels (auto-download if missing):")

    try:
        from libephemeris.state import _get_data_dir

        data_directory = _get_data_dir()
        kernel_by_tier = {
            "base": "de440s.bsp",
            "medium": "de440.bsp",
            "extended": "de441.bsp",
        }

        for tier in tiers:
            filename = kernel_by_tier.get(tier, f"<unknown tier {tier}>")
            path = os.path.join(data_directory, filename)

            if os.path.exists(path):
                report_ok(f"{filename} ({tier})")
            else:
                report_warning(
                    f"{filename} ({tier})",
                    f"absent from {data_directory}; skyfield will download",
                )

    except Exception as error:
        report_warning(
            f"kernel directory resolution ({type(error).__name__}: {error})"
        )


def check_horizons_network() -> None:
    """Warn, rather than fail, because cached kernels may be sufficient."""

    if "net" not in ENABLED_FLAGS:
        return

    print("\nNetwork (JPL Horizons reachability):")
    url = "https://ssd.jpl.nasa.gov/"

    try:
        with urllib.request.urlopen(url, timeout=5) as response:
            report_ok(f"reachable ({response.status} {url})")
    except Exception as error:
        report_warning(
            f"unreachable ({type(error).__name__})",
            "SPK auto-download/Horizons may fail",
        )


def print_fix_command(command: str, note: str = "") -> None:
    if command in printed_fix_commands:
        return

    printed_fix_commands.add(command)
    suffix = f"   {_D}# {note}{_R}" if note else ""
    print(f"  {_D}$ {_R}{command}{suffix}")


def print_remediation() -> None:
    print(
        f"{_B}{_Y}To fix, run the commands below, then re-run "
        f"'./regenerate-leb.sh --doctor':{_R}\n"
    )

    if "core" in required_extras:
        print_fix_command('uv pip install -e ".[dev]"', "core runtime deps")

    if "nbody" in required_extras or assist_data_missing:
        # assist's C extension needs rebound.h during compilation. Installing
        # rebound 4.x first and disabling build isolation makes that header
        # visible; a naive isolated install of .[nbody] can fail.
        print_fix_command(
            'uv pip install "rebound>=4.4.11,<5.0.0"',
            "rebound 5.x drops the C header assist needs",
        )
        print_fix_command(
            "uv pip install setuptools wheel",
            "--no-build-isolation needs these present",
        )
        print_fix_command(
            "uv pip install assist --no-build-isolation",
            "builds against rebound.h from rebound 4.x",
        )

    if assist_data_missing:
        print_fix_command(
            'uv run python -c "from libephemeris.rebound_integration import '
            'download_assist_data as download; '
            'download(planets_de441=True)"',
            "fetch ASSIST kernels (~2.6 GB including deep-time DE441)",
        )


def main() -> int:
    check_python_modules(MODULE_SPECIFICATION)
    check_vendored_spk_reader()
    check_assist_data()
    check_planetary_kernels(SELECTED_TIERS)
    check_horizons_network()

    print()
    if hard_failure_count == 0:
        print(
            f"{_G}RESULT: all required items present "
            f"({warning_count} warning(s)){_R}"
        )
        return 0

    print_remediation()
    print()
    print(
        f"{_B}{_RED}RESULT: {hard_failure_count} required item(s) MISSING, "
        f"{warning_count} warning(s){_R}"
    )
    return 2


raise SystemExit(main())
PY

  exit_status=$?
  if ((exit_status == 0)); then
    return 0
  fi

  if ((report_only)) || ((DRY_RUN)); then
    echo
    echo "Preflight found issues; continuing only because this is doctor/dry-run mode."
    return "$exit_status"
  fi

  echo
  fatal "preflight dependency check failed. Install the missing items and rerun." \
    "Use './regenerate-leb.sh --doctor' to re-check without regenerating."
}

# ==============================================================================
# Final result and main workflow
# ==============================================================================

print_final_result() {
  local index

  echo
  if ! has_failures; then
    printf '%s✓ SUCCESS%s  LEB regeneration completed\n' "$C_BOLD$C_GREEN" "$C_RESET"
    print_rule
    printf '%s  log: %s%s\n' "$C_DIM" "$LOG_FILE" "$C_RESET"
    return 0
  fi

  printf '%s✗ FAILED%s  one or more steps did not complete:\n' "$C_BOLD$C_RED" "$C_RESET"
  for index in "${!FAILED_COMMANDS[@]}"; do
    printf '%s  • %s%s %s(exit %d)%s\n' \
      "$C_RED" "${FAILED_COMMANDS[$index]}" "$C_RESET" "$C_DIM" "${FAILED_STATUSES[$index]}" "$C_RESET"
  done
  print_rule
  printf '%s  log: %s%s\n' "$C_DIM" "$LOG_FILE" "$C_RESET"
  return 1
}

main() {
  local tier

  parse_arguments "$@"
  resolve_configuration

  mkdir -p "$LOG_DIR"
  trap cleanup_lock EXIT

  # Environment validation happens before the lock and before any output backup.
  select_python_interpreter

  if ((DOCTOR_ONLY)); then
    run_preflight --report
    exit $?
  fi

  if ! ((SKIP_PREFLIGHT)); then
    # In dry-run mode a failing preflight returns non-zero but intentionally does
    # not stop the rest of the plan from being printed.
    run_preflight
  fi

  acquire_lock
  print_configuration

  refresh_refitted_tno_spk_cache | tee -a "$LOG_FILE"

  for tier in "${SELECTED_TIERS[@]}"; do
    if ((ENABLE_LEB1)); then
      process_leb1_for_tier "$tier"
      has_failures && break
    fi

    if ((ENABLE_LEB2)); then
      process_leb2_for_tier "$tier"
      has_failures && break
    fi

    install_outputs_for_tier "$tier"
    has_failures && break
  done

  print_final_result
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
