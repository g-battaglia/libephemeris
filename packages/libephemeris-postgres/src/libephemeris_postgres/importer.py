# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Resumable upload of pinned LEB2 files, followed by streaming verification."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from libephemeris import CoefficientSourceError
from libephemeris.download import DATA_FILES

from .config import admin_dsn
from .source import BLOCK_SIZE


def upload_files(paths: list[str | Path], *, dsn: str | None = None) -> None:
    """Upload missing blocks; publish only after the stored SHA-256 matches."""
    import psycopg

    for raw_path in paths:
        path = Path(raw_path)
        sha = DATA_FILES.get(path.name, {}).get("sha256")
        if not path.name.endswith(".leb2") or not sha:
            raise CoefficientSourceError("Upload requires a pinned LEB2 artifact")
        with path.open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != sha:
                raise CoefficientSourceError(
                    "Local LEB2 checksum does not match its pin"
                )
        size = path.stat().st_size
        try:
            with psycopg.connect(admin_dsn(dsn)) as conn:
                # One session lock spans batch commits and concurrent resumptions.
                conn.execute("SELECT pg_advisory_lock(hashtextextended(%s, 0))", (sha,))
                conn.execute(
                    "INSERT INTO libephemeris.files (sha256,name,size) VALUES (%s,%s,%s) "
                    "ON CONFLICT DO NOTHING",
                    (sha, path.name, size),
                )
                row = conn.execute(
                    "SELECT name,size,complete FROM libephemeris.files WHERE sha256=%s",
                    (sha,),
                ).fetchone()
                if row is None or row[:2] != (path.name, size):
                    raise CoefficientSourceError("Stored LEB2 identity mismatch")
                if row[2]:
                    continue
                present = {
                    row[0]
                    for row in conn.execute(
                        "SELECT block_no FROM libephemeris.blocks WHERE sha256=%s",
                        (sha,),
                    )
                }
                conn.commit()
                with path.open("rb") as stream:
                    for first in range(0, (size + BLOCK_SIZE - 1) // BLOCK_SIZE, 256):
                        with conn.cursor().copy(
                            "COPY libephemeris.blocks (sha256,block_no,data) "
                            "FROM STDIN (FORMAT BINARY)"
                        ) as copy:
                            copy.set_types(["text", "int4", "bytea"])
                            for number in range(
                                first,
                                min(first + 256, (size + BLOCK_SIZE - 1) // BLOCK_SIZE),
                            ):
                                data = stream.read(BLOCK_SIZE)
                                if number not in present:
                                    copy.write_row((sha, number, data))
                        conn.commit()
                _verify(conn, sha, size)
                conn.execute(
                    "UPDATE libephemeris.files SET complete=true WHERE sha256=%s",
                    (sha,),
                )
                conn.commit()
                conn.execute("ANALYZE libephemeris.files, libephemeris.blocks")
        except CoefficientSourceError:
            raise
        except Exception:
            raise CoefficientSourceError("PostgreSQL upload failed") from None


def _verify(conn: Any, sha: str, size: int) -> None:
    digest = hashlib.sha256()
    total = 0
    with conn.cursor(name="verify_blocks", binary=True) as cur:
        cur.execute(
            "SELECT block_no,data FROM libephemeris.blocks WHERE sha256=%s ORDER BY block_no",
            (sha,),
        )
        for expected, (number, data) in enumerate(cur):
            if number != expected or len(data) != min(BLOCK_SIZE, size - total):
                raise CoefficientSourceError("Stored LEB2 block layout mismatch")
            digest.update(data)
            total += len(data)
    if total != size or digest.hexdigest() != sha:
        raise CoefficientSourceError("Stored LEB2 checksum mismatch")
