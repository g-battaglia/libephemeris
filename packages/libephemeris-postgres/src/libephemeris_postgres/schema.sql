-- SPDX-License-Identifier: AGPL-3.0-only
CREATE SCHEMA IF NOT EXISTS libephemeris;
CREATE TABLE IF NOT EXISTS libephemeris.files (
  sha256 text PRIMARY KEY CHECK (sha256 ~ '^[0-9a-f]{64}$'),
  name text NOT NULL,
  size bigint NOT NULL CHECK (size > 0),
  complete boolean NOT NULL DEFAULT false
);
CREATE TABLE IF NOT EXISTS libephemeris.blocks (
  sha256 text NOT NULL REFERENCES libephemeris.files,
  block_no integer NOT NULL CHECK (block_no >= 0),
  data bytea NOT NULL CHECK (octet_length(data) BETWEEN 1 AND 65536),
  PRIMARY KEY (sha256, block_no)
);
ALTER TABLE libephemeris.blocks ALTER COLUMN data SET STORAGE EXTERNAL;
