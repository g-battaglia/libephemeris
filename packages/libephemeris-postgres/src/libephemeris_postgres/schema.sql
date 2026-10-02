-- SPDX-License-Identifier: AGPL-3.0-only
CREATE SCHEMA IF NOT EXISTS libephemeris;
CREATE TABLE IF NOT EXISTS libephemeris.schema_version (
  version integer PRIMARY KEY CHECK (version = 1)
);
INSERT INTO libephemeris.schema_version VALUES (1) ON CONFLICT DO NOTHING;

CREATE TABLE IF NOT EXISTS libephemeris.datasets (
  dataset_id uuid PRIMARY KEY,
  tier text NOT NULL CHECK (tier IN ('base', 'medium', 'extended')),
  complete boolean NOT NULL DEFAULT false,
  libephemeris_version text NOT NULL,
  created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE IF NOT EXISTS libephemeris.artifacts (
  dataset_id uuid NOT NULL REFERENCES libephemeris.datasets,
  artifact_no smallint NOT NULL,
  name text NOT NULL,
  sha256 text NOT NULL,
  jd_start float8 NOT NULL,
  jd_end float8 NOT NULL,
  delta_t_jd float8[] NOT NULL,
  delta_t_val float8[] NOT NULL,
  PRIMARY KEY (dataset_id, artifact_no), UNIQUE (dataset_id, name)
);
CREATE TABLE IF NOT EXISTS libephemeris.series (
  dataset_id uuid NOT NULL,
  artifact_no smallint NOT NULL,
  body_id integer NOT NULL,
  coord_type smallint NOT NULL,
  segment_count integer NOT NULL CHECK (segment_count > 0),
  jd_start float8 NOT NULL,
  jd_end float8 NOT NULL CHECK (jd_end > jd_start),
  interval_days float8 NOT NULL CHECK (interval_days > 0),
  degree smallint NOT NULL CHECK (degree BETWEEN 0 AND 256),
  components smallint NOT NULL CHECK (components > 0),
  page_size integer NOT NULL CHECK (page_size > 0),
  PRIMARY KEY (dataset_id, artifact_no, body_id),
  FOREIGN KEY (dataset_id, artifact_no) REFERENCES libephemeris.artifacts
);
CREATE TABLE IF NOT EXISTS libephemeris.pages (
  dataset_id uuid NOT NULL,
  artifact_no smallint NOT NULL,
  body_id integer NOT NULL,
  page_no integer NOT NULL CHECK (page_no >= 0),
  coeffs bytea NOT NULL,
  PRIMARY KEY (dataset_id, artifact_no, body_id, page_no),
  FOREIGN KEY (dataset_id, artifact_no, body_id)
    REFERENCES libephemeris.series
);
CREATE TABLE IF NOT EXISTS libephemeris.stars (
  dataset_id uuid NOT NULL,
  artifact_no smallint NOT NULL,
  star_id integer NOT NULL,
  ra_j2000 float8 NOT NULL,
  dec_j2000 float8 NOT NULL,
  pm_ra float8 NOT NULL,
  pm_dec float8 NOT NULL,
  parallax float8 NOT NULL,
  rv float8 NOT NULL,
  magnitude float8 NOT NULL,
  PRIMARY KEY (dataset_id, artifact_no, star_id)
);
