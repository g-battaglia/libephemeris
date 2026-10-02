-- SPDX-License-Identifier: AGPL-3.0-only
-- Project-native, language-neutral coefficient storage. Provision explicitly.
CREATE SCHEMA IF NOT EXISTS libephemeris;
CREATE TABLE IF NOT EXISTS libephemeris.schema_version (
    version integer PRIMARY KEY CHECK (version = 1)
);
INSERT INTO libephemeris.schema_version VALUES (1) ON CONFLICT DO NOTHING;
CREATE TABLE IF NOT EXISTS libephemeris.datasets (
    dataset_id uuid PRIMARY KEY,
    tier text NOT NULL CHECK (tier IN ('base', 'medium', 'extended')),
    published boolean NOT NULL DEFAULT false,
    manifest jsonb NOT NULL
);
CREATE TABLE IF NOT EXISTS libephemeris.series (
    dataset_id uuid NOT NULL REFERENCES libephemeris.datasets,
    body_id integer NOT NULL,
    coord_type integer NOT NULL CHECK (coord_type BETWEEN 0 AND 4),
    segment_count integer NOT NULL CHECK (segment_count > 0),
    jd_start double precision NOT NULL CHECK (jd_start > '-Infinity'::float8 AND jd_start < 'Infinity'::float8),
    jd_end double precision NOT NULL CHECK (jd_end >= jd_start AND jd_end < 'Infinity'::float8),
    interval_days double precision NOT NULL CHECK (interval_days > 0 AND interval_days < 'Infinity'::float8),
    degree integer NOT NULL CHECK (degree BETWEEN 0 AND 256),
    components integer NOT NULL CHECK (
        (body_id = -1 AND components = 2) OR
        (body_id <> -1 AND components = 3)
    ),
    PRIMARY KEY (dataset_id, body_id)
);
CREATE TABLE IF NOT EXISTS libephemeris.segments (
    dataset_id uuid NOT NULL,
    body_id integer NOT NULL,
    segment_index integer NOT NULL CHECK (segment_index >= 0),
    payload bytea NOT NULL,
    PRIMARY KEY (dataset_id, body_id, segment_index),
    FOREIGN KEY (dataset_id, body_id) REFERENCES libephemeris.series
);
CREATE TABLE IF NOT EXISTS libephemeris.delta_t (
    dataset_id uuid NOT NULL REFERENCES libephemeris.datasets,
    jd double precision NOT NULL,
    days double precision NOT NULL,
    PRIMARY KEY (dataset_id, jd)
);
CREATE TABLE IF NOT EXISTS libephemeris.stars (
    dataset_id uuid NOT NULL REFERENCES libephemeris.datasets,
    star_id integer NOT NULL,
    ra double precision NOT NULL,
    dec double precision NOT NULL,
    pm_ra double precision NOT NULL,
    pm_dec double precision NOT NULL,
    parallax double precision NOT NULL,
    rv double precision NOT NULL,
    magnitude double precision NOT NULL,
    PRIMARY KEY (dataset_id, star_id)
);
CREATE TABLE IF NOT EXISTS libephemeris.sections (
    dataset_id uuid NOT NULL REFERENCES libephemeris.datasets,
    artifact_index integer NOT NULL,
    section_id integer NOT NULL,
    payload bytea NOT NULL,
    PRIMARY KEY (dataset_id, artifact_index, section_id)
);
-- Published versions are immutable, including their scientific children.
CREATE OR REPLACE FUNCTION libephemeris.guard_publication() RETURNS trigger
LANGUAGE plpgsql AS $$
DECLARE id uuid; is_published boolean;
BEGIN
    IF TG_TABLE_NAME = 'datasets' THEN
        IF TG_OP = 'INSERT' THEN
            IF NEW.published THEN
                RAISE EXCEPTION 'Construct a dataset before publishing it';
            END IF;
        ELSIF OLD.published THEN
            RAISE EXCEPTION 'Published ephemeris datasets are immutable';
        END IF;
        -- Publication verifies the complete scientific inventory in one scan.
        -- Together with the parent-row locks below this prevents a concurrent
        -- writer from adding/removing segments while publication is checked.
        IF TG_OP = 'UPDATE' AND NEW.published THEN
            IF NOT EXISTS (SELECT 1 FROM libephemeris.series WHERE dataset_id = NEW.dataset_id AND body_id <> -1) THEN
                RAISE EXCEPTION 'Dataset contains no body series';
            END IF;
            IF EXISTS (
                SELECT 1 FROM libephemeris.series s
                LEFT JOIN libephemeris.segments p USING (dataset_id, body_id)
                WHERE s.dataset_id = NEW.dataset_id
                GROUP BY s.body_id, s.segment_count, s.degree, s.components
                HAVING count(p.segment_index) <> s.segment_count
                    OR min(p.segment_index) <> 0
                    OR max(p.segment_index) <> s.segment_count - 1
                    OR bool_or(octet_length(p.payload) <> (s.degree + 1) * s.components * 8)
            ) THEN
                RAISE EXCEPTION 'Coefficient inventory is incomplete or malformed';
            END IF;
        END IF;
    ELSE
        IF TG_OP = 'DELETE' THEN id := OLD.dataset_id; ELSE id := NEW.dataset_id; END IF;
        SELECT published INTO is_published FROM libephemeris.datasets
            WHERE dataset_id = id FOR SHARE;
        IF is_published THEN
            RAISE EXCEPTION 'Published ephemeris datasets are immutable';
        END IF;
        IF TG_OP = 'UPDATE' AND OLD.dataset_id <> NEW.dataset_id THEN
            RAISE EXCEPTION 'Dataset identity cannot change';
        END IF;
    END IF;
    IF TG_OP = 'DELETE' THEN RETURN OLD; ELSE RETURN NEW; END IF;
END $$;
DO $$
DECLARE tab text;
BEGIN
    FOREACH tab IN ARRAY ARRAY['datasets','series','segments','delta_t','stars','sections']
    LOOP
        IF NOT EXISTS (SELECT 1 FROM pg_trigger WHERE tgname = 'immutable_' || tab
            AND tgrelid = ('libephemeris.' || tab)::regclass) THEN
            EXECUTE format('CREATE TRIGGER %I BEFORE INSERT OR UPDATE OR DELETE ON libephemeris.%I FOR EACH ROW EXECUTE FUNCTION libephemeris.guard_publication()', 'immutable_' || tab, tab);
        END IF;
    END LOOP;
END $$;
