"""
Tests for SPK body name mapping.

These tests verify:
- SPK_BODY_NAME_MAP contains correct mappings for all supported minor bodies
- Helper functions correctly retrieve Horizons IDs and NAIF IDs
- Mapping values are consistent with NAIF ID constants
"""

import pytest

import libephemeris as eph
from libephemeris.constants import (
    # Body IDs
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    VARUNA,
    HYGIEA,
    DAVIDA,
    INTERAMNIA,
    BENNU,
    # NAIF IDs
    NAIF_CHIRON,
    NAIF_PHOLUS,
    NAIF_CERES,
    NAIF_PALLAS,
    NAIF_JUNO,
    NAIF_VESTA,
    NAIF_ERIS,
    NAIF_SEDNA,
    NAIF_HAUMEA,
    NAIF_MAKEMAKE,
    NAIF_IXION,
    NAIF_ORCUS,
    NAIF_QUAOAR,
    NAIF_ASTEROID_OFFSET,
    # Mapping and functions
    SPK_BODY_NAME_MAP,
    SPK_AUTO_DOWNLOAD_BLOCKED,
    REQUIRED_SPK_BODIES,
    get_horizons_id,
    get_naif_id_from_ipl,
    get_spk_body_info_from_map,
    is_spk_auto_download_blocked,
)


class TestSpkBodyNameMap:
    """Test the SPK_BODY_NAME_MAP dictionary."""

    def test_map_exists_and_not_empty(self):
        """SPK_BODY_NAME_MAP exists and contains entries."""
        assert SPK_BODY_NAME_MAP is not None
        assert len(SPK_BODY_NAME_MAP) > 0

    def test_map_contains_centaurs(self):
        """Map contains centaurs (Chiron, Pholus)."""
        assert CHIRON in SPK_BODY_NAME_MAP
        assert PHOLUS in SPK_BODY_NAME_MAP

    def test_map_contains_main_belt_asteroids(self):
        """Map contains main belt asteroids (Ceres, Pallas, Juno, Vesta)."""
        assert CERES in SPK_BODY_NAME_MAP
        assert PALLAS in SPK_BODY_NAME_MAP
        assert JUNO in SPK_BODY_NAME_MAP
        assert VESTA in SPK_BODY_NAME_MAP

    def test_map_contains_tnos(self):
        """Map contains trans-Neptunian objects."""
        assert ERIS in SPK_BODY_NAME_MAP
        assert SEDNA in SPK_BODY_NAME_MAP
        assert HAUMEA in SPK_BODY_NAME_MAP
        assert MAKEMAKE in SPK_BODY_NAME_MAP
        assert IXION in SPK_BODY_NAME_MAP
        assert ORCUS in SPK_BODY_NAME_MAP
        assert QUAOAR in SPK_BODY_NAME_MAP
        assert VARUNA in SPK_BODY_NAME_MAP

    def test_map_values_are_tuples(self):
        """Map values are tuples of (horizons_id, naif_id)."""
        for ipl, value in SPK_BODY_NAME_MAP.items():
            assert isinstance(value, tuple), f"Value for {ipl} is not a tuple"
            assert len(value) == 2, f"Value for {ipl} does not have 2 elements"
            horizons_id, naif_id = value
            assert isinstance(horizons_id, str), (
                f"Horizons ID for {ipl} is not a string"
            )
            assert isinstance(naif_id, int), f"NAIF ID for {ipl} is not an integer"


class TestChironMapping:
    """Test Chiron-specific mapping values."""

    def test_chiron_horizons_id(self):
        """Chiron's Horizons ID is '2060'."""
        horizons_id, _ = SPK_BODY_NAME_MAP[CHIRON]
        assert horizons_id == "2060"

    def test_chiron_naif_id(self):
        """Chiron's NAIF ID matches NAIF_CHIRON constant."""
        _, naif_id = SPK_BODY_NAME_MAP[CHIRON]
        assert naif_id == NAIF_CHIRON
        assert naif_id == 2060 + NAIF_ASTEROID_OFFSET


class TestMainBeltAsteroidMappings:
    """Test main belt asteroid mapping values."""

    def test_ceres_mapping(self):
        """Ceres mapping uses name syntax to bypass JPL major body index."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[CERES]
        assert horizons_id == "Ceres;"
        assert naif_id == NAIF_CERES
        assert naif_id == 1 + NAIF_ASTEROID_OFFSET

    def test_pallas_mapping(self):
        """Pallas mapping uses name syntax to bypass JPL major body index."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[PALLAS]
        assert horizons_id == "Pallas;"
        assert naif_id == NAIF_PALLAS

    def test_juno_mapping(self):
        """Juno mapping uses name syntax to bypass JPL major body index."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[JUNO]
        assert horizons_id == "Juno;"
        assert naif_id == NAIF_JUNO

    def test_vesta_mapping(self):
        """Vesta mapping uses name syntax to bypass JPL major body index."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[VESTA]
        assert horizons_id == "Vesta;"
        assert naif_id == NAIF_VESTA


class TestTnoMappings:
    """Test trans-Neptunian object mapping values."""

    def test_eris_mapping(self):
        """Eris mapping is correct (asteroid #136199)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[ERIS]
        assert horizons_id == "136199"
        assert naif_id == NAIF_ERIS
        assert naif_id == 136199 + NAIF_ASTEROID_OFFSET

    def test_sedna_mapping(self):
        """Sedna mapping is correct (asteroid #90377)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[SEDNA]
        assert horizons_id == "90377"
        assert naif_id == NAIF_SEDNA

    def test_haumea_mapping(self):
        """Haumea mapping is correct (asteroid #136108)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[HAUMEA]
        assert horizons_id == "136108"
        assert naif_id == NAIF_HAUMEA

    def test_makemake_mapping(self):
        """Makemake mapping is correct (asteroid #136472)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[MAKEMAKE]
        assert horizons_id == "136472"
        assert naif_id == NAIF_MAKEMAKE

    def test_quaoar_mapping(self):
        """Quaoar mapping is correct (asteroid #50000)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[QUAOAR]
        assert horizons_id == "50000"
        assert naif_id == NAIF_QUAOAR

    def test_varuna_mapping(self):
        """Varuna mapping is correct (asteroid #20000)."""
        horizons_id, naif_id = SPK_BODY_NAME_MAP[VARUNA]
        assert horizons_id == "20000"
        assert naif_id == 20000 + NAIF_ASTEROID_OFFSET


class TestGetHorizonsId:
    """Test the get_horizons_id() helper function."""

    def test_get_horizons_id_chiron(self):
        """get_horizons_id returns '2060' for Chiron."""
        assert get_horizons_id(CHIRON) == "2060"

    def test_get_horizons_id_eris(self):
        """get_horizons_id returns '136199' for Eris."""
        assert get_horizons_id(ERIS) == "136199"

    def test_get_horizons_id_ceres(self):
        """get_horizons_id returns 'Ceres;' for Ceres (name syntax)."""
        assert get_horizons_id(CERES) == "Ceres;"

    def test_get_horizons_id_unknown(self):
        """get_horizons_id returns None for unknown body ID."""
        assert get_horizons_id(999999) is None

    def test_get_horizons_id_all_mapped_bodies(self):
        """get_horizons_id works for all mapped bodies."""
        for ipl in SPK_BODY_NAME_MAP:
            result = get_horizons_id(ipl)
            assert result is not None, f"get_horizons_id returned None for ipl={ipl}"
            assert result == SPK_BODY_NAME_MAP[ipl][0]


class TestGetNaifIdFromIpl:
    """Test the get_naif_id_from_ipl() helper function."""

    def test_get_naif_id_chiron(self):
        """get_naif_id_from_ipl returns NAIF_CHIRON for Chiron."""
        assert get_naif_id_from_ipl(CHIRON) == NAIF_CHIRON

    def test_get_naif_id_eris(self):
        """get_naif_id_from_ipl returns NAIF_ERIS for Eris."""
        assert get_naif_id_from_ipl(ERIS) == NAIF_ERIS

    def test_get_naif_id_ceres(self):
        """get_naif_id_from_ipl returns NAIF_CERES for Ceres."""
        assert get_naif_id_from_ipl(CERES) == NAIF_CERES

    def test_get_naif_id_unknown(self):
        """get_naif_id_from_ipl returns None for unknown body ID."""
        assert get_naif_id_from_ipl(999999) is None

    def test_get_naif_id_all_mapped_bodies(self):
        """get_naif_id_from_ipl works for all mapped bodies."""
        for ipl in SPK_BODY_NAME_MAP:
            result = get_naif_id_from_ipl(ipl)
            assert result is not None, (
                f"get_naif_id_from_ipl returned None for ipl={ipl}"
            )
            assert result == SPK_BODY_NAME_MAP[ipl][1]


class TestGetSpkBodyInfoFromMap:
    """Test the get_spk_body_info_from_map() helper function."""

    def test_get_info_chiron(self):
        """get_spk_body_info_from_map returns correct tuple for Chiron."""
        result = get_spk_body_info_from_map(CHIRON)
        assert result == ("2060", NAIF_CHIRON)

    def test_get_info_eris(self):
        """get_spk_body_info_from_map returns correct tuple for Eris."""
        result = get_spk_body_info_from_map(ERIS)
        assert result == ("136199", NAIF_ERIS)

    def test_get_info_unknown(self):
        """get_spk_body_info_from_map returns None for unknown body ID."""
        assert get_spk_body_info_from_map(999999) is None

    def test_get_info_all_mapped_bodies(self):
        """get_spk_body_info_from_map works for all mapped bodies."""
        for ipl, expected in SPK_BODY_NAME_MAP.items():
            result = get_spk_body_info_from_map(ipl)
            assert result is not None, (
                f"get_spk_body_info_from_map returned None for ipl={ipl}"
            )
            assert result == expected


class TestNaifIdConsistency:
    """Test that NAIF IDs in the map match the NAIF_* constants."""

    @pytest.mark.parametrize(
        "ipl,naif_constant",
        [
            (CHIRON, NAIF_CHIRON),
            (PHOLUS, NAIF_PHOLUS),
            (CERES, NAIF_CERES),
            (PALLAS, NAIF_PALLAS),
            (JUNO, NAIF_JUNO),
            (VESTA, NAIF_VESTA),
            (ERIS, NAIF_ERIS),
            (SEDNA, NAIF_SEDNA),
            (HAUMEA, NAIF_HAUMEA),
            (MAKEMAKE, NAIF_MAKEMAKE),
            (IXION, NAIF_IXION),
            (ORCUS, NAIF_ORCUS),
            (QUAOAR, NAIF_QUAOAR),
        ],
    )
    def test_naif_id_matches_constant(self, ipl, naif_constant):
        """NAIF ID in map matches the corresponding NAIF_* constant."""
        _, naif_id = SPK_BODY_NAME_MAP[ipl]
        assert naif_id == naif_constant


class TestModuleExports:
    """Test that mapping and functions are exported from libephemeris."""

    def test_map_exported(self):
        """SPK_BODY_NAME_MAP is exported from libephemeris."""
        assert hasattr(eph, "SPK_BODY_NAME_MAP")
        assert eph.SPK_BODY_NAME_MAP is SPK_BODY_NAME_MAP

    def test_get_horizons_id_exported(self):
        """get_horizons_id is exported from libephemeris."""
        assert hasattr(eph, "get_horizons_id")
        assert eph.get_horizons_id(CHIRON) == "2060"

    def test_get_naif_id_from_ipl_exported(self):
        """get_naif_id_from_ipl is exported from libephemeris."""
        assert hasattr(eph, "get_naif_id_from_ipl")
        assert eph.get_naif_id_from_ipl(CHIRON) == NAIF_CHIRON

    def test_get_spk_body_info_from_map_exported(self):
        """get_spk_body_info_from_map is exported from libephemeris."""
        assert hasattr(eph, "get_spk_body_info_from_map")
        assert eph.get_spk_body_info_from_map(CHIRON) == ("2060", NAIF_CHIRON)


class TestHorizonsIdFormat:
    """Test that Horizons IDs are in the expected format for JPL Horizons API."""

    _MAJOR_BODY_INDEX_ASTEROIDS = {
        CERES,
        PALLAS,
        JUNO,
        VESTA,
        HYGIEA,
        DAVIDA,
        INTERAMNIA,
        BENNU,
    }

    def test_horizons_ids_are_numeric_strings(self):
        """All Horizons IDs are numeric strings, except major body index asteroids.

        Major body index asteroids (Ceres, Pallas, Juno, Vesta) use name syntax
        (e.g., "Ceres;") to bypass JPL Horizons restriction that refuses SPK
        generation for bodies in its major body index.
        """
        for ipl, (horizons_id, _) in SPK_BODY_NAME_MAP.items():
            if ipl in self._MAJOR_BODY_INDEX_ASTEROIDS:
                assert horizons_id.endswith(";"), (
                    f"Major body index asteroid {ipl} should use name syntax"
                )
            else:
                assert horizons_id.isdigit(), (
                    f"Horizons ID '{horizons_id}' for ipl={ipl} is not a numeric string"
                )

    def test_horizons_ids_match_asteroid_numbers(self):
        """Horizons IDs match expected asteroid catalog numbers or name syntax."""
        expected = {
            CHIRON: "2060",  # 2060 Chiron
            PHOLUS: "5145",  # 5145 Pholus
            CERES: "Ceres;",  # 1 Ceres (name syntax)
            PALLAS: "Pallas;",  # 2 Pallas (name syntax)
            JUNO: "Juno;",  # 3 Juno (name syntax)
            VESTA: "Vesta;",  # 4 Vesta (name syntax)
            ERIS: "136199",  # 136199 Eris
            SEDNA: "90377",  # 90377 Sedna
            HAUMEA: "136108",  # 136108 Haumea
            MAKEMAKE: "136472",  # 136472 Makemake
            IXION: "28978",  # 28978 Ixion
            ORCUS: "90482",  # 90482 Orcus
            QUAOAR: "50000",  # 50000 Quaoar
            VARUNA: "20000",  # 20000 Varuna
        }
        for ipl, expected_id in expected.items():
            horizons_id, _ = SPK_BODY_NAME_MAP[ipl]
            assert horizons_id == expected_id, (
                f"Horizons ID for ipl={ipl}: expected '{expected_id}', got '{horizons_id}'"
            )


class TestSpkAutoDownloadBlocked:
    """Verify the SPK_AUTO_DOWNLOAD_BLOCKED constant and helper."""

    def test_bennu_is_blocked(self):
        """Bennu is in the blocked set."""
        assert BENNU in SPK_AUTO_DOWNLOAD_BLOCKED
        assert is_spk_auto_download_blocked(BENNU)

    def test_blocked_bodies_are_in_name_map(self):
        """All blocked bodies must also be in SPK_BODY_NAME_MAP."""
        for ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            assert ipl in SPK_BODY_NAME_MAP, (
                f"Body {ipl} is in SPK_AUTO_DOWNLOAD_BLOCKED but not SPK_BODY_NAME_MAP"
            )

    def test_blocked_bodies_not_in_required(self):
        """Blocked bodies must not be in REQUIRED_SPK_BODIES."""
        for ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            assert ipl not in REQUIRED_SPK_BODIES, (
                f"Body {ipl} is blocked but also in REQUIRED_SPK_BODIES"
            )

    def test_is_spk_downloadable_returns_false_for_blocked(self):
        """is_spk_downloadable() returns False for blocked bodies."""
        from libephemeris.minor_bodies import is_spk_downloadable

        for ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            assert not is_spk_downloadable(ipl), (
                f"is_spk_downloadable should return False for blocked body {ipl}"
            )

    def test_get_horizons_id_still_works_for_blocked(self):
        """Blocked bodies still have valid Horizons IDs for manual use."""
        for ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            assert get_horizons_id(ipl) is not None

    def test_non_blocked_body_returns_false(self):
        """Non-blocked body returns False from helper."""
        assert not is_spk_auto_download_blocked(CHIRON)

    def test_exported_from_package(self):
        """SPK_AUTO_DOWNLOAD_BLOCKED is exported from libephemeris."""
        assert hasattr(eph, "SPK_AUTO_DOWNLOAD_BLOCKED")
        assert eph.SPK_AUTO_DOWNLOAD_BLOCKED is SPK_AUTO_DOWNLOAD_BLOCKED
        assert hasattr(eph, "is_spk_auto_download_blocked")
        assert eph.is_spk_auto_download_blocked is is_spk_auto_download_blocked
