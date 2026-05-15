"""
Tests for automatic SPK kernel download for major asteroids.

Tests verify:
- Major asteroid identification functions work correctly
- SPK availability checking works
- Auto-download function handles various scenarios gracefully
- Integration with calc_minor_body_heliocentric works
"""

import pytest
from unittest.mock import patch, MagicMock
import libephemeris as eph
from libephemeris import state
from libephemeris.constants import (
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    CHIRON,
    ERIS,
    SEDNA,
    NAIF_CERES,
    NAIF_CHIRON,
)
from libephemeris.minor_bodies import (
    is_major_asteroid,
    get_major_asteroid_info,
    list_major_asteroids,
    is_spk_available_for_body,
    auto_download_asteroid_spk,
    ensure_major_asteroid_spk,
    MAJOR_ASTEROID_SPK_INFO,
    calc_minor_body_heliocentric,
)


class TestIsMajorAsteroid:
    """Test is_major_asteroid() function."""

    def test_ceres_is_major(self):
        """Ceres should be identified as a major asteroid."""
        assert is_major_asteroid(CERES) is True

    def test_pallas_is_major(self):
        """Pallas should be identified as a major asteroid."""
        assert is_major_asteroid(PALLAS) is True

    def test_juno_is_major(self):
        """Juno should be identified as a major asteroid."""
        assert is_major_asteroid(JUNO) is True

    def test_vesta_is_major(self):
        """Vesta should be identified as a major asteroid."""
        assert is_major_asteroid(VESTA) is True

    def test_chiron_is_major(self):
        """Chiron should be identified as a major asteroid."""
        assert is_major_asteroid(CHIRON) is True

    def test_eris_is_not_major(self):
        """Eris (TNO) should not be identified as a major asteroid."""
        assert is_major_asteroid(ERIS) is False

    def test_sedna_is_not_major(self):
        """Sedna (TNO) should not be identified as a major asteroid."""
        assert is_major_asteroid(SEDNA) is False

    def test_invalid_id_is_not_major(self):
        """Invalid body ID should not be identified as a major asteroid."""
        assert is_major_asteroid(999999) is False
        assert is_major_asteroid(-1) is False


class TestGetMajorAsteroidInfo:
    """Test get_major_asteroid_info() function."""

    def test_ceres_info(self):
        """Should return correct info for Ceres."""
        info = get_major_asteroid_info(CERES)
        assert info is not None
        ast_num, horizons_id, naif_id, name = info
        assert ast_num == 1
        assert horizons_id == "Ceres;"
        assert naif_id == NAIF_CERES
        assert name == "Ceres"

    def test_chiron_info(self):
        """Should return correct info for Chiron."""
        info = get_major_asteroid_info(CHIRON)
        assert info is not None
        ast_num, horizons_id, naif_id, name = info
        assert ast_num == 2060
        assert horizons_id == "2060"
        assert naif_id == NAIF_CHIRON
        assert name == "Chiron"

    def test_eris_returns_none(self):
        """Eris (TNO) should return None."""
        assert get_major_asteroid_info(ERIS) is None

    def test_invalid_id_returns_none(self):
        """Invalid body ID should return None."""
        assert get_major_asteroid_info(999999) is None


class TestListMajorAsteroids:
    """Test list_major_asteroids() function."""

    def test_returns_list(self):
        """Should return a list."""
        result = list_major_asteroids()
        assert isinstance(result, list)

    def test_contains_expected_bodies(self):
        """Should contain all expected major asteroids."""
        result = list_major_asteroids()
        body_ids = [item[0] for item in result]
        names = [item[1] for item in result]

        assert CERES in body_ids
        assert PALLAS in body_ids
        assert JUNO in body_ids
        assert VESTA in body_ids
        assert CHIRON in body_ids

        assert "Ceres" in names
        assert "Pallas" in names
        assert "Juno" in names
        assert "Vesta" in names
        assert "Chiron" in names

    def test_correct_tuple_format(self):
        """Each item should be (body_id, name) tuple."""
        result = list_major_asteroids()
        for item in result:
            assert isinstance(item, tuple)
            assert len(item) == 2
            assert isinstance(item[0], int)
            assert isinstance(item[1], str)


class TestIsSpkAvailableForBody:
    """Test is_spk_available_for_body() function."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    def test_returns_false_when_not_registered(self):
        """Should return False when no SPK is registered."""
        assert is_spk_available_for_body(CERES) is False

    def test_returns_true_when_registered(self):
        """Should return True when SPK is registered."""
        state._SPK_BODY_MAP[CERES] = ("/path/to/ceres.bsp", NAIF_CERES)
        assert is_spk_available_for_body(CERES) is True


class TestAutoDownloadAsteroidSpk:
    """Test auto_download_asteroid_spk() function."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    def test_returns_none_for_non_major_asteroid(self):
        """Should return None for bodies not in any supported map."""
        result = auto_download_asteroid_spk(999999)
        assert result is None

    def test_returns_none_for_invalid_body(self):
        """Should return None for invalid body IDs."""
        result = auto_download_asteroid_spk(999999)
        assert result is None

    @patch("libephemeris.spk_auto._check_astroquery_available")
    def test_returns_path_when_already_registered(self, mock_check):
        """Should return path when SPK is already registered."""
        mock_check.return_value = True
        expected_path = "/path/to/ceres.bsp"
        state._SPK_BODY_MAP[CERES] = (expected_path, NAIF_CERES)

        result = auto_download_asteroid_spk(CERES)
        assert result == expected_path

    @pytest.mark.skip(reason="Implementation now uses direct HTTP, not astroquery")
    @patch("libephemeris.spk_auto._check_astroquery_available")
    def test_returns_none_when_astroquery_unavailable(self, mock_check):
        """Should return None when astroquery is not available."""
        mock_check.return_value = False
        result = auto_download_asteroid_spk(CERES)
        assert result is None

    @pytest.mark.skip(
        reason="Implementation now uses direct HTTP via download_and_register_spk"
    )
    @patch("libephemeris.spk_auto._check_astroquery_available")
    @patch("libephemeris.spk_auto.auto_get_spk")
    def test_calls_auto_get_spk_with_correct_params(self, mock_auto_get, mock_check):
        """Should call auto_get_spk with correct parameters."""
        mock_check.return_value = True
        mock_auto_get.return_value = "/path/to/ceres.bsp"

        result = auto_download_asteroid_spk(CERES)

        mock_auto_get.assert_called_once()
        call_kwargs = mock_auto_get.call_args.kwargs
        assert call_kwargs["body_id"] == "1"
        assert call_kwargs["ipl"] == CERES
        assert call_kwargs["naif_id"] == NAIF_CERES
        assert "jd_start" in call_kwargs
        assert "jd_end" in call_kwargs
        assert result == "/path/to/ceres.bsp"

    @pytest.mark.skip(
        reason="Implementation now uses direct HTTP via download_and_register_spk"
    )
    @patch("libephemeris.spk_auto._check_astroquery_available")
    @patch("libephemeris.spk_auto.auto_get_spk")
    def test_uses_custom_date_range(self, mock_auto_get, mock_check):
        """Should use custom date range when provided."""
        mock_check.return_value = True
        mock_auto_get.return_value = "/path/to/ceres.bsp"

        jd_start = 2451545.0
        jd_end = 2455197.5

        auto_download_asteroid_spk(CERES, jd_start=jd_start, jd_end=jd_end)

        call_kwargs = mock_auto_get.call_args.kwargs
        assert call_kwargs["jd_start"] == jd_start
        assert call_kwargs["jd_end"] == jd_end

    @pytest.mark.skip(
        reason="Implementation now uses direct HTTP via download_and_register_spk"
    )
    @patch("libephemeris.spk_auto._check_astroquery_available")
    @patch("libephemeris.spk_auto.auto_get_spk")
    def test_handles_download_exception_gracefully(self, mock_auto_get, mock_check):
        """Should return None when download fails."""
        mock_check.return_value = True
        mock_auto_get.side_effect = Exception("Network error")

        result = auto_download_asteroid_spk(CERES)
        assert result is None


class TestEnsureMajorAsteroidSpk:
    """Test ensure_major_asteroid_spk() function."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    def test_returns_true_when_already_available(self):
        """Should return True immediately when SPK is already registered."""
        state._SPK_BODY_MAP[CERES] = ("/path/to/ceres.bsp", NAIF_CERES)
        result = ensure_major_asteroid_spk(CERES)
        assert result is True

    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    def test_calls_download_when_not_available(self, mock_download):
        """Should call auto_download when SPK is not available."""
        mock_download.return_value = "/path/to/ceres.bsp"

        result = ensure_major_asteroid_spk(CERES)

        mock_download.assert_called_once()
        assert result is True

    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    def test_returns_false_when_download_fails(self, mock_download):
        """Should return False when download fails."""
        mock_download.return_value = None

        result = ensure_major_asteroid_spk(CERES)

        assert result is False


class TestMajorAsteroidSpkInfo:
    """Test MAJOR_ASTEROID_SPK_INFO dictionary."""

    def test_has_expected_keys(self):
        """Dictionary should have all expected body IDs."""
        assert CERES in MAJOR_ASTEROID_SPK_INFO
        assert PALLAS in MAJOR_ASTEROID_SPK_INFO
        assert JUNO in MAJOR_ASTEROID_SPK_INFO
        assert VESTA in MAJOR_ASTEROID_SPK_INFO
        assert CHIRON in MAJOR_ASTEROID_SPK_INFO

    def test_value_format(self):
        """Each value should be (asteroid_number, horizons_id, naif_id, name)."""
        for body_id, info in MAJOR_ASTEROID_SPK_INFO.items():
            assert isinstance(info, tuple)
            assert len(info) == 4
            ast_num, horizons_id, naif_id, name = info
            assert isinstance(ast_num, int)
            assert isinstance(horizons_id, str)
            assert isinstance(naif_id, int)
            assert isinstance(name, str)

    def test_naif_id_convention(self):
        """NAIF IDs should follow the convention: asteroid_number + 2000000."""
        for body_id, info in MAJOR_ASTEROID_SPK_INFO.items():
            ast_num, _, naif_id, _ = info
            expected_naif = ast_num + 2000000
            assert naif_id == expected_naif, (
                f"NAIF ID mismatch for body {body_id}: "
                f"expected {expected_naif}, got {naif_id}"
            )


class TestCalcMinorBodyHeliocentricWithSpk:
    """Test calc_minor_body_heliocentric with SPK integration."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    def test_falls_back_to_keplerian_when_no_spk(self):
        """Should use Keplerian when no SPK is registered."""
        jd = 2451545.0  # J2000.0
        lon, lat, dist = calc_minor_body_heliocentric(CERES, jd)

        # Just verify we get reasonable values
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_use_spk_false_bypasses_spk(self):
        """Should bypass SPK when use_spk=False."""
        # Even if SPK is registered, use_spk=False should use Keplerian
        state._SPK_BODY_MAP[CERES] = ("/path/to/ceres.bsp", NAIF_CERES)

        jd = 2451545.0
        # This should not raise an error even though the SPK path is fake
        # because use_spk=False bypasses SPK entirely
        lon, lat, dist = calc_minor_body_heliocentric(CERES, jd, use_spk=False)

        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_invalid_body_raises_error(self):
        """Should raise ValueError for invalid body ID."""
        with pytest.raises(ValueError, match="illegal planet number"):
            calc_minor_body_heliocentric(999999, 2451545.0)


class TestExportsFromMain:
    """Test that functions are exported from main libephemeris module."""

    def test_auto_download_asteroid_spk_exported(self):
        """auto_download_asteroid_spk should be exported."""
        assert hasattr(eph, "auto_download_asteroid_spk")
        assert callable(eph.auto_download_asteroid_spk)

    def test_is_spk_available_for_body_exported(self):
        """is_spk_available_for_body should be exported."""
        assert hasattr(eph, "is_spk_available_for_body")
        assert callable(eph.is_spk_available_for_body)

    def test_ensure_major_asteroid_spk_exported(self):
        """ensure_major_asteroid_spk should be exported."""
        assert hasattr(eph, "ensure_major_asteroid_spk")
        assert callable(eph.ensure_major_asteroid_spk)

    def test_is_major_asteroid_exported(self):
        """is_major_asteroid should be exported."""
        assert hasattr(eph, "is_major_asteroid")
        assert callable(eph.is_major_asteroid)

    def test_list_major_asteroids_exported(self):
        """list_major_asteroids should be exported."""
        assert hasattr(eph, "list_major_asteroids")
        assert callable(eph.list_major_asteroids)

    def test_major_asteroid_spk_info_exported(self):
        """MAJOR_ASTEROID_SPK_INFO should be exported."""
        assert hasattr(eph, "MAJOR_ASTEROID_SPK_INFO")
        assert isinstance(eph.MAJOR_ASTEROID_SPK_INFO, dict)


class TestEnsureMajorAsteroidSpkLogging:
    """Test logging behavior in ensure_major_asteroid_spk()."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    @pytest.fixture
    def enable_propagation(self):
        """Enable log propagation and DEBUG level temporarily so caplog captures logs."""
        import logging

        from libephemeris.logging_config import get_logger

        # Ensure logger is configured by calling get_logger()
        logger = get_logger()

        original_propagate = logger.propagate
        original_level = logger.level
        original_handler_levels = [h.level for h in logger.handlers]

        logger.propagate = True
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)

        yield

        logger.propagate = original_propagate
        logger.setLevel(original_level)
        for handler, level in zip(logger.handlers, original_handler_levels):
            handler.setLevel(level)

    def test_logs_debug_when_checking_availability(self, caplog, enable_propagation):
        """Should log DEBUG message when checking SPK availability."""
        import logging

        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            state._SPK_BODY_MAP[CERES] = ("/path/to/ceres.bsp", NAIF_CERES)
            ensure_major_asteroid_spk(CERES)

        # Check that checking message was logged
        assert any(
            "Checking SPK availability for body" in record.message
            and "Ceres" in record.message
            for record in caplog.records
        )

    def test_logs_debug_when_already_cached(self, caplog, enable_propagation):
        """Should log DEBUG message when SPK is already cached."""
        import logging

        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            state._SPK_BODY_MAP[CERES] = ("/path/to/ceres.bsp", NAIF_CERES)
            ensure_major_asteroid_spk(CERES)

        # Check that cache hit message was logged at DEBUG level
        cache_records = [
            record
            for record in caplog.records
            if "already cached" in record.message and "Ceres" in record.message
        ]
        assert len(cache_records) >= 1
        assert cache_records[0].levelno == logging.DEBUG

    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    def test_logs_info_when_downloading(
        self, mock_download, caplog, enable_propagation
    ):
        """Should log INFO message when download is needed."""
        import logging

        mock_download.return_value = "/path/to/ceres.bsp"

        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            ensure_major_asteroid_spk(CERES)

        # Check that download message was logged at INFO level
        download_records = [
            record
            for record in caplog.records
            if "not cached, downloading" in record.message and "Ceres" in record.message
        ]
        assert len(download_records) >= 1
        assert download_records[0].levelno == logging.INFO

    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    def test_logs_body_name_not_just_id(
        self, mock_download, caplog, enable_propagation
    ):
        """Log messages should include body name (e.g., 'Chiron'), not just ID."""
        import logging

        mock_download.return_value = "/path/to/chiron.bsp"

        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            ensure_major_asteroid_spk(CHIRON)

        # All log messages mentioning this body should include the name
        relevant_records = [
            record
            for record in caplog.records
            if "Chiron" in record.message or str(CHIRON) in record.message
        ]
        assert len(relevant_records) >= 1
        # Ensure "Chiron" appears (not just the numeric ID)
        assert any("Chiron" in record.message for record in relevant_records)


class TestRequiredSpkBodies:
    """Test REQUIRED_SPK_BODIES constant."""

    def test_is_frozenset(self):
        """REQUIRED_SPK_BODIES should be an immutable frozenset."""
        from libephemeris.constants import REQUIRED_SPK_BODIES

        assert isinstance(REQUIRED_SPK_BODIES, frozenset)

    def test_has_expected_bodies(self):
        """REQUIRED_SPK_BODIES should contain all expected body IDs."""
        from libephemeris.constants import (
            REQUIRED_SPK_BODIES,
            CHIRON,
            CERES,
            PALLAS,
            JUNO,
            VESTA,
            PHOLUS,
            NESSUS,
            ERIS,
        )

        # All major asteroids from MAJOR_ASTEROID_SPK_INFO
        assert CERES in REQUIRED_SPK_BODIES
        assert PALLAS in REQUIRED_SPK_BODIES
        assert JUNO in REQUIRED_SPK_BODIES
        assert VESTA in REQUIRED_SPK_BODIES
        assert CHIRON in REQUIRED_SPK_BODIES

        # Additional centaurs and dwarf planets
        assert PHOLUS in REQUIRED_SPK_BODIES
        assert NESSUS in REQUIRED_SPK_BODIES
        assert ERIS in REQUIRED_SPK_BODIES

    def test_all_bodies_in_spk_body_name_map(self):
        """All REQUIRED_SPK_BODIES should have entries in SPK_BODY_NAME_MAP."""
        from libephemeris.constants import REQUIRED_SPK_BODIES, SPK_BODY_NAME_MAP

        for body_id in REQUIRED_SPK_BODIES:
            assert body_id in SPK_BODY_NAME_MAP, (
                f"Body {body_id} is in REQUIRED_SPK_BODIES but not in SPK_BODY_NAME_MAP"
            )

    def test_exported_from_main_module(self):
        """REQUIRED_SPK_BODIES should be exported from main libephemeris module."""
        assert hasattr(eph, "REQUIRED_SPK_BODIES")
        assert isinstance(eph.REQUIRED_SPK_BODIES, frozenset)

    def test_count_is_expected(self):
        """REQUIRED_SPK_BODIES should have 8 bodies."""
        from libephemeris.constants import REQUIRED_SPK_BODIES

        assert len(REQUIRED_SPK_BODIES) == 8

    def test_major_asteroids_are_subset(self):
        """All bodies in MAJOR_ASTEROID_SPK_INFO should be in REQUIRED_SPK_BODIES."""
        from libephemeris.constants import REQUIRED_SPK_BODIES

        for body_id in MAJOR_ASTEROID_SPK_INFO.keys():
            assert body_id in REQUIRED_SPK_BODIES, (
                f"Body {body_id} is in MAJOR_ASTEROID_SPK_INFO but not in REQUIRED_SPK_BODIES"
            )


class TestEnsureMajorAsteroidSpkForNonMajor:
    """Test ensure_major_asteroid_spk with non-major asteroids (REQUIRED_SPK_BODIES extras)."""

    @pytest.fixture(autouse=True)
    def clear_spk_state(self):
        """Clear SPK state before and after each test."""
        orig_map = dict(state._SPK_BODY_MAP)
        state._SPK_BODY_MAP.clear()
        yield
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_map)

    def test_returns_true_when_already_cached_for_pholus(self):
        """Should return True when SPK is already cached for non-major asteroid."""
        from libephemeris.constants import PHOLUS, NAIF_PHOLUS

        state._SPK_BODY_MAP[PHOLUS] = ("/path/to/pholus.bsp", NAIF_PHOLUS)

        result = ensure_major_asteroid_spk(PHOLUS)

        assert result is True

    def test_returns_true_when_already_cached_for_eris(self):
        """Should return True when SPK is already cached for Eris."""
        from libephemeris.constants import ERIS, NAIF_ERIS

        state._SPK_BODY_MAP[ERIS] = ("/path/to/eris.bsp", NAIF_ERIS)

        result = ensure_major_asteroid_spk(ERIS)

        assert result is True

    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    def test_returns_false_when_astroquery_not_available_for_pholus(self, mock_auto):
        """Should return False when astroquery is not available."""
        from libephemeris.constants import PHOLUS

        # auto_download_asteroid_spk returns None for non-major asteroids
        mock_auto.return_value = None

        # Mock spk_auto._check_astroquery_available to return False
        with patch(
            "libephemeris.spk_auto._check_astroquery_available", return_value=False
        ):
            result = ensure_major_asteroid_spk(PHOLUS)

        assert result is False

    @pytest.mark.skip(
        reason="Old architecture: auto_download_asteroid_spk now handles all bodies"
    )
    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    @patch("libephemeris.spk_auto._check_astroquery_available")
    @patch("libephemeris.spk_auto.auto_get_spk")
    def test_uses_spk_body_name_map_for_pholus(
        self, mock_auto_get, mock_check, mock_auto_download
    ):
        """Should use SPK_BODY_NAME_MAP for non-major asteroids."""
        from libephemeris.constants import PHOLUS, NAIF_PHOLUS

        mock_auto_download.return_value = None
        mock_check.return_value = True
        mock_auto_get.return_value = "/path/to/pholus.bsp"

        result = ensure_major_asteroid_spk(PHOLUS)

        assert result is True
        mock_auto_get.assert_called_once()
        call_kwargs = mock_auto_get.call_args[1]
        assert call_kwargs["body_id"] == "5145"
        assert call_kwargs["ipl"] == PHOLUS
        assert call_kwargs["naif_id"] == NAIF_PHOLUS

    @pytest.mark.skip(
        reason="Old architecture: auto_download_asteroid_spk now handles all bodies"
    )
    @patch("libephemeris.minor_bodies.auto_download_asteroid_spk")
    @patch("libephemeris.spk_auto._check_astroquery_available")
    @patch("libephemeris.spk_auto.auto_get_spk")
    def test_uses_spk_body_name_map_for_eris(
        self, mock_auto_get, mock_check, mock_auto_download
    ):
        """Should use SPK_BODY_NAME_MAP for Eris."""
        from libephemeris.constants import ERIS, NAIF_ERIS

        mock_auto_download.return_value = None
        mock_check.return_value = True
        mock_auto_get.return_value = "/path/to/eris.bsp"

        result = ensure_major_asteroid_spk(ERIS)

        assert result is True
        mock_auto_get.assert_called_once()
        call_kwargs = mock_auto_get.call_args[1]
        assert call_kwargs["body_id"] == "136199"
        assert call_kwargs["ipl"] == ERIS
        assert call_kwargs["naif_id"] == NAIF_ERIS
