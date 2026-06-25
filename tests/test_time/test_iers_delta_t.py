"""
Tests for IERS observed Delta T download and lookup functionality.

This module tests the automatic download of observed Delta T values from IERS
for recent dates (1973-present).
"""

import pytest
import os

from libephemeris import julday
from libephemeris.iers_data import (
    download_delta_t_data,
    get_observed_delta_t,
    get_observed_delta_t_data_range,
    is_observed_delta_t_available,
    get_delta_t_iers,
    clear_iers_cache,
    set_iers_cache_dir,
    set_iers_auto_download,
    get_iers_cache_info,
    load_iers_data,
    _parse_delta_t_data,
    _get_delta_t_cache_path,
    DeltaTDataPoint,
)


class TestDeltaTDataPointDataclass:
    """Test the DeltaTDataPoint dataclass."""

    @pytest.mark.unit
    def test_delta_t_data_point_creation(self):
        """Test creating a DeltaTDataPoint instance."""
        point = DeltaTDataPoint(
            mjd=51544.0,
            year=2000,
            month=1,
            day=1,
            delta_t=63.83,
        )
        assert point.mjd == 51544.0
        assert point.year == 2000
        assert point.month == 1
        assert point.day == 1
        assert point.delta_t == 63.83


class TestParseDeltaTData:
    """Test parsing of IERS Delta T data file format."""

    @pytest.fixture
    def sample_delta_t_file(self, tmp_path):
        """Create a sample Delta T data file for testing."""
        content = """ 1973  2  1  43.4724
 1973  3  1  43.5648
 1973  4  1  43.6737
 2000  1  1  63.8285
 2000  2  1  63.8557
 2020  1  1  69.3612
 2020  6  1  69.4386
"""
        file_path = tmp_path / "deltat.data"
        file_path.write_text(content)
        return str(file_path)

    @pytest.mark.unit
    def test_parse_delta_t_data_valid_file(self, sample_delta_t_file):
        """Test parsing a valid Delta T data file."""
        data = _parse_delta_t_data(sample_delta_t_file)
        assert len(data) == 7
        assert all(isinstance(p, DeltaTDataPoint) for p in data)

    @pytest.mark.unit
    def test_parse_delta_t_data_sorted(self, sample_delta_t_file):
        """Test that parsed data is sorted by MJD."""
        data = _parse_delta_t_data(sample_delta_t_file)
        mjds = [p.mjd for p in data]
        assert mjds == sorted(mjds)

    @pytest.mark.unit
    def test_parse_delta_t_data_values(self, sample_delta_t_file):
        """Test specific parsed values."""
        data = _parse_delta_t_data(sample_delta_t_file)
        # Find the 2000-01-01 entry
        jan2000 = [p for p in data if p.year == 2000 and p.month == 1][0]
        assert jan2000.day == 1
        assert jan2000.delta_t == pytest.approx(63.8285, abs=0.0001)

    @pytest.mark.unit
    def test_parse_delta_t_data_skips_comments(self, tmp_path):
        """Test that comment lines are skipped."""
        content = """# This is a comment
 2000  1  1  63.8285
# Another comment
 2000  2  1  63.8557
"""
        file_path = tmp_path / "deltat.data"
        file_path.write_text(content)
        data = _parse_delta_t_data(str(file_path))
        assert len(data) == 2

    @pytest.mark.unit
    def test_parse_delta_t_data_skips_malformed_lines(self, tmp_path):
        """Test that malformed lines are skipped."""
        content = """ 2000  1  1  63.8285
malformed line
 2000  2  1
 2000  3  1  63.9075
invalid year month day delta_t
"""
        file_path = tmp_path / "deltat.data"
        file_path.write_text(content)
        data = _parse_delta_t_data(str(file_path))
        # Should only get the valid entries
        assert len(data) == 2


class TestDeltaTDownload:
    """Test Delta T data download functionality."""

    @pytest.fixture(autouse=True)
    def setup_and_cleanup(self, tmp_path):
        """Setup temporary cache directory and cleanup after tests."""
        # Store original state
        original_dir = None
        try:
            from libephemeris.iers_data import _IERS_CACHE_DIR

            original_dir = _IERS_CACHE_DIR
        except ImportError:
            pass

        # Clear any existing cache before test
        clear_iers_cache()

        # Set up temporary cache directory
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        yield

        # Cleanup
        clear_iers_cache()
        set_iers_cache_dir(original_dir)

    @pytest.mark.network
    def test_download_delta_t_data(self):
        """Test downloading Delta T data from IERS."""
        # Enable auto-download for this test
        set_iers_auto_download(True)
        try:
            path = download_delta_t_data(force=True)
            assert os.path.exists(path)
            assert path.endswith("deltat.data")
            # Check file has content
            assert os.path.getsize(path) > 1000  # Should be at least 1KB
        except ConnectionError:
            pytest.skip("Network not available")
        finally:
            set_iers_auto_download(False)

    @pytest.mark.network
    def test_download_delta_t_data_uses_cache(self):
        """Test that cached file is used when available."""
        set_iers_auto_download(True)
        try:
            # First download
            path1 = download_delta_t_data(force=True)
            mtime1 = os.path.getmtime(path1)

            # Second download should use cache
            path2 = download_delta_t_data(force=False)
            mtime2 = os.path.getmtime(path2)

            assert path1 == path2
            assert mtime1 == mtime2
        except ConnectionError:
            pytest.skip("Network not available")
        finally:
            set_iers_auto_download(False)

    @pytest.mark.unit
    def test_get_delta_t_cache_path(self):
        """Test cache path generation."""
        path = _get_delta_t_cache_path()
        assert path.endswith("deltat.data")
        assert self.temp_cache in path


class TestObservedDeltaTLookup:
    """Test observed Delta T value lookup and interpolation."""

    @pytest.fixture(autouse=True)
    def setup_data(self, tmp_path):
        """Set up test data."""
        clear_iers_cache()
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        # Create sample data file
        content = """ 1973  2  1  43.4724
 1973  3  1  43.5648
 2000  1  1  63.8285
 2000  2  1  63.8557
 2000  3  1  63.8804
 2020  1  1  69.3612
 2020  6  1  69.4386
 2020 12  1  69.3630
"""
        os.makedirs(os.path.join(self.temp_cache), exist_ok=True)
        file_path = os.path.join(self.temp_cache, "deltat.data")
        with open(file_path, "w") as f:
            f.write(content)

        # Load the data
        load_iers_data()

        yield

        clear_iers_cache()
        set_iers_cache_dir(None)

    @pytest.mark.unit
    def test_get_observed_delta_t_exact_date(self):
        """Test getting Delta T for an exact date in the data."""
        # Jan 1, 2000
        jd = julday(2000, 1, 1, 0.0)
        dt = get_observed_delta_t(jd)
        assert dt is not None
        assert dt == pytest.approx(63.8285, abs=0.001)

    @pytest.mark.unit
    def test_get_observed_delta_t_interpolated(self):
        """Test interpolation between data points."""
        # Mid-January 2000 (between Jan 1 and Feb 1)
        jd = julday(2000, 1, 15, 12.0)
        dt = get_observed_delta_t(jd)
        assert dt is not None
        # Should be between 63.8285 and 63.8557
        assert 63.8285 < dt < 63.8557

    @pytest.mark.unit
    def test_get_observed_delta_t_out_of_range_before(self):
        """Test that None is returned for dates before data range."""
        jd = julday(1970, 1, 1, 0.0)  # Before 1973
        dt = get_observed_delta_t(jd)
        assert dt is None

    @pytest.mark.unit
    def test_get_observed_delta_t_out_of_range_after(self):
        """Test that None is returned for dates after data range."""
        jd = julday(2025, 1, 1, 0.0)  # After our test data
        dt = get_observed_delta_t(jd)
        assert dt is None

    @pytest.mark.unit
    def test_is_observed_delta_t_available_in_range(self):
        """Test checking if observed data is available for a date in range."""
        jd = julday(2000, 6, 15, 12.0)
        assert is_observed_delta_t_available(jd) is True

    @pytest.mark.unit
    def test_is_observed_delta_t_available_out_of_range(self):
        """Test checking if observed data is available for a date out of range."""
        jd = julday(1960, 1, 1, 0.0)
        assert is_observed_delta_t_available(jd) is False

    @pytest.mark.unit
    def test_get_observed_delta_t_data_range(self):
        """Test getting the observed Delta T data range."""
        data_range = get_observed_delta_t_data_range()
        assert data_range is not None
        jd_start, jd_end = data_range
        # Should cover 1973 to 2020 based on our test data
        assert jd_start < julday(1973, 3, 1, 0.0)
        assert jd_end > julday(2020, 11, 1, 0.0)


class TestDeltaTIERSIntegration:
    """Test the integration of observed Delta T with get_delta_t_iers."""

    @pytest.fixture(autouse=True)
    def setup_data(self, tmp_path):
        """Set up test data."""
        clear_iers_cache()
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        # Create sample Delta T data file
        delta_t_content = """ 2000  1  1  63.8285
 2000  6  1  63.9691
 2001  1  1  64.0908
"""
        delta_t_path = os.path.join(self.temp_cache, "deltat.data")
        with open(delta_t_path, "w") as f:
            f.write(delta_t_content)

        # Load the data
        load_iers_data()

        yield

        clear_iers_cache()
        set_iers_cache_dir(None)

    @pytest.mark.unit
    def test_get_delta_t_iers_uses_observed(self):
        """Test that get_delta_t_iers uses observed data when available."""
        jd = julday(2000, 1, 1, 0.0)
        dt = get_delta_t_iers(jd)
        assert dt is not None
        # Should match our observed data
        assert dt == pytest.approx(63.8285, abs=0.01)

    @pytest.mark.unit
    def test_get_delta_t_iers_returns_none_out_of_range(self):
        """Test that get_delta_t_iers returns None for dates outside data."""
        jd = julday(1960, 1, 1, 0.0)
        dt = get_delta_t_iers(jd)
        # Should return None as we don't have UT1-UTC data either
        assert dt is None


class TestCacheInfo:
    """Test cache information includes Delta T data."""

    @pytest.fixture(autouse=True)
    def setup_data(self, tmp_path):
        """Set up test data."""
        clear_iers_cache()
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        # Create sample Delta T data file
        delta_t_content = """ 2000  1  1  63.8285
 2020  1  1  69.3612
"""
        delta_t_path = os.path.join(self.temp_cache, "deltat.data")
        with open(delta_t_path, "w") as f:
            f.write(delta_t_content)

        # Load the data
        load_iers_data()

        yield

        clear_iers_cache()
        set_iers_cache_dir(None)

    @pytest.mark.unit
    def test_cache_info_includes_delta_t(self):
        """Test that cache info includes Delta T information."""
        info = get_iers_cache_info()
        assert "delta_t_exists" in info
        assert "delta_t_age_days" in info
        assert "delta_t_entries" in info
        assert "delta_t_range_start_jd" in info
        assert "delta_t_range_end_jd" in info

    @pytest.mark.unit
    def test_cache_info_delta_t_values(self):
        """Test cache info Delta T values are correct."""
        info = get_iers_cache_info()
        assert info["delta_t_exists"] is True
        assert info["delta_t_entries"] == 2
        assert info["delta_t_range_start_jd"] is not None
        assert info["delta_t_range_end_jd"] is not None


class TestDeltaTNetworkIntegration:
    """Integration tests that require network access."""

    @pytest.fixture(autouse=True)
    def setup_and_cleanup(self, tmp_path):
        """Setup temporary cache directory and cleanup after tests."""
        clear_iers_cache()
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        yield

        clear_iers_cache()
        set_iers_cache_dir(None)

    @pytest.mark.network
    @pytest.mark.slow
    def test_full_delta_t_download_and_lookup(self):
        """Test downloading real IERS data and looking up Delta T."""
        set_iers_auto_download(True)
        try:
            # Download the actual data
            download_delta_t_data(force=True)
            load_iers_data()

            # Check that we have data
            info = get_iers_cache_info()
            assert info["delta_t_exists"] is True
            assert info["delta_t_entries"] > 100  # Should have many entries

            # Look up a known date (J2000)
            jd = julday(2000, 1, 1, 0.0)
            dt = get_observed_delta_t(jd)
            assert dt is not None
            # Delta T at J2000 should be around 63.8 seconds
            assert dt == pytest.approx(63.8285, abs=0.1)

            # Look up a recent date
            jd_recent = julday(2020, 1, 1, 0.0)
            dt_recent = get_observed_delta_t(jd_recent)
            assert dt_recent is not None
            # Delta T in 2020 should be around 69 seconds
            assert 68 < dt_recent < 70

        except ConnectionError:
            pytest.skip("Network not available")
        finally:
            set_iers_auto_download(False)

    @pytest.mark.network
    def test_observed_delta_t_accuracy(self):
        """Test that observed Delta T values are accurate."""
        set_iers_auto_download(True)
        try:
            download_delta_t_data(force=True)
            load_iers_data()

            # Known Delta T values from IERS (approximate)
            known_values = [
                # (year, month, day, expected_delta_t)
                (2000, 1, 1, 63.83),
                (2010, 1, 1, 66.07),
                (2015, 1, 1, 67.64),
            ]

            for year, month, day, expected in known_values:
                jd = julday(year, month, day, 0.0)
                dt = get_observed_delta_t(jd)
                if dt is not None:
                    assert dt == pytest.approx(expected, abs=0.2), (
                        f"Delta T for {year}-{month:02d}-{day:02d}: "
                        f"expected {expected}, got {dt}"
                    )

        except ConnectionError:
            pytest.skip("Network not available")
        finally:
            set_iers_auto_download(False)


class TestDeltatExWithIERS:
    """Test that deltat_ex() also uses IERS Delta T when enabled."""

    @pytest.fixture(autouse=True)
    def setup_data(self, tmp_path):
        """Set up test data."""
        clear_iers_cache()
        self.temp_cache = str(tmp_path / "iers_cache")
        os.makedirs(self.temp_cache, exist_ok=True)
        set_iers_cache_dir(self.temp_cache)

        # Create sample Delta T data file
        delta_t_content = """ 2000  1  1  63.8285
 2000  6  1  63.9691
 2001  1  1  64.0908
"""
        delta_t_path = os.path.join(self.temp_cache, "deltat.data")
        with open(delta_t_path, "w") as f:
            f.write(delta_t_content)

        # Load the data
        load_iers_data()

        yield

        clear_iers_cache()
        set_iers_cache_dir(None)
        # Reset IERS enabled state
        from libephemeris import set_iers_delta_t_enabled

        set_iers_delta_t_enabled(False)

    @pytest.mark.unit
    def test_swe_deltat_ex_uses_iers_when_enabled(self):
        """Test that deltat_ex uses IERS data when enabled."""
        from libephemeris import deltat_ex, set_iers_delta_t_enabled

        jd = julday(2000, 1, 1, 0.0)

        # Enable IERS Delta T
        set_iers_delta_t_enabled(True)

        dt = deltat_ex(jd)
        dt_seconds = dt * 86400

        # Should match our observed data (63.8285 seconds)
        assert dt_seconds == pytest.approx(63.8285, abs=0.01)

    @pytest.mark.unit
    def test_swe_deltat_and_swe_deltat_ex_consistent(self):
        """Test that deltat and deltat_ex return the same value."""
        from libephemeris import deltat, deltat_ex, set_iers_delta_t_enabled

        jd = julday(2000, 1, 1, 0.0)

        # Enable IERS Delta T
        set_iers_delta_t_enabled(True)

        dt = deltat(jd)
        dt_ex = deltat_ex(jd)

        # Both should return identical values
        assert dt == pytest.approx(dt_ex, abs=1e-10)


class TestIERSPublicAPIExports:
    """Test that IERS functions are properly exported from the main module."""

    @pytest.mark.unit
    def test_iers_configuration_functions_exported(self):
        """Test that IERS configuration functions are exported."""
        from libephemeris import (
            set_iers_delta_t_enabled,
            get_iers_delta_t_enabled,
            set_iers_cache_dir,
            get_iers_cache_dir,
            set_iers_auto_download,
            get_iers_auto_download,
        )

        # Just verify they're callable
        assert callable(set_iers_delta_t_enabled)
        assert callable(get_iers_delta_t_enabled)
        assert callable(set_iers_cache_dir)
        assert callable(get_iers_cache_dir)
        assert callable(set_iers_auto_download)
        assert callable(get_iers_auto_download)

    @pytest.mark.unit
    def test_iers_download_functions_exported(self):
        """Test that IERS download functions are exported."""
        from libephemeris import (
            download_iers_finals,
            download_leap_seconds,
            download_delta_t_data,
            load_iers_data,
        )

        # Just verify they're callable
        assert callable(download_iers_finals)
        assert callable(download_leap_seconds)
        assert callable(download_delta_t_data)
        assert callable(load_iers_data)

    @pytest.mark.unit
    def test_iers_lookup_functions_exported(self):
        """Test that IERS lookup functions are exported."""
        from libephemeris import (
            get_observed_delta_t,
            get_observed_delta_t_data_range,
            is_observed_delta_t_available,
            get_delta_t_iers,
            get_ut1_utc,
            get_tai_utc,
            get_iers_data_range,
            is_iers_data_available,
        )

        # Just verify they're callable
        assert callable(get_observed_delta_t)
        assert callable(get_observed_delta_t_data_range)
        assert callable(is_observed_delta_t_available)
        assert callable(get_delta_t_iers)
        assert callable(get_ut1_utc)
        assert callable(get_tai_utc)
        assert callable(get_iers_data_range)
        assert callable(is_iers_data_available)

    @pytest.mark.unit
    def test_iers_cache_functions_exported(self):
        """Test that IERS cache management functions are exported."""
        import libephemeris

        # Just verify they're accessible as attributes
        assert callable(libephemeris.clear_iers_cache)
        assert callable(libephemeris.delete_iers_cache_files)
        assert callable(libephemeris.get_iers_cache_info)
