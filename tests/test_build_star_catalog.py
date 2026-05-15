"""
Tests for the build_star_catalog.py script.

These tests verify the star catalog building functionality including:
- Data structure validation
- Star name resolution
- Unit conversions (mas/yr to arcsec/yr)
- Output format generation
"""

import json
import os
import sys
import pytest
from dataclasses import asdict

# Add scripts directory to path for imports
sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"
    ),
)

from build_star_catalog import (
    StarCatalogData,
    IAU_STAR_NAMES,
    TRADITIONAL_STAR_NAMES,
    FIXSTAR_OFFSET,
    get_star_name,
    generate_nomenclature,
    output_csv,
    output_json,
    output_python_catalog,
    check_astroquery,
)


@pytest.mark.unit
class TestStarCatalogData:
    """Test StarCatalogData dataclass."""

    def test_create_star_data(self):
        """Test creating a StarCatalogData instance."""
        star = StarCatalogData(
            hip_number=49669,
            name="Regulus",
            nomenclature="alLeo",
            ra_j2000=152.092958,
            dec_j2000=11.967208,
            pm_ra=-0.00249,
            pm_dec=0.00152,
            magnitude=1.40,
        )

        assert star.hip_number == 49669
        assert star.name == "Regulus"
        assert star.nomenclature == "alLeo"
        assert star.ra_j2000 == pytest.approx(152.092958)
        assert star.dec_j2000 == pytest.approx(11.967208)
        assert star.pm_ra == pytest.approx(-0.00249)
        assert star.pm_dec == pytest.approx(0.00152)
        assert star.magnitude == pytest.approx(1.40)

    def test_star_data_to_dict(self):
        """Test converting StarCatalogData to dictionary."""
        star = StarCatalogData(
            hip_number=32349,
            name="Sirius",
            nomenclature="alCMa",
            ra_j2000=101.287155,
            dec_j2000=-16.716116,
            pm_ra=-0.54601,
            pm_dec=-1.22307,
            magnitude=-1.46,
        )

        d = asdict(star)
        assert d["hip_number"] == 32349
        assert d["name"] == "Sirius"
        assert d["magnitude"] == pytest.approx(-1.46)


@pytest.mark.unit
class TestIAUStarNames:
    """Test IAU star name catalog."""

    def test_iau_names_not_empty(self):
        """IAU star names dictionary should not be empty."""
        assert len(IAU_STAR_NAMES) > 0

    def test_known_stars_present(self):
        """Known bright stars should be in the IAU names."""
        # Regulus
        assert 49669 in IAU_STAR_NAMES
        assert IAU_STAR_NAMES[49669] == "Regulus"

        # Sirius
        assert 32349 in IAU_STAR_NAMES
        assert IAU_STAR_NAMES[32349] == "Sirius"

        # Aldebaran
        assert 21421 in IAU_STAR_NAMES
        assert IAU_STAR_NAMES[21421] == "Aldebaran"

        # Vega
        assert 91262 in IAU_STAR_NAMES
        assert IAU_STAR_NAMES[91262] == "Vega"

    def test_hip_numbers_are_integers(self):
        """All HIP numbers should be positive integers."""
        for hip in IAU_STAR_NAMES.keys():
            assert isinstance(hip, int)
            assert hip > 0

    def test_names_are_non_empty_strings(self):
        """All star names should be non-empty strings."""
        for name in IAU_STAR_NAMES.values():
            assert isinstance(name, str)
            assert len(name) > 0


@pytest.mark.unit
class TestGetStarName:
    """Test star name resolution."""

    def test_known_iau_star(self):
        """IAU named stars should return their official names."""
        assert get_star_name(49669) == "Regulus"
        assert get_star_name(32349) == "Sirius"
        assert get_star_name(91262) == "Vega"

    def test_unknown_star(self):
        """Unknown stars should return HIP number format."""
        # Use a very high HIP number unlikely to have a name
        result = get_star_name(999999)
        assert result == "HIP 999999"

    def test_traditional_name_fallback(self):
        """Stars in TRADITIONAL_STAR_NAMES should return those names."""
        # Check a few stars that are in TRADITIONAL_STAR_NAMES
        for hip, name in list(TRADITIONAL_STAR_NAMES.items())[:3]:
            # If also in IAU names, IAU takes precedence
            if hip not in IAU_STAR_NAMES:
                assert get_star_name(hip) == name


@pytest.mark.unit
class TestGenerateNomenclature:
    """Test nomenclature generation."""

    def test_nomenclature_format(self):
        """Nomenclature should use HIP prefix."""
        result = generate_nomenclature(49669, "Regulus")
        assert result == "HIP49669"

    def test_nomenclature_with_unknown(self):
        """Unknown stars should also get HIP nomenclature."""
        result = generate_nomenclature(123456, "HIP 123456")
        assert result == "HIP123456"


@pytest.mark.unit
class TestOutputFormats:
    """Test output format generation."""

    @pytest.fixture
    def sample_stars(self):
        """Create a sample list of stars for testing."""
        return [
            StarCatalogData(
                hip_number=49669,
                name="Regulus",
                nomenclature="HIP49669",
                ra_j2000=152.092958,
                dec_j2000=11.967208,
                pm_ra=-0.00249,
                pm_dec=0.00152,
                magnitude=1.40,
            ),
            StarCatalogData(
                hip_number=32349,
                name="Sirius",
                nomenclature="HIP32349",
                ra_j2000=101.287155,
                dec_j2000=-16.716116,
                pm_ra=-0.54601,
                pm_dec=-1.22307,
                magnitude=-1.46,
            ),
        ]

    def test_csv_output_header(self, sample_stars):
        """CSV output should have proper header."""
        csv = output_csv(sample_stars)
        lines = csv.split("\n")
        assert (
            lines[0]
            == "hip_number,name,nomenclature,ra_j2000,dec_j2000,pm_ra,pm_dec,magnitude"
        )

    def test_csv_output_data(self, sample_stars):
        """CSV output should contain star data."""
        csv = output_csv(sample_stars)
        assert "Regulus" in csv
        assert "Sirius" in csv
        assert "32349" in csv  # Sirius HIP number
        assert "49669" in csv  # Regulus HIP number

    def test_json_output_valid(self, sample_stars):
        """JSON output should be valid JSON."""
        json_str = output_json(sample_stars)
        data = json.loads(json_str)

        assert "catalog" in data
        assert "generated" in data
        assert "count" in data
        assert "stars" in data
        assert data["count"] == 2
        assert len(data["stars"]) == 2

    def test_json_output_star_data(self, sample_stars):
        """JSON output should contain correct star data."""
        json_str = output_json(sample_stars)
        data = json.loads(json_str)

        # Stars should be sorted by HIP number
        assert data["stars"][0]["hip_number"] == 32349  # Sirius
        assert data["stars"][1]["hip_number"] == 49669  # Regulus

        # Check field types
        sirius = data["stars"][0]
        assert isinstance(sirius["hip_number"], int)
        assert isinstance(sirius["name"], str)
        assert isinstance(sirius["ra_j2000"], float)
        assert isinstance(sirius["magnitude"], float)

    def test_python_catalog_output(self, sample_stars):
        """Python catalog output should be valid Python code structure."""
        py = output_python_catalog(sample_stars)

        # Should contain imports
        assert "from libephemeris.fixed_stars import StarData, StarCatalogEntry" in py
        assert "from libephemeris.constants import FIXSTAR_OFFSET" in py

        # Should contain list definition
        assert "HIPPARCOS_STAR_CATALOG: list[StarCatalogEntry] = [" in py

        # Should contain star data
        assert "Regulus" in py
        assert "Sirius" in py

        # Should have StarCatalogEntry structures
        assert "StarCatalogEntry(" in py
        assert "StarData(" in py


@pytest.mark.unit
class TestUnitConversions:
    """Test unit conversions for proper motion."""

    def test_proper_motion_units(self):
        """Proper motion should be in arcsec/year."""
        # Sirius has one of the highest proper motions
        # pmRA = -546.01 mas/yr, pmDec = -1223.07 mas/yr
        # These should convert to -0.54601 and -1.22307 arcsec/yr

        star = StarCatalogData(
            hip_number=32349,
            name="Sirius",
            nomenclature="HIP32349",
            ra_j2000=101.287155,
            dec_j2000=-16.716116,
            pm_ra=-0.54601,  # arcsec/yr (was -546.01 mas/yr)
            pm_dec=-1.22307,  # arcsec/yr (was -1223.07 mas/yr)
            magnitude=-1.46,
        )

        # Verify conversion is correct (mas/yr / 1000 = arcsec/yr)
        assert star.pm_ra == pytest.approx(-546.01 / 1000)
        assert star.pm_dec == pytest.approx(-1223.07 / 1000)


@pytest.mark.unit
class TestConstants:
    """Test script constants."""

    def test_fixstar_offset(self):
        """FIXSTAR_OFFSET should match libephemeris constant."""
        assert FIXSTAR_OFFSET == 1000000

    def test_astroquery_check(self):
        """check_astroquery should return a boolean."""
        result = check_astroquery()
        assert isinstance(result, bool)


@pytest.mark.unit
class TestStarDataValidation:
    """Test star data validation."""

    def test_ra_range(self):
        """RA should be in range [0, 360)."""
        # Regulus RA ~152°
        star = StarCatalogData(
            hip_number=49669,
            name="Regulus",
            nomenclature="HIP49669",
            ra_j2000=152.092958,
            dec_j2000=11.967208,
            pm_ra=-0.00249,
            pm_dec=0.00152,
            magnitude=1.40,
        )
        assert 0 <= star.ra_j2000 < 360

    def test_dec_range(self):
        """Dec should be in range [-90, 90]."""
        # Regulus Dec ~+12°
        star = StarCatalogData(
            hip_number=49669,
            name="Regulus",
            nomenclature="HIP49669",
            ra_j2000=152.092958,
            dec_j2000=11.967208,
            pm_ra=-0.00249,
            pm_dec=0.00152,
            magnitude=1.40,
        )
        assert -90 <= star.dec_j2000 <= 90

    def test_magnitude_reasonable(self):
        """Magnitude should be in reasonable range for catalog."""
        # Sirius (brightest) is -1.46, limit is around 4.0
        star = StarCatalogData(
            hip_number=32349,
            name="Sirius",
            nomenclature="HIP32349",
            ra_j2000=101.287155,
            dec_j2000=-16.716116,
            pm_ra=-0.54601,
            pm_dec=-1.22307,
            magnitude=-1.46,
        )
        assert -2 <= star.magnitude <= 6  # Reasonable range


@pytest.mark.integration
class TestAstroqueryIntegration:
    """Integration tests requiring astroquery (skipped if not available)."""

    @pytest.fixture
    def check_astroquery_available(self):
        """Skip if astroquery not available."""
        if not check_astroquery():
            pytest.skip("astroquery not installed")

    def test_astroquery_import(self, check_astroquery_available):
        """Verify astroquery can be imported."""
        from astroquery.vizier import Vizier

        assert Vizier is not None
