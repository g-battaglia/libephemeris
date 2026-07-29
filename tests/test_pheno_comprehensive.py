class TestUranusPhotometricLatitudeBackends:
    """The Uranus sub-Earth photometric latitude (Mallama & Hilton 2018
    Eq. 15) is applied on BOTH backends: the Skyfield path used to leave the
    ecliptic coordinates at zero for Uranus, degenerating the term to a
    constant and splitting the backends by tens of millimagnitudes."""

    def test_backends_agree_on_uranus_magnitude(self):
        import libephemeris as le

        prev = le.get_calc_mode()
        try:
            for y in (1950, 2030, 2070):
                jd = le.julday(y, 9, 15, 0.0)
                le.set_calc_mode("leb")
                m_leb = le.pheno_ut(jd, le.URANUS, le.FLG_SWIEPH)[4]
                le.set_calc_mode("skyfield")
                m_sky = le.pheno_ut(jd, le.URANUS, le.FLG_SWIEPH)[4]
                assert abs(m_leb - m_sky) < 0.001
        finally:
            le.set_calc_mode(prev)
