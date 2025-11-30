"""
Thread safety test for EphemerisContext.

This test verifies that multiple threads can use separate EphemerisContext
instances concurrently without interference.
"""

import concurrent.futures
from libephemeris import EphemerisContext, SE_SUN, SE_MOON, SEFLG_TOPOCTR


def test_concurrent_contexts():
    """Test that multiple contexts can be used concurrently."""
    
    def calculate_for_location(worker_id, lon, lat):
        """Worker function that creates its own context."""
        ctx = EphemerisContext()
        ctx.set_topo(lon, lat, 0)
        
        # Calculate Sun position (topocentric)
        jd = 2451545.0  # J2000.0
        sun_pos, _ = ctx.calc_ut(jd, SE_SUN, SEFLG_TOPOCTR)
        moon_pos, _ = ctx.calc_ut(jd, SE_MOON, SEFLG_TOPOCTR)
        
        return {
            'worker_id': worker_id,
            'lon': lon,
            'lat': lat,
            'sun_lon': sun_pos[0],
            'moon_lon': moon_pos[0]
        }
    
    # Test locations with different coordinates
    locations = [
        (0, 0.0, 0.0),      # Equator
        (1, 12.5, 41.9),    # Rome
        (2, -0.1, 51.5),    # London
        (3, 139.7, 35.7),   # Tokyo
        (4, -74.0, 40.7),   # New York
        (5, 151.2, -33.9),  # Sydney
        (6, -122.4, 37.8),  # San Francisco
        (7, 2.3, 48.9),     # Paris
        (8, 13.4, 52.5),    # Berlin
        (9, -43.2, -22.9),  # Rio de Janeiro
    ]
    
    print("Testing concurrent calculations with 10 threads...")
    
    # Execute calculations concurrently
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(calculate_for_location, wid, lon, lat) 
                   for wid, lon, lat in locations]
        results = [f.result() for f in futures]
    
    # Verify each thread got different results based on their location
    print("\nResults:")
    for r in sorted(results, key=lambda x: x['worker_id']):
        print(f"Worker {r['worker_id']}: "
              f"Lon={r['lon']:7.2f}°, Lat={r['lat']:6.2f}° → "
              f"Sun={r['sun_lon']:6.2f}°, Moon={r['moon_lon']:6.2f}°")
    
    # Verify all Sun positions are NOT identical (topocentric varies by location)
    sun_positions = [r['sun_lon'] for r in results]
    unique_sun = len(set(round(pos, 4) for pos in sun_positions))
    
    print(f"\n✓ Thread safety test passed!")
    print(f"  {len(results)} threads executed concurrently")
    print(f"  {unique_sun} unique Sun positions (topocentric variation)")
    
    assert unique_sun > 1, "Expected different topocentric positions for different locations"
    
    return results


def test_resource_sharing():
    """Test that contexts share expensive resources."""
    
    print("\nTesting resource sharing...")
    
    ctx1 = EphemerisContext()
    ctx2 = EphemerisContext()
    
    # Both should use same SpiceKernel instance
    planets1 = ctx1.get_planets()
    planets2 = ctx2.get_planets()
    
    # Check if they're the same object (shared)
    is_shared = planets1 is planets2
    
    print(f"  Planets ephemeris shared: {is_shared}")
    print(f"  Timescale shared: {ctx1.get_timescale() is ctx2.get_timescale()}")
    
    assert is_shared, "Ephemeris should be shared between contexts"
    
    print("✓ Resource sharing test passed!")


def test_isolation():
    """Test that contexts are isolated from each other."""
    
    print("\nTesting context isolation...")
    
    ctx1 = EphemerisContext()
    ctx2 = EphemerisContext()
    
    # Set different topocentric locations
    ctx1.set_topo(12.5, 41.9, 0)  # Rome
    ctx2.set_topo(-0.1, 51.5, 0)  # London
    
    # Verify they're independent
    topo1 = ctx1.get_topo()
    topo2 = ctx2.get_topo()
    
    assert topo1.latitude.degrees != topo2.latitude.degrees
    assert topo1.longitude.degrees != topo2.longitude.degrees
    
    # Set different sidereal modes
    ctx1.set_sid_mode(1)  # Lahiri
    ctx2.set_sid_mode(0)  # Fagan-Bradley
    
    assert ctx1.get_sid_mode() != ctx2.get_sid_mode()
    
    print("  ✓ Topo locations isolated")
    print("  ✓ Sidereal modes isolated")
    print("✓ Context isolation test passed!")


if __name__ == "__main__":
    print("=" * 60)
    print("EphemerisContext Thread Safety Tests")
    print("=" * 60)
    
    test_resource_sharing()
    test_isolation()
    test_concurrent_contexts()
    
    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✓")
    print("=" * 60)
