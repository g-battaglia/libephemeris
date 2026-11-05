//! Allocation regression for successful warmed native positions and clocks/houses.
//!
//! This binary contains one test so unrelated tests cannot pollute the process
//! allocator counters. Acquisition, expected outputs and warmups are outside
//! the counted regions. Python object allocation is not measured.

#[allow(
    dead_code,
    reason = "Shared test data helpers also serve native unit tests."
)]
#[path = "../src/foundation/test_data.rs"]
mod test_data;
use libephemeris_core::{
    config::Runtime,
    constants::{FLG_EQUATORIAL, FLG_SPEED, FLG_SWIEPH, FLG_XYZ, MARS},
    context::{EphemerisContext, NetworkPolicy, ProviderMode},
    leb::{IndexedLimits, LebFile, LebReadPolicy},
    native::{Position, Session},
};
use std::{
    alloc::{GlobalAlloc, Layout, System},
    hint::{black_box, spin_loop},
    path::PathBuf,
    ptr::NonNull,
    sync::atomic::{AtomicUsize, Ordering},
};

/// Operation counts are monotonic and global across every thread. Requiring
/// zero operations also requires zero allocated, released or resized bytes.
#[derive(Debug, Default, PartialEq, Eq)]
struct AllocationCounts {
    allocations: usize,
    zeroed_allocations: usize,
    deallocations: usize,
    reallocations: usize,
}

impl AllocationCounts {
    fn change_from(self, before: Self) -> Self {
        let delta = |current: usize, previous: usize| {
            // Counter saturation must fail the test, never certify zero work.
            // This check runs outside allocator callbacks, where panic is safe.
            assert!(current < usize::MAX && previous < usize::MAX);
            current
                .checked_sub(previous)
                .expect("monotonic allocation counter")
        };
        Self {
            allocations: delta(self.allocations, before.allocations),
            zeroed_allocations: delta(self.zeroed_allocations, before.zeroed_allocations),
            deallocations: delta(self.deallocations, before.deallocations),
            reallocations: delta(self.reallocations, before.reallocations),
        }
    }
}

struct CountingAllocator {
    allocations: AtomicUsize,
    zeroed_allocations: AtomicUsize,
    deallocations: AtomicUsize,
    reallocations: AtomicUsize,
}

impl CountingAllocator {
    const fn new() -> Self {
        Self {
            allocations: AtomicUsize::new(0),
            zeroed_allocations: AtomicUsize::new(0),
            deallocations: AtomicUsize::new(0),
            reallocations: AtomicUsize::new(0),
        }
    }

    fn record(counter: &AtomicUsize) {
        // No allocation, callback, reset, thread-local filter or panic. Count
        // attempts too, regardless of the underlying allocator's result.
        let _ = counter.fetch_update(Ordering::SeqCst, Ordering::SeqCst, |value| {
            Some(value.saturating_add(1))
        });
    }

    fn counts(&self) -> AllocationCounts {
        AllocationCounts {
            allocations: self.allocations.load(Ordering::SeqCst),
            zeroed_allocations: self.zeroed_allocations.load(Ordering::SeqCst),
            deallocations: self.deallocations.load(Ordering::SeqCst),
            reallocations: self.reallocations.load(Ordering::SeqCst),
        }
    }
}

// SAFETY: System retains ownership of every block. Each method forwards the
// caller's pointer, layout and size unchanged and returns System's result.
// Atomic bookkeeping never accesses memory in those blocks, extends a lifetime,
// introduces aliases, allocates or unwinds. The caller retains GlobalAlloc's
// length, alignment, ownership and lifetime obligations.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        Self::record(&self.allocations);
        // SAFETY: the caller's allocation layout is forwarded unchanged.
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        Self::record(&self.zeroed_allocations);
        // SAFETY: System receives the caller's exact valid allocation layout.
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        Self::record(&self.deallocations);
        // SAFETY: the caller owns this live System block with this exact layout;
        // no references are retained or used by the bookkeeping.
        unsafe { System.dealloc(pointer, layout) }
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        Self::record(&self.reallocations);
        // SAFETY: pointer, original layout and nonzero new size retain the
        // caller's ownership/aliasing contract and are forwarded unchanged.
        unsafe { System.realloc(pointer, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator::new();

/// Own a real block for positive controls without optimizer-elidable containers.
struct ProbeAllocation {
    pointer: NonNull<u8>,
    layout: Layout,
}

impl ProbeAllocation {
    fn new(zeroed: bool) -> Self {
        let layout = Layout::from_size_align(64, 8).unwrap();
        // SAFETY: the nonzero layout is valid. This object takes sole ownership
        // of the returned block; resize/Drop use the same allocator and layout.
        let pointer = unsafe {
            if zeroed {
                ALLOCATOR.alloc_zeroed(layout)
            } else {
                ALLOCATOR.alloc(layout)
            }
        };
        Self {
            pointer: NonNull::new(pointer).expect("positive-control allocation"),
            layout,
        }
    }

    fn resize(&mut self, size: usize) {
        let layout = Layout::from_size_align(size, self.layout.align()).unwrap();
        assert!(size > 0);
        // SAFETY: self uniquely owns this live block, no references to its bytes
        // exist, and the original layout and nonzero new size are exact. If
        // reallocation fails, self still owns the original pointer for Drop.
        let pointer = unsafe { ALLOCATOR.realloc(self.pointer.as_ptr(), self.layout, size) };
        self.pointer = NonNull::new(pointer).expect("positive-control reallocation");
        self.layout = layout;
    }
}

impl Drop for ProbeAllocation {
    fn drop(&mut self) {
        // SAFETY: this is the sole owner of a still-live System allocation;
        // resize updates both pointer and layout together after success.
        unsafe { ALLOCATOR.dealloc(self.pointer.as_ptr(), self.layout) }
    }
}

fn measure_operation<T>(operation: impl FnOnce() -> T) -> (T, AllocationCounts) {
    let before = ALLOCATOR.counts();
    let value = operation();
    (value, ALLOCATOR.counts().change_from(before))
}

struct WorkerCompletion<'a>(&'a AtomicUsize);

impl Drop for WorkerCompletion<'_> {
    fn drop(&mut self) {
        // A failed positive-control allocation must reach join and fail the
        // test, rather than leave the measuring thread waiting indefinitely.
        self.0.store(3, Ordering::SeqCst);
    }
}

fn assert_allocator_controls() {
    // An ordinary Rust container also proves this wrapper is installed as the
    // global allocator; direct calibration calls alone would not establish it.
    let (container, counts) = measure_operation(|| black_box(Vec::<u8>::with_capacity(64)));
    assert_eq!(
        counts,
        AllocationCounts {
            allocations: 1,
            ..AllocationCounts::default()
        }
    );
    let (_, counts) = measure_operation(|| drop(container));
    assert_eq!(
        counts,
        AllocationCounts {
            deallocations: 1,
            ..AllocationCounts::default()
        }
    );

    for zeroed in [false, true] {
        let (mut block, counts) = measure_operation(|| ProbeAllocation::new(zeroed));
        assert_eq!(
            counts,
            AllocationCounts {
                allocations: usize::from(!zeroed),
                zeroed_allocations: usize::from(zeroed),
                ..AllocationCounts::default()
            }
        );
        // Growing, unchanged-size and shrinking reallocations each count as an
        // event; equal sizes or net byte balance must not mask activity.
        for size in [128, 128, 32] {
            let (_, counts) = measure_operation(|| block.resize(size));
            assert_eq!(
                counts,
                AllocationCounts {
                    reallocations: 1,
                    ..AllocationCounts::default()
                }
            );
        }
        // The block predates this measured region: its release must still count.
        let (_, counts) = measure_operation(|| drop(block));
        assert_eq!(
            counts,
            AllocationCounts {
                deallocations: 1,
                ..AllocationCounts::default()
            }
        );
    }

    let phase = AtomicUsize::new(0);
    let wait_for = |wanted| {
        while phase.load(Ordering::SeqCst) != wanted {
            spin_loop();
        }
    };
    std::thread::scope(|scope| {
        let worker = scope.spawn(|| {
            phase.store(1, Ordering::SeqCst);
            wait_for(2);
            {
                let _completion = WorkerCompletion(&phase);
                drop(black_box(Vec::<u8>::with_capacity(64)));
            }
            wait_for(4);
        });
        wait_for(1);
        let before = ALLOCATOR.counts();
        phase.store(2, Ordering::SeqCst);
        wait_for(3);
        let after = ALLOCATOR.counts();
        // Release/join outside the measured interval, before any assertion can
        // unwind. Only the worker performs heap operations within this interval.
        phase.store(4, Ordering::SeqCst);
        worker.join().unwrap();
        assert_eq!(
            after.change_from(before),
            AllocationCounts {
                allocations: 1,
                deallocations: 1,
                ..AllocationCounts::default()
            }
        );
    });
}

fn assert_exact(actual: &Position, expected: &Position) {
    assert_eq!(
        actual.values.map(f64::to_bits),
        expected.values.map(f64::to_bits)
    );
    assert_eq!(actual.returned_flags, expected.returned_flags);
    assert_eq!(actual.diagnostic, expected.diagnostic);
    assert_eq!(actual.warnings, expected.warnings);
    assert_eq!(actual.messages, expected.messages);
}

fn counted_calls(mut calculate: impl FnMut() -> Position, expected: &Position) -> AllocationCounts {
    for _ in 0..32 {
        assert_exact(&calculate(), expected);
    }
    let before = ALLOCATOR.counts();
    for _ in 0..512 {
        let actual = calculate();
        assert_exact(&actual, expected);
        // Drop all owned output within the counted region as well.
        black_box(actual);
    }
    ALLOCATOR.counts().change_from(before)
}

/// Measure the complete returned value, including its destruction. The scalar
/// clock and fixed-layout SDK house outputs have no caller-owned heap payload.
fn counted_warmed_values<T>(
    mut calculate: impl FnMut() -> T,
    expected: &T,
    assert_same: impl Fn(&T, &T),
) -> AllocationCounts {
    for _ in 0..32 {
        assert_same(&calculate(), expected);
    }
    let before = ALLOCATOR.counts();
    for _ in 0..512 {
        let actual = calculate();
        assert_same(&actual, expected);
        black_box(actual);
    }
    ALLOCATOR.counts().change_from(before)
}

fn assert_sidereal_and_house_allocations() {
    // A sealed, configured native session needs no image for these supplied
    // sidereal angles or tropical P/O/G houses. Constant Delta-T makes both
    // represented clocks identical and avoids an unrelated data acquisition.
    for mode in [ProviderMode::NativeJpl, ProviderMode::SealedLeb] {
        let session = Session::from_context(
            EphemerisContext::new()
                .with_mode(mode)
                .with_network_policy(Some(NetworkPolicy::Sealed))
                .with_iers_delta_t(Some(false))
                .with_user_delta_t(Some(0.0)),
        );
        let dates = [
            ("J2000", 2_451_545.0),
            ("2300", 2_451_545.0 + 3.0 * 36_525.0),
            ("locator node switch", 2_451_545.0 + 5.0 * 36_525.0),
            ("extended start", -3_100_015.5),
            ("extended end", 8_000_016.5),
            ("model start", 2_451_545.0 - 2000.0 * 36_525.0),
            ("model end", 2_451_545.0 + 2000.0 * 36_525.0),
        ];
        for (name, jd) in dates {
            let expected = session.sidtime0(jd, 23.4, 0.0).unwrap();
            assert!(expected.is_finite());
            let counts = counted_warmed_values(
                || {
                    session
                        .sidtime0(black_box(jd), black_box(23.4), black_box(0.0))
                        .unwrap()
                },
                &expected,
                |actual, expected| assert_eq!(actual.to_bits(), expected.to_bits()),
            );
            assert_eq!(
                counts,
                AllocationCounts::default(),
                "sidtime0, {mode:?}, {name}"
            );

            // Equatorial latitude keeps all three constructions outside polar
            // refusals. Call the actual SDK paths, including their internal
            // container conversions; do not replace them with the ARMC kernels.
            for system in *b"PO" {
                let expected = session.houses(jd, 0.0, 12.0, system).unwrap();
                assert!(expected.cusps.iter().all(|value| value.is_finite()));
                assert!(expected.angles.iter().all(|value| value.is_finite()));
                let counts = counted_warmed_values(
                    || {
                        session
                            .houses(
                                black_box(jd),
                                black_box(0.0),
                                black_box(12.0),
                                black_box(system),
                            )
                            .unwrap()
                    },
                    &expected,
                    |actual, expected| {
                        assert_eq!(
                            actual.cusps.map(f64::to_bits),
                            expected.cusps.map(f64::to_bits)
                        );
                        assert_eq!(
                            actual.angles.map(f64::to_bits),
                            expected.angles.map(f64::to_bits)
                        );
                    },
                );
                assert_eq!(
                    counts,
                    AllocationCounts::default(),
                    "houses {}, {mode:?}, {name}",
                    char::from(system)
                );
            }

            // G has a distinct fixed 36-sector SDK result. It also reports rates;
            // its complete acquisition and all returned fields remain measured.
            let expected = session.houses_gauquelin(jd, 0, 0.0, 12.0).unwrap();
            assert!(expected.sectors.iter().all(|value| value.is_finite()));
            assert!(expected.angles.iter().all(|value| value.is_finite()));
            let counts = counted_warmed_values(
                || {
                    session
                        .houses_gauquelin(black_box(jd), 0, black_box(0.0), black_box(12.0))
                        .unwrap()
                },
                &expected,
                |actual, expected| {
                    assert_eq!(
                        actual.sectors.map(f64::to_bits),
                        expected.sectors.map(f64::to_bits)
                    );
                    assert_eq!(
                        actual.angles.map(f64::to_bits),
                        expected.angles.map(f64::to_bits)
                    );
                    assert_eq!(
                        actual.sector_speeds.map(f64::to_bits),
                        expected.sector_speeds.map(f64::to_bits)
                    );
                    assert_eq!(
                        actual.angle_speeds.map(f64::to_bits),
                        expected.angle_speeds.map(f64::to_bits)
                    );
                },
            );
            assert_eq!(
                counts,
                AllocationCounts::default(),
                "Gauquelin houses, {mode:?}, {name}"
            );
        }
    }
}

#[test]
fn warmed_session_and_runtime_positions_have_no_temporary_allocations() {
    assert_allocator_controls();
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(PathBuf::from)
        .unwrap_or_else(test_data::core_leb_path);
    for indexed in [false, true] {
        let limits = IndexedLimits::default();
        let file = if indexed {
            LebFile::open_indexed(&path, limits)
        } else {
            LebFile::open(&path)
        }
        .expect("provision the identified validation asset before this regression");
        assert!(test_data::admitted_core_hash(file.sha256().unwrap()));
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_iers_delta_t(Some(false));
        context.set_tier_name("base").unwrap();
        context
            .set_leb_read_policy(Some(if indexed {
                LebReadPolicy::Indexed(limits)
            } else {
                LebReadPolicy::Owned
            }))
            .unwrap();
        context.attach_leb_file(file);
        let session = Session::from_context(context.clone());
        let runtime = Runtime::from_context(context);
        for is_ut in [false, true] {
            for output in [0, FLG_EQUATORIAL, FLG_XYZ, FLG_EQUATORIAL | FLG_XYZ] {
                for speed in [0, FLG_SPEED] {
                    let flags = FLG_SWIEPH | output | speed;
                    let calculate = |use_runtime| match (use_runtime, is_ut) {
                        (false, false) => session.calc(2451545.0, MARS, flags),
                        (false, true) => session.calc_ut(2451545.0, MARS, flags),
                        (true, false) => runtime.calc(2451545.0, MARS, flags),
                        (true, true) => runtime.calc_ut(2451545.0, MARS, flags),
                    };
                    let expected = calculate(false).unwrap();
                    assert!(expected.diagnostic.is_empty());
                    assert!(expected.warnings.is_empty());
                    assert!(expected.messages.is_empty());
                    for use_runtime in [false, true] {
                        let stats = counted_calls(|| calculate(use_runtime).unwrap(), &expected);
                        assert_eq!(
                            stats,
                            AllocationCounts::default(),
                            "indexed={indexed}, runtime={use_runtime}, UT={is_ut}, flags={flags}"
                        );
                    }
                }
            }
        }
    }
    assert_sidereal_and_house_allocations();
}
