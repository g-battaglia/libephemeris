//! Reuse owned Session/kernel snapshots in a bounded synchronous worker pool.
//! Run `cargo run --example bounded_native_workers -- /path/to/de440.bsp`.
//! The Sun calculation uses the primary kernel's existing Sun/Earth segments.
//! An async server can await its own reply adapter. Dropping a reply receiver
//! abandons the result, while queued/running work retains its Session and its
//! worker/queue place until completion. Shutdown drains all submitted work.
use libephemeris_core::{
    constants::{FLG_SPEED, FLG_SWIEPH, SUN},
    context::{EphemerisContext, NetworkPolicy, ProviderMode},
    errors::Error,
    native::{Position, Session},
    spk::{SpkRegistry, register_all_naif_targets},
    time::DeltaTModel,
};
use std::{
    error::Error as StdError,
    sync::{Arc, Mutex, mpsc},
    thread::{self, JoinHandle},
};

/// Owned settings/mappings and reply; no kernel is copied by submission.
struct Calculation {
    session: Arc<Session>,
    jd_ut: f64,
    body: i32,
    flags: i32,
    reply: mpsc::SyncSender<Result<Position, Error>>,
}

/// Fixed workers and equally bounded queued jobs. Computation runs without the
/// receiver mutex, including when the result's consumer has disconnected.
struct NativeWorkers {
    sender: Option<mpsc::SyncSender<Calculation>>,
    threads: Vec<JoinHandle<()>>,
}
impl NativeWorkers {
    fn new(workers: usize) -> Self {
        assert!(workers > 0, "at least one native worker is required");
        let (sender, receiver) = mpsc::sync_channel::<Calculation>(workers);
        let receiver = Arc::new(Mutex::new(receiver));
        let threads = (0..workers)
            .map(|_| {
                let receiver = Arc::clone(&receiver);
                thread::spawn(move || {
                    loop {
                        let job = receiver.lock().expect("owned queue lock").recv();
                        let Ok(job) = job else { break };
                        let result = job.session.calc_ut(job.jd_ut, job.body, job.flags);
                        // A dropped receiver never cancels the native calculation.
                        let _ = job.reply.send(result);
                    }
                })
            })
            .collect();
        Self {
            sender: Some(sender),
            threads,
        }
    }

    /// Submission can block until the bounded queue has space. An async server
    /// uses its admission/executor adapter rather than blocking its event loop.
    fn submit(
        &self,
        session: Arc<Session>,
        jd_ut: f64,
        body: i32,
        flags: i32,
    ) -> mpsc::Receiver<Result<Position, Error>> {
        let (reply, result) = mpsc::sync_channel(1);
        self.sender
            .as_ref()
            .expect("open worker pool")
            .send(Calculation {
                session,
                jd_ut,
                body,
                flags,
                reply,
            })
            .unwrap_or_else(|_| panic!("live worker receiver"));
        result
    }
}
impl Drop for NativeWorkers {
    fn drop(&mut self) {
        // Disconnect admission, then finish both queued and running jobs.
        self.sender.take();
        for worker in self.threads.drain(..) {
            if worker.join().is_err() {
                eprintln!("native example worker panicked");
            }
        }
    }
}

fn main() -> Result<(), Box<dyn StdError>> {
    let image = std::env::args().nth(1).ok_or("provide an SPK path")?;
    let mut registry = SpkRegistry::new();
    registry.register_image(&image, std::fs::read(&image)?)?;
    register_all_naif_targets(&mut registry, &image)?;
    let mut context = EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_network_policy(Some(NetworkPolicy::Sealed))
        .with_auto_spk_download(Some(false))
        .with_delta_t_model(Some(DeltaTModel::Smh2016))
        .with_iers_delta_t(Some(false));
    *context.spk_registry_mut() = registry;
    let session = Arc::new(Session::from_context(context));
    let workers = NativeWorkers::new(4);
    let replies: Vec<_> = (0..8)
        .map(|day| {
            workers.submit(
                Arc::clone(&session),
                2451545.0 + f64::from(day),
                SUN,
                FLG_SWIEPH | FLG_SPEED,
            )
        })
        .collect();
    for reply in replies {
        let position = reply.recv()??;
        println!("{:?}, flags={}", position.values, position.returned_flags);
    }
    Ok(())
}
