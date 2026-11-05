//! Policy-gated native HTTP/TLS transport.
//!
//! [`HttpRequest`] validates HTTP/HTTPS URLs. [`TransportConfig`] supplies
//! connection and overall deadlines plus an in-memory body limit. [`fetch`]
//! checks the effective network policy before socket work and returns typed
//! URL, transport, status and size errors. Query/download callers own retries.
//!
//! TLS certificate verification and redirect handling use the configured
//! native client. Request-specific headers are retained; the default user
//! agent is added only when none is supplied. Large assets use the streaming
//! download layer rather than the bounded in-memory response.
//!
//! Timeouts are seconds, body limits are bytes, and this module interprets
//! no astronomical frame or epoch. Hermetic tests use loopback servers.

use std::io::{Read, Write};
use std::sync::{Arc, Mutex};

use std::time::Duration;
use ureq::ResponseExt as _;

use crate::context::{EffectiveNetworkPolicy, EphemerisContext, require_network};
use crate::errors::{Error, Result};

/// Default transport timeout per call site, in seconds.
///
/// Matches the release default (`timeout: int = 30`) of
/// `download_file` and the Horizons client (further call-site defaults
/// belong to their owning lanes). SPK generation (120 s) and the ASSIST bootstrap (300 s)
/// pass their own values via [`TransportConfig::with_timeout_secs`].
pub const DEFAULT_TIMEOUT_SECS: f64 = 30.0;

/// Default in-memory response-body bound for the transport.
///
/// 32 MiB comfortably holds the largest in-memory release payload (a
/// base64 SPK document inside a Horizons JSON response) while keeping a
/// hostile server from exhausting RAM. Bulk downloads stream to a file
/// instead (the download layer).
pub const DEFAULT_MAX_BODY_BYTES: usize = 32 * 1024 * 1024;

/// Default `User-Agent` header (`libephemeris-rust/<crate version>`).
pub const DEFAULT_USER_AGENT: &str = concat!("libephemeris-core/", env!("CARGO_PKG_VERSION"));

/// HTTP method for a transport [`HttpRequest`].
///
/// Only the methods the release call sites use are exposed. Every
/// release `open_url` call observed at the pin is a `GET`; `POST` is
/// carried because the transport surface (`docs/rust.md`) allows
/// it and a future Horizons query form may need it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HttpMethod {
    /// Retrieve a resource (all release call sites at the pin).
    Get,
    /// Submit a byte body (transport capability; no release caller yet).
    Post,
}

/// One policy-gated HTTP request for the transport.
///
/// Private fields with a checked constructor: the URL scheme/host rule is
/// enforced at construction, so a malformed URL never reaches a socket.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HttpRequest {
    url: String,
    method: HttpMethod,
    headers: Vec<(String, String)>,
    body: Vec<u8>,
}

impl HttpRequest {
    /// Build a `GET` request, validating the URL.
    pub fn get(url: &str) -> Result<Self> {
        validate_url(url)?;
        Ok(Self {
            url: url.to_string(),
            method: HttpMethod::Get,
            headers: Vec::new(),
            body: Vec::new(),
        })
    }

    /// Build a `POST` request with a byte body, validating the URL.
    pub fn post(url: &str, body: Vec<u8>) -> Result<Self> {
        validate_url(url)?;
        Ok(Self {
            url: url.to_string(),
            method: HttpMethod::Post,
            headers: Vec::new(),
            body,
        })
    }

    /// Append a request header (stored; sent in order).
    ///
    /// Values that the HTTP layer rejects surface at [`fetch`] time as
    /// [`Error::HttpTransport`]; ordinary release headers (ASCII
    /// `User-Agent`, `Content-Type`) always pass.
    pub fn header(mut self, name: &str, value: &str) -> Self {
        self.headers.push((name.to_string(), value.to_string()));
        self
    }

    /// Request URL.
    pub fn url(&self) -> &str {
        &self.url
    }

    /// Request method.
    pub fn method(&self) -> HttpMethod {
        self.method
    }

    /// Request headers in send order.
    pub fn headers(&self) -> &[(String, String)] {
        &self.headers
    }

    /// `POST` body bytes (empty for `GET`).
    pub fn body(&self) -> &[u8] {
        &self.body
    }
}

/// Accept exactly `http://host/...` and `https://host/...` with a
/// non-empty host. This is a deliberately small gate, not a full URI
/// parser: anything subtler (user-info, ports, escapes) passes through to
/// the HTTP layer, which applies its own strict parse before any I/O.
fn validate_url(url: &str) -> Result<()> {
    let rest = if let Some(rest) = url.strip_prefix("http://") {
        rest
    } else if let Some(rest) = url.strip_prefix("https://") {
        rest
    } else {
        return Err(Error::HttpInvalidUrl {
            url: url.to_string(),
        });
    };
    let host = rest.split(['/', '?', '#']).next().unwrap_or("");
    if host.is_empty() {
        return Err(Error::HttpInvalidUrl {
            url: url.to_string(),
        });
    }
    Ok(())
}

/// Explicit timeouts and bounds for one transport [`fetch`] call.
///
/// Private fields with checked constructors: a non-finite or
/// non-positive timeout is [`Error::InvalidTimeout`], a zero body bound
/// is [`Error::InvalidBodyLimit`]. Nothing is silently clamped.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportConfig {
    connect_timeout: Duration,
    overall_timeout: Duration,
    max_body_bytes: usize,
    user_agent: String,
}

impl TransportConfig {
    /// Default config: 30 s budgets, 32 MiB body bound,
    /// crate-version `User-Agent`.
    pub fn new() -> Result<Self> {
        Self::with_timeout_secs(DEFAULT_TIMEOUT_SECS)
    }

    /// Config with both time budgets set from one release-style `timeout`
    /// value (seconds), keeping the default body bound and `User-Agent`.
    ///
    /// Mirrors the release call sites that spell a single `timeout`
    /// (30/120/300 s): the connect budget covers socket setup plus the TLS
    /// handshake, the overall budget covers the whole exchange end to end
    /// (headers plus body).
    pub fn with_timeout_secs(timeout_secs: f64) -> Result<Self> {
        Ok(Self {
            connect_timeout: checked_timeout(timeout_secs)?,
            overall_timeout: checked_timeout(timeout_secs)?,
            max_body_bytes: DEFAULT_MAX_BODY_BYTES,
            user_agent: DEFAULT_USER_AGENT.to_string(),
        })
    }

    /// Override the connect (socket + TLS handshake) budget.
    pub fn with_connect_timeout_secs(mut self, secs: f64) -> Result<Self> {
        self.connect_timeout = checked_timeout(secs)?;
        Ok(self)
    }

    /// Override the overall end-to-end (headers plus body) budget.
    pub fn with_overall_timeout_secs(mut self, secs: f64) -> Result<Self> {
        self.overall_timeout = checked_timeout(secs)?;
        Ok(self)
    }

    /// Override the in-memory response-body bound (must be >= 1 byte).
    pub fn with_max_body_bytes(mut self, limit_bytes: usize) -> Result<Self> {
        if limit_bytes == 0 {
            return Err(Error::InvalidBodyLimit { limit_bytes });
        }
        self.max_body_bytes = limit_bytes;
        Ok(self)
    }

    /// Override the injected `User-Agent` value.
    pub fn with_user_agent(mut self, user_agent: &str) -> Self {
        self.user_agent = user_agent.to_string();
        self
    }

    /// Connect budget (socket setup plus TLS handshake).
    pub fn connect_timeout(&self) -> Duration {
        self.connect_timeout
    }

    /// Overall end-to-end budget (headers plus body), also reported in
    /// [`Error::HttpTimeout`].
    pub fn overall_timeout(&self) -> Duration {
        self.overall_timeout
    }

    /// Overall budget in seconds (for error reporting).
    pub fn overall_timeout_secs(&self) -> f64 {
        self.overall_timeout.as_secs_f64()
    }

    /// In-memory response-body bound in bytes.
    pub fn max_body_bytes(&self) -> usize {
        self.max_body_bytes
    }

    /// `User-Agent` injected when the request carries none.
    pub fn user_agent(&self) -> &str {
        &self.user_agent
    }
}

impl Default for TransportConfig {
    fn default() -> Self {
        Self::with_timeout_secs(DEFAULT_TIMEOUT_SECS).expect("default transport timeout is valid")
    }
}

/// Finite, positive seconds become a [`Duration`]; anything else is
/// [`Error::InvalidTimeout`]. `Duration::from_secs_f64` would panic on
/// non-finite/negative input, so the check precedes the conversion.
fn checked_timeout(secs: f64) -> Result<Duration> {
    if !secs.is_finite() || secs <= 0.0 {
        return Err(Error::InvalidTimeout { seconds: secs });
    }
    Ok(Duration::from_secs_f64(secs))
}

/// One bounded HTTP response from the transport.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HttpResponse {
    status: u16,
    headers: Vec<(String, String)>,
    body: Vec<u8>,
}

impl HttpResponse {
    /// Numeric status code (2xx here; other codes are [`Error::HttpStatus`]).
    pub fn status(&self) -> u16 {
        self.status
    }

    /// Response headers in received order (names lowercased).
    pub fn headers(&self) -> &[(String, String)] {
        &self.headers
    }

    /// Raw response-body bytes (within the configured bound).
    pub fn body(&self) -> &[u8] {
        &self.body
    }

    /// Response body as UTF-8 text (e.g. for the Horizons/JSON query callers).
    ///
    /// A non-UTF-8 body is [`Error::HttpTransport`], never lossy text:
    /// silently replacing bytes would corrupt API payloads.
    pub fn body_text(&self) -> Result<&str> {
        std::str::from_utf8(&self.body).map_err(|e| Error::HttpTransport {
            purpose: "decode response body as text".to_string(),
            detail: e.to_string(),
        })
    }
}

/// Perform one policy-gated HTTP exchange (the `open_url` equivalent).
///
/// The sealed-policy check runs first: under
/// [`EffectiveNetworkPolicy::Sealed`] this returns
/// [`Error::NetworkSealed`] without touching a socket. Otherwise exactly
/// one request/response exchange runs under `config` and the bounded body
/// is returned. No retry, no cache, no fallback — those are caller-side
/// (the query/download layers).
pub fn fetch(
    request: &HttpRequest,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<HttpResponse> {
    let mut response = open_response(request, purpose, effective, config, false, false)?;
    let limit = config.max_body_bytes;
    let mut body = Vec::new();
    let mut chunk = [0u8; 65536];
    loop {
        let remaining = limit.saturating_sub(body.len());
        let capacity = chunk.len().min(remaining.saturating_add(1));
        let count = response.read(&mut chunk[..capacity])?;
        if count == 0 {
            break;
        }
        if count > remaining {
            return Err(Error::HttpResponseTooLarge {
                purpose: purpose.to_string(),
                limit_bytes: limit,
            });
        }
        body.extend_from_slice(&chunk[..count]);
    }
    Ok(HttpResponse {
        status: response.status,
        headers: response.headers,
        body,
    })
}

/// Owned HTTP response whose body is read incrementally without a memory cap.
///
/// Headers and status are available before any body read. The caller provides
/// each destination buffer and owns its byte/progress/validation policy. The
/// connection and timeout state are owned by the response, which can move to a
/// worker thread. Closing or dropping it releases that resource once; no
/// astronomical data, callbacks, foreign objects or borrowed context are held.
pub struct StreamingResponse {
    status: u16,
    headers: Vec<(String, String)>,
    reader: Option<ureq::BodyReader<'static>>,
    purpose: String,
    config: TransportConfig,
    received_head: Option<ReceivedResponseHead>,
}

/// Original received HTTP status/header bytes and final response URL.
///
/// This metadata is available only from [`open_stream_with_status`]. It is
/// observed from the provider's actual plaintext header buffer before protocol
/// consumption, without reading body bytes or changing HTTP/TLS processing.
/// Within the provider's accepted HTTP domain, original header spelling,
/// duplicate fields and Latin-1 status reason bytes are retained. This metadata
/// does not widen that provider's accepted header or informational-status domain.
#[derive(Debug)]
pub struct ReceivedResponseHead {
    url: String,
    status_line: Vec<u8>,
    headers: Vec<u8>,
}

impl ReceivedResponseHead {
    /// Final response URL after the transport's ordinary redirect handling.
    pub fn url(&self) -> &str {
        &self.url
    }

    /// Original status line, including its terminating CRLF.
    pub fn status_line(&self) -> &[u8] {
        &self.status_line
    }

    /// Original header fields through the terminating blank line.
    pub fn header_bytes(&self) -> &[u8] {
        &self.headers
    }
}

impl StreamingResponse {
    /// Received HTTP status code. The checked opener rejects failure statuses.
    pub fn status(&self) -> u16 {
        self.status
    }

    /// Received response headers, with lowercased names.
    pub fn headers(&self) -> &[(String, String)] {
        &self.headers
    }

    /// Actual status/header bytes retained by [`open_stream_with_status`].
    ///
    /// Checked [`open_stream`] responses do not allocate this metadata.
    pub fn received_head(&self) -> Option<&ReceivedResponseHead> {
        self.received_head.as_ref()
    }

    /// Read at most `buffer.len` bytes, returning zero at the body end.
    ///
    /// Body I/O faults retain the exchange's purpose and timeout budget. A
    /// closed response refuses with `HttpTransport`; no partial body is
    /// substituted. An empty buffer consumes no bytes and returns zero.
    pub fn read(&mut self, buffer: &mut [u8]) -> Result<usize> {
        let reader = self.reader.as_mut().ok_or_else(|| Error::HttpTransport {
            purpose: self.purpose.clone(),
            detail: "response body is closed".to_owned(),
        })?;
        reader
            .read(buffer)
            .map_err(|error| map_body_error(&error, &self.purpose, &self.config))
    }

    /// Release the body connection. Repeated closes are harmless.
    pub fn close(&mut self) {
        self.reader.take();
    }
}

/// Open one policy-gated HTTP exchange and retain its body for streaming.
///
/// The policy check precedes socket work. Request headers, TLS, redirects,
/// status/error mapping match [`fetch`]. Connect timeout is unchanged, but
/// the overall timeout value becomes an independent timeout for each socket
/// read/write, as in a streaming socket response. A healthy bulk transfer may
/// take longer than that value in total. `max_body_bytes` applies only to
/// `fetch`. No body byte is read before returning headers;
/// the caller owns chunk size, progress, retries and durable promotion. A
/// caller that has already consumed its public policy guard supplies `Allow`.
pub fn open_stream(
    request: &HttpRequest,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<StreamingResponse> {
    open_response(request, purpose, effective, config, true, false)
}

/// Open a streaming response while retaining HTTP failure status and its body.
///
/// Policy, request, TLS, redirect and per-socket timeout behavior match
/// [`open_stream`]. HTTP status alone is not converted to [`Error::HttpStatus`]:
/// the caller receives the actual final URL, raw status/header metadata and the
/// same owned unread body for either success or failure. This supports genuine
/// public response/error owners that can read or explicitly close that body.
/// Transport/protocol faults still return typed errors. No retry, body pre-read,
/// status-message substitution or language-runtime dependency is introduced.
/// Informational responses, accepted header syntax and content decoding retain
/// the native HTTP provider's semantics. Adapters that must consume the actual
/// wire protocol with their own HTTP implementation use [`connect_response_wire`].
pub fn open_stream_with_status(
    request: &HttpRequest,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<StreamingResponse> {
    open_response(request, purpose, effective, config, true, true)
}

/// Select the timeout contract without changing request/TLS/body ownership.
/// Buffered fetch retains its total exchange budget; bulk streaming refreshes
/// the I/O budget at each actual transport read/write.
fn open_response(
    request: &HttpRequest,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    streaming: bool,
    retain_http_status: bool,
) -> Result<StreamingResponse> {
    require_network(effective, purpose)?;

    let observed_head = retain_http_status.then(|| Arc::new(Mutex::new(None)));
    let agent_config = ureq::config::Config::builder()
        .http_status_as_error(!retain_http_status)
        .timeout_connect(Some(config.connect_timeout))
        .timeout_global((!streaming).then_some(config.overall_timeout))
        .build();
    let agent = if streaming {
        ureq::Agent::with_parts(
            agent_config,
            SocketTimeoutConnector {
                inner: ureq::unversioned::transport::DefaultConnector::default(),
                timeout: config.overall_timeout,
                observed_head: observed_head.clone(),
            },
            ureq::unversioned::resolver::DefaultResolver::default(),
        )
    } else {
        ureq::Agent::new_with_config(agent_config)
    };

    let has_user_agent = request
        .headers
        .iter()
        .any(|(name, _)| name.eq_ignore_ascii_case("user-agent"));

    let response = match request.method {
        HttpMethod::Get => {
            let mut call = agent.get(request.url.as_str());
            if !has_user_agent {
                call = call.header("User-Agent", config.user_agent.as_str());
            }
            for (name, value) in &request.headers {
                call = call.header(name.as_str(), value.as_str());
            }
            call.call()
        }
        HttpMethod::Post => {
            let mut send = agent.post(request.url.as_str());
            if !has_user_agent {
                send = send.header("User-Agent", config.user_agent.as_str());
            }
            for (name, value) in &request.headers {
                send = send.header(name.as_str(), value.as_str());
            }
            send.send(request.body.as_slice())
        }
    };
    let response = response.map_err(|e| map_call_error(&e, purpose, config))?;

    let status = response.status().as_u16();
    let received_head = if let Some(observed) = observed_head {
        let raw = observed
            .lock()
            .map_err(|_| Error::HttpTransport {
                purpose: purpose.to_owned(),
                detail: "response header observation lock poisoned".to_owned(),
            })?
            .take()
            .ok_or_else(|| Error::HttpTransport {
                purpose: purpose.to_owned(),
                detail: "received response header bytes are unavailable".to_owned(),
            })?;
        // The provider has already validated the complete HTTP response. Only
        // split the retained bytes at its original status-line boundary.
        let line_end =
            raw.iter()
                .position(|&byte| byte == b'\n')
                .ok_or_else(|| Error::HttpTransport {
                    purpose: purpose.to_owned(),
                    detail: "received response status line is incomplete".to_owned(),
                })?
                + 1;
        Some(ReceivedResponseHead {
            url: response.get_uri().to_string(),
            status_line: raw[..line_end].to_vec(),
            headers: raw[line_end..].to_vec(),
        })
    } else {
        None
    };
    let headers: Vec<(String, String)> = response
        .headers()
        .iter()
        .map(|(name, value)| {
            (
                name.to_string(),
                String::from_utf8_lossy(value.as_bytes()).into_owned(),
            )
        })
        .collect();

    let (_, body) = response.into_parts();
    Ok(StreamingResponse {
        status,
        headers,
        reader: Some(body.into_reader()),
        purpose: purpose.to_owned(),
        config: config.clone(),
        received_head,
    })
}

/// One exchange's optional actual header publication; it owns no body bytes
/// or language-runtime values and is shared only with its native transports.
type HeadObservation = Arc<Mutex<Option<Vec<u8>>>>;

/// Default verified TCP/TLS connector with a bulk socket-I/O timeout adapter.
/// It observes optional received metadata without changing certificates, proxy
/// selection, buffers, redirects or pool state.
#[derive(Debug)]
struct SocketTimeoutConnector {
    inner: ureq::unversioned::transport::DefaultConnector,
    timeout: Duration,
    observed_head: Option<HeadObservation>,
}

impl ureq::unversioned::transport::Connector for SocketTimeoutConnector {
    type Out = SocketTimeoutTransport;

    /// Wrap only an actually connected transport, retaining normal refusals.
    fn connect(
        &self,
        details: &ureq::unversioned::transport::ConnectionDetails<'_>,
        chained: Option<()>,
    ) -> std::result::Result<Option<Self::Out>, ureq::Error> {
        self.inner.connect(details, chained).map(|transport| {
            transport.map(|inner| SocketTimeoutTransport {
                inner,
                timeout: self.timeout,
                observed_head: self.observed_head.clone(),
                awaiting_head: true,
            })
        })
    }
}

/// Owned connected transport. Each delegated I/O starts the same timeout anew;
/// cached buffered input consumes no socket wait. The wrapper supplies no data
/// or retry and retains the underlying transport's close/TLS behavior.
#[derive(Debug)]
struct SocketTimeoutTransport {
    inner: Box<dyn ureq::unversioned::transport::Transport>,
    timeout: Duration,
    observed_head: Option<HeadObservation>,
    awaiting_head: bool,
}

impl SocketTimeoutTransport {
    /// Observe complete header bytes before the provider consumes its buffer.
    ///
    /// This copies only the actual header framing; the provider remains the
    /// sole HTTP parser and validates status, headers, redirects and bodies.
    /// Interim responses retain the waiting state until the final response.
    fn observe_head(&mut self) {
        if !self.awaiting_head || self.observed_head.is_none() {
            return;
        }
        let input = self.inner.buffers().input();
        if !input.starts_with(b"HTTP/") {
            return;
        }
        let mut end = 0;
        let mut complete = false;
        for line in input.split_inclusive(|&byte| byte == b'\n') {
            end += line.len();
            if line == b"\r\n" || line == b"\n" {
                complete = true;
                break;
            }
        }
        if !complete {
            return;
        }
        let first = input.split(|&byte| byte == b'\n').next().unwrap_or(&[]);
        let status = first
            .split(|byte| byte.is_ascii_whitespace())
            .filter(|word| !word.is_empty())
            .nth(1)
            .unwrap_or(&[]);
        let interim = status.first() == Some(&b'1');
        if let Some(observed) = &self.observed_head {
            if let Ok(mut stored) = observed.lock() {
                *stored = Some(input[..end].to_vec());
            }
        }
        self.awaiting_head = interim;
    }

    /// Reuse the operation's timeout classification with a fresh I/O duration.
    fn deadline(
        &self,
        original: ureq::unversioned::transport::NextTimeout,
    ) -> ureq::unversioned::transport::NextTimeout {
        ureq::unversioned::transport::NextTimeout {
            after: self.timeout.into(),
            reason: original.reason,
        }
    }
}

impl ureq::unversioned::transport::Transport for SocketTimeoutTransport {
    /// Observe requested header metadata, then borrow the actual provider buffers.
    /// No body byte is copied or consumed by this observation.
    fn buffers(&mut self) -> &mut dyn ureq::unversioned::transport::Buffers {
        self.observe_head();
        self.inner.buffers()
    }

    /// Send the acquired request bytes under one fresh write timeout.
    fn transmit_output(
        &mut self,
        amount: usize,
        timeout: ureq::unversioned::transport::NextTimeout,
    ) -> std::result::Result<(), ureq::Error> {
        self.awaiting_head = true;
        self.inner.transmit_output(amount, self.deadline(timeout))
    }

    /// Receive under one fresh read timeout, preserving provider progress.
    fn await_input(
        &mut self,
        timeout: ureq::unversioned::transport::NextTimeout,
    ) -> std::result::Result<bool, ureq::Error> {
        self.inner.await_input(self.deadline(timeout))
    }

    /// Delegate connection health and reuse without extending its lifetime.
    fn is_open(&mut self) -> bool {
        self.inner.is_open()
    }

    /// Preserve the provider's actual verified TLS state.
    fn is_tls(&self) -> bool {
        self.inner.is_tls()
    }
}

/// An owned verified connection whose HTTP bytes have not been interpreted.
///
/// The native connector owns TCP, proxy and verified TLS setup. The caller owns
/// request serialization and response framing. Reads preserve every received
/// byte, including informational status lines, folded fields and compressed
/// bodies; no HTTP parser, redirect or decompressor runs on this connection.
/// The value owns its native resources and may move to a blocking worker.
pub struct ResponseWire {
    transport: Option<Box<dyn ureq::unversioned::transport::Transport>>,
    purpose: String,
    config: TransportConfig,
    deadline: Option<Duration>,
}

impl ResponseWire {
    /// Receive at most the requested bytes without HTTP/body transformation.
    /// Each blocking receive gets the configured independent socket-I/O budget.
    /// A closed owner refuses; an empty destination consumes no socket bytes.
    pub fn read(&mut self, destination: &mut [u8]) -> Result<usize> {
        let transport = self
            .transport
            .as_mut()
            .ok_or_else(|| Error::HttpTransport {
                purpose: self.purpose.clone(),
                detail: "response wire is closed".to_owned(),
            })?;
        if destination.is_empty() {
            return Ok(0);
        }
        loop {
            let buffered = transport.buffers().input();
            if !buffered.is_empty() {
                let count = destination.len().min(buffered.len());
                destination[..count].copy_from_slice(&buffered[..count]);
                transport.buffers().input_consume(count);
                return Ok(count);
            }
            let progress = transport
                .await_input(ureq::unversioned::transport::NextTimeout {
                    after: wire_deadline(self.deadline),
                    reason: ureq::Timeout::RecvBody,
                })
                .map_err(|error| map_call_error(&error, &self.purpose, &self.config))?;
            if !progress {
                return Ok(0);
            }
        }
    }

    /// Send all caller-serialized request bytes under native socket timeouts.
    /// The provider sends each filled output buffer without adding HTTP fields,
    /// changing method/body data or interpreting a received response.
    pub fn write_all(&mut self, mut bytes: &[u8]) -> Result<()> {
        let transport = self
            .transport
            .as_mut()
            .ok_or_else(|| Error::HttpTransport {
                purpose: self.purpose.clone(),
                detail: "response wire is closed".to_owned(),
            })?;
        while !bytes.is_empty() {
            let output = transport.buffers().output();
            let count = output.len().min(bytes.len());
            if count == 0 {
                return Err(Error::HttpTransport {
                    purpose: self.purpose.clone(),
                    detail: "response wire has no output capacity".to_owned(),
                });
            }
            output[..count].copy_from_slice(&bytes[..count]);
            transport
                .transmit_output(
                    count,
                    ureq::unversioned::transport::NextTimeout {
                        after: wire_deadline(self.deadline),
                        reason: ureq::Timeout::SendBody,
                    },
                )
                .map_err(|error| map_call_error(&error, &self.purpose, &self.config))?;
            bytes = &bytes[count..];
        }
        Ok(())
    }

    /// Release the actual native connection once. Repeated close is harmless.
    pub fn close(&mut self) {
        self.transport.take();
    }
}

/// Connect a policy-gated native HTTP wire without sending a request.
/// The caller owns response framing. Native proxy/TLS selection, deadlines
/// and connection failures use the shared transport configuration.
pub fn connect_response_wire(
    url: &str,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<ResponseWire> {
    connect_response_wire_inner(
        url,
        purpose,
        effective,
        config,
        Some(config.connect_timeout),
        Some(config.overall_timeout),
    )
}

/// Connect a native wire with an optional independent socket deadline.
///
/// `None` keeps blocking connect/read/write operations unbounded. A finite
/// duration supplies the deadline for each operation independently. This
/// additive acquisition boundary does not change the positive, finite checked
/// contract of [`TransportConfig`], [`fetch`] or [`open_stream`]. The caller
/// owns language-specific timeout conversion before supplying a duration.
/// Policy and HTTP(S) URL admission still precede connection acquisition.
pub fn connect_response_wire_with_deadline(
    url: &str,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    deadline: Option<Duration>,
) -> Result<ResponseWire> {
    let mut config = TransportConfig::default();
    if let Some(deadline) = deadline {
        config.connect_timeout = deadline;
        config.overall_timeout = deadline;
    }
    connect_response_wire_inner(url, purpose, effective, &config, deadline, deadline)
}

/// Share connection ownership without substituting an unbounded deadline with
/// a large numeric duration. The checked entry retains distinct connect/I/O
/// budgets; the acquired entry supplies its actual optional socket budget.
fn connect_response_wire_inner(
    url: &str,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    connect_deadline: Option<Duration>,
    deadline: Option<Duration>,
) -> Result<ResponseWire> {
    require_network(effective, purpose)?;
    validate_url(url)?;
    use ureq::unversioned::transport::Connector;
    let inner: WireInner = if let Some(timeout) = deadline {
        Box::new(
            SocketTimeoutConnector {
                inner: ureq::unversioned::transport::DefaultConnector::default(),
                timeout,
                observed_head: None,
            }
            .chain(WireBoxConnector),
        )
    } else {
        Box::new(ureq::unversioned::transport::DefaultConnector::default())
    };
    connect_response_wire_parts(
        url,
        purpose,
        config,
        connect_deadline,
        deadline,
        inner,
        ureq::unversioned::resolver::DefaultResolver::default(),
    )
}

/// Transfer the connection from the selected native connector before the
/// provider serializes any HTTP request. The resolver and connector are paired.
fn connect_response_wire_parts(
    url: &str,
    purpose: &str,
    config: &TransportConfig,
    connect_deadline: Option<Duration>,
    deadline: Option<Duration>,
    inner: WireInner,
    resolver: impl ureq::unversioned::resolver::Resolver,
) -> Result<ResponseWire> {
    let handoff = Arc::new(Mutex::new(None));
    let agent = ureq::Agent::with_parts(
        ureq::config::Config::builder()
            .timeout_connect(connect_deadline)
            .build(),
        WireConnector {
            inner,
            handoff: Arc::clone(&handoff),
            connecting: std::sync::atomic::AtomicUsize::new(0),
        },
        resolver,
    );
    // The connector transfers its actual verified connection before returning
    // control to the provider. The provider therefore serializes/sends nothing.
    let result = agent.get(url).call();
    let transport = handoff
        .lock()
        .map_err(|_| Error::HttpTransport {
            purpose: purpose.to_owned(),
            detail: "response wire ownership lock poisoned".to_owned(),
        })?
        .take();
    match transport {
        Some(transport) => Ok(ResponseWire {
            transport: Some(transport),
            purpose: purpose.to_owned(),
            config: config.clone(),
            deadline,
        }),
        None => Err(match result {
            Err(error) => map_call_error(&error, purpose, config),
            Ok(_) => Error::HttpTransport {
                purpose: purpose.to_owned(),
                detail: "native response wire was not acquired".to_owned(),
            },
        }),
    }
}

/// An actually acquired TCP connection and its socket timeout mode.
///
/// This owner contains no language-runtime values. `None` means blocking I/O;
/// a zero duration means nonblocking I/O; positive values are independent
/// receive/send budgets. Closing the owner closes the one owned socket.
pub struct AcquiredResponseTcp {
    stream: std::net::TcpStream,
    deadline: Option<Duration>,
}

impl AcquiredResponseTcp {
    /// Take an exclusively owned TCP socket after the caller's real connection
    /// primitives have completed. Normalize the descriptor's blocking flag for
    /// Rust I/O without changing the acquired timeout mode or socket identity.
    pub fn new(stream: std::net::TcpStream, deadline: Option<Duration>) -> std::io::Result<Self> {
        stream.set_nonblocking(deadline == Some(Duration::ZERO))?;
        Ok(Self { stream, deadline })
    }
}

/// Call-owned socket acquisition at the native TCP connector boundary.
///
/// Implementations resolve/connect the actual host and port and transfer an
/// exclusively owned socket. All raw operand conversion and refused operations
/// belong to the implementation. The native chain retains proxy tunneling,
/// verified TLS and subsequent resource ownership. No callback is retained by
/// a returned wire after connection acquisition completes.
pub trait ResponseTcpAcquisition: std::fmt::Debug + Send + Sync + 'static {
    /// Acquire the actual TCP connection, preserving the implementation's DNS,
    /// timeout-conversion, address fallback and connection-error order.
    fn connect(&self, host: &str, port: u16) -> std::io::Result<AcquiredResponseTcp>;
}

/// Connect a wire using a call-owned TCP acquisition, then native TLS/proxy
/// processing. This entry performs no preliminary DNS lookup or timeout
/// conversion: the acquired connector consumes those actual operations once.
/// Checked native openers and their resolver/timeout contracts are unchanged.
pub fn connect_response_wire_acquired(
    url: &str,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    acquisition: Arc<dyn ResponseTcpAcquisition>,
) -> Result<ResponseWire> {
    use ureq::unversioned::transport::Connector;
    require_network(effective, purpose)?;
    validate_url(url)?;
    let inner =
        ().chain(ureq::unversioned::transport::ConnectProxyConnector::default())
            .chain(AcquiredTcpConnector { acquisition })
            .chain(ureq::unversioned::transport::RustlsConnector::default())
            .chain(WireBoxConnector);
    connect_response_wire_parts(
        url,
        purpose,
        &TransportConfig::default(),
        None,
        None,
        Box::new(inner),
        AcquiredTcpResolver,
    )
}

/// DNS is deliberately deferred to the acquired TCP operation. There is no
/// resolved-address substitute: the paired connector consumes the actual URI
/// and performs its own resolution, as a remote-resolving proxy connector does.
#[derive(Debug)]
struct AcquiredTcpResolver;
impl ureq::unversioned::resolver::Resolver for AcquiredTcpResolver {
    fn resolve(
        &self,
        _: &ureq::http::Uri,
        _: &ureq::config::Config,
        _: ureq::unversioned::transport::NextTimeout,
    ) -> std::result::Result<ureq::unversioned::resolver::ResolvedSocketAddrs, ureq::Error> {
        Ok(self.empty())
    }
}

/// Pass a completed proxy tunnel through unchanged, otherwise acquire exactly
/// the underlying TCP host/port requested by the native chain.
#[derive(Debug)]
struct AcquiredTcpConnector {
    acquisition: Arc<dyn ResponseTcpAcquisition>,
}
impl<In: ureq::unversioned::transport::Transport> ureq::unversioned::transport::Connector<In>
    for AcquiredTcpConnector
{
    type Out = ureq::unversioned::transport::Either<In, AcquiredTcpTransport>;
    fn connect(
        &self,
        details: &ureq::unversioned::transport::ConnectionDetails<'_>,
        chained: Option<In>,
    ) -> std::result::Result<Option<Self::Out>, ureq::Error> {
        if let Some(transport) = chained {
            return Ok(Some(ureq::unversioned::transport::Either::A(transport)));
        }
        let host = details
            .uri
            .host()
            .ok_or_else(|| ureq::Error::BadUri("missing TCP host".to_owned()))?;
        let host = host
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .unwrap_or(host);
        let port = details
            .uri
            .port_u16()
            .unwrap_or(if details.needs_tls() { 443 } else { 80 });
        let acquired = self.acquisition.connect(host, port)?;
        let buffers = ureq::unversioned::transport::LazyBuffers::new(
            details.config.input_buffer_size(),
            details.config.output_buffer_size(),
        );
        Ok(Some(ureq::unversioned::transport::Either::B(
            AcquiredTcpTransport { acquired, buffers },
        )))
    }
}

/// Erase only connector types; the actual native transport owner is retained.
#[derive(Debug)]
struct WireBoxConnector;
impl<In: ureq::unversioned::transport::Transport> ureq::unversioned::transport::Connector<In>
    for WireBoxConnector
{
    type Out = Box<dyn ureq::unversioned::transport::Transport>;
    fn connect(
        &self,
        _: &ureq::unversioned::transport::ConnectionDetails<'_>,
        chained: Option<In>,
    ) -> std::result::Result<Option<Self::Out>, ureq::Error> {
        Ok(chained.map(|transport| Box::new(transport) as Self::Out))
    }
}

/// Native buffered I/O over the actual transferred socket. Its acquired mode
/// survives TLS framing and each subsequent request/response operation.
struct AcquiredTcpTransport {
    acquired: AcquiredResponseTcp,
    buffers: ureq::unversioned::transport::LazyBuffers,
}
impl std::fmt::Debug for AcquiredTcpTransport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AcquiredTcpTransport")
            .finish_non_exhaustive()
    }
}
impl ureq::unversioned::transport::Transport for AcquiredTcpTransport {
    fn buffers(&mut self) -> &mut dyn ureq::unversioned::transport::Buffers {
        &mut self.buffers
    }
    fn transmit_output(
        &mut self,
        amount: usize,
        timeout: ureq::unversioned::transport::NextTimeout,
    ) -> std::result::Result<(), ureq::Error> {
        use ureq::unversioned::transport::Buffers;
        let deadline = self.acquired.deadline.filter(|d| !d.is_zero());
        self.acquired.stream.set_write_timeout(deadline)?;
        self.acquired
            .stream
            .write_all(&self.buffers.output()[..amount])
            .map_err(|error| acquired_io_error(error, self.acquired.deadline, timeout.reason))
    }
    fn await_input(
        &mut self,
        timeout: ureq::unversioned::transport::NextTimeout,
    ) -> std::result::Result<bool, ureq::Error> {
        use ureq::unversioned::transport::Buffers;
        let deadline = self.acquired.deadline.filter(|d| !d.is_zero());
        self.acquired.stream.set_read_timeout(deadline)?;
        let count = self
            .acquired
            .stream
            .read(self.buffers.input_append_buf())
            .map_err(|error| acquired_io_error(error, self.acquired.deadline, timeout.reason))?;
        self.buffers.input_appended(count);
        Ok(count != 0)
    }
    fn is_open(&mut self) -> bool {
        true
    }
}

/// A positive acquired socket budget maps an elapsed socket wait to a native
/// timeout; nonblocking WouldBlock remains the actual I/O error.
fn acquired_io_error(
    error: std::io::Error,
    deadline: Option<Duration>,
    reason: ureq::Timeout,
) -> ureq::Error {
    if deadline.is_some_and(|d| !d.is_zero())
        && matches!(
            error.kind(),
            std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut
        )
    {
        ureq::Error::Timeout(reason)
    } else {
        error.into()
    }
}

/// Represent no deadline with the provider's explicit unbounded value.
fn wire_deadline(deadline: Option<Duration>) -> ureq::unversioned::transport::time::Duration {
    deadline.map_or(
        ureq::unversioned::transport::time::Duration::NotHappening,
        Into::into,
    )
}

/// Native connector type erasure retains the same owned transport allocation.
type WireInner = Box<
    dyn ureq::unversioned::transport::Connector<
            (),
            Out = Box<dyn ureq::unversioned::transport::Transport>,
        >,
>;

/// Per-call transport transfer; it owns native resources only, never callbacks.
type WireHandoff = Arc<Mutex<Option<Box<dyn ureq::unversioned::transport::Transport>>>>;

/// Private control transfer recognized only together with the actual wire owner.
/// It is never exposed as a public error or used as a replacement response.
#[derive(Debug)]
struct WireTransferred;
impl std::fmt::Display for WireTransferred {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("native response wire ownership transferred")
    }
}
impl std::error::Error for WireTransferred {}

/// Acquire native verification/timeout state, then transfer before HTTP I/O.
#[derive(Debug)]
struct WireConnector {
    inner: WireInner,
    handoff: WireHandoff,
    connecting: std::sync::atomic::AtomicUsize,
}

/// One nested connector invocation, including the provider's CONNECT proxy
/// setup. Only the outer origin connection may leave the connector chain.
struct WireConnectScope<'a>(&'a std::sync::atomic::AtomicUsize);
impl Drop for WireConnectScope<'_> {
    fn drop(&mut self) {
        self.0.fetch_sub(1, std::sync::atomic::Ordering::AcqRel);
    }
}

impl ureq::unversioned::transport::Connector for WireConnector {
    type Out = Box<dyn ureq::unversioned::transport::Transport>;
    fn connect(
        &self,
        details: &ureq::unversioned::transport::ConnectionDetails<'_>,
        chained: Option<()>,
    ) -> std::result::Result<Option<Self::Out>, ureq::Error> {
        // This agent performs one connect operation. Its native CONNECT proxy
        // recursively invokes the same connector to acquire the proxy socket;
        // that intermediate owner must reach the provider's tunnel/TLS setup.
        let depth = self
            .connecting
            .fetch_add(1, std::sync::atomic::Ordering::AcqRel);
        let _scope = WireConnectScope(&self.connecting);
        let transport = self.inner.connect(details, chained)?;
        if depth != 0 {
            return Ok(transport);
        }
        if let Some(transport) = transport {
            *self
                .handoff
                .lock()
                .map_err(|_| ureq::Error::Other(Box::new(WireTransferred)))? = Some(transport);
            return Err(ureq::Error::Other(Box::new(WireTransferred)));
        }
        Ok(None)
    }
}

/// Perform one policy-gated HTTP exchange under a context's effective
/// network policy (sealed-LEB contexts refuse before any socket work).
pub fn fetch_with_context(
    ctx: &EphemerisContext,
    request: &HttpRequest,
    purpose: &str,
    config: &TransportConfig,
) -> Result<HttpResponse> {
    fetch(request, purpose, ctx.effective_network_policy(), config)
}

/// Map a `ureq` call-phase error to the typed taxonomy.
///
/// `StatusCode` keeps its numeric code (the `HTTPError` equivalent);
/// `Timeout` keeps the configured overall budget; the size guard is
/// enforced separately in [`fetch`] so `BodyExceedsLimit` cannot arrive
/// here. Everything else (DNS, refused connections, TLS verification,
/// protocol faults, redirect refusal) is `HttpTransport` with a short
/// provider message — the `URLError`/`HTTPException`/`OSError` surface
/// the release callers catch.
///
/// Shared `pub(crate)` with the streaming downloader
/// ([`crate::download`]), which performs its own exchange because the
/// bounded in-memory [`fetch`] cannot hold multi-hundred-megabyte assets.
pub(crate) fn map_call_error(e: &ureq::Error, purpose: &str, config: &TransportConfig) -> Error {
    match e {
        ureq::Error::StatusCode(status) => Error::HttpStatus {
            status: *status,
            purpose: purpose.to_string(),
        },
        ureq::Error::Timeout(_) => Error::HttpTimeout {
            purpose: purpose.to_string(),
            seconds: config.overall_timeout_secs(),
        },
        other => Error::HttpTransport {
            purpose: purpose.to_string(),
            detail: short_detail(other),
        },
    }
}

/// Map a body-phase I/O error. A wait that outlasts the budget is a
/// timeout; anything else is a transport failure. No response payload is
/// ever embedded in the detail. Shared `pub(crate)` with [`crate::download`].
pub(crate) fn map_body_error(e: &std::io::Error, purpose: &str, config: &TransportConfig) -> Error {
    let typed_timeout = e
        .get_ref()
        .and_then(|cause| cause.downcast_ref::<ureq::Error>())
        .is_some_and(|cause| matches!(cause, ureq::Error::Timeout(_)));
    if e.kind() == std::io::ErrorKind::TimedOut || typed_timeout {
        Error::HttpTimeout {
            purpose: purpose.to_string(),
            seconds: config.overall_timeout_secs(),
        }
    } else {
        Error::HttpTransport {
            purpose: purpose.to_string(),
            detail: e.to_string(),
        }
    }
}

/// One-line provider message for the `HttpTransport` detail field.
fn short_detail(e: &ureq::Error) -> String {
    let text = e.to_string();
    match text.find('\n') {
        Some(idx) => text[..idx].to_string(),
        None => text,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{NetworkPolicy, ProviderMode};
    use std::net::TcpListener;
    use std::time::Duration;

    /// One-shot loopback HTTP server: serves a single connection with a
    /// canned status/body (after an optional delay), records the raw
    /// request head plus body for assertions, then exits.
    struct Fixture {
        url: String,
        received: std::sync::mpsc::Receiver<Vec<u8>>,
        _thread: std::thread::JoinHandle<()>,
    }

    /// Serve one connection: read the request head (and any declared body),
    /// optionally sleep, then write `status_line` + body.
    fn serve_once(status_line: &str, body: &[u8], delay: Duration, extra_headers: &str) -> Fixture {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind loopback");
        let addr = listener.local_addr().expect("loopback addr");
        let (tx, rx) = std::sync::mpsc::channel();
        let status_line = status_line.to_string();
        let body = body.to_vec();
        let extra_headers = extra_headers.to_string();
        let thread = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("accept");
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .expect("read timeout");
            let mut raw = Vec::new();
            let mut one = [0u8; 1];
            // Read until the end of the request head.
            loop {
                match stream.read(&mut one) {
                    Ok(0) => break,
                    Ok(_) => {
                        raw.push(one[0]);
                        if raw.ends_with(b"\r\n\r\n") {
                            break;
                        }
                        if raw.len() > 65536 {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            // Drain a declared request body so POST framing is complete.
            let head = String::from_utf8_lossy(&raw).into_owned();
            let content_length = head
                .lines()
                .find_map(|line| {
                    let (name, value) = line.split_once(':')?;
                    if name.trim().eq_ignore_ascii_case("content-length") {
                        value.trim().parse::<usize>().ok()
                    } else {
                        None
                    }
                })
                .unwrap_or(0);
            let mut rest = vec![0u8; content_length];
            if content_length > 0 {
                let _ = stream.read_exact(&mut rest);
                raw.extend_from_slice(&rest);
            }
            let _ = tx.send(raw);
            if !delay.is_zero() {
                std::thread::sleep(delay);
            }
            let head = format!(
                "{status_line}\r\nContent-Length: {}\r\n{extra_headers}\r\n",
                body.len()
            );
            let _ = stream.write_all(head.as_bytes());
            let _ = stream.write_all(&body);
            let _ = stream.flush();
        });
        Fixture {
            url: format!("http://{addr}/vectors"),
            received: rx,
            _thread: thread,
        }
    }

    fn recv_head(fixture: &Fixture) -> String {
        let raw = fixture
            .received
            .recv_timeout(Duration::from_secs(10))
            .expect("server saw the request");
        String::from_utf8_lossy(&raw).into_owned()
    }

    fn allow() -> EffectiveNetworkPolicy {
        EffectiveNetworkPolicy::Allow
    }

    #[test]
    fn status_preserving_stream_keeps_actual_head_and_body_without_message_substitution() {
        for (code, reason) in [
            (200, "Noncanonical success"),
            (503, "Maintenance: actual phrase"),
        ] {
            let fixture = serve_once(
                &format!("HTTP/1.1 {code} {reason}"),
                b"body!",
                Duration::ZERO,
                "X-MiXeD: first\r\nX-MiXeD: second\r\nConnection: close\r\n",
            );
            let request = HttpRequest::get(&fixture.url).unwrap();
            let mut response = open_stream_with_status(
                &request,
                "retained actual response",
                allow(),
                &TransportConfig::default(),
            )
            .unwrap();
            assert_eq!(response.status(), code);
            let head = response.received_head().unwrap();
            assert_eq!(head.url(), fixture.url);
            assert_eq!(
                head.status_line(),
                format!("HTTP/1.1 {code} {reason}\r\n").as_bytes()
            );
            assert_eq!(head.header_bytes(), b"Content-Length: 5\r\nX-MiXeD: first\r\nX-MiXeD: second\r\nConnection: close\r\n\r\n");
            let mut body = Vec::new();
            let mut buffer = [0; 2];
            loop {
                let count = response.read(&mut buffer).unwrap();
                if count == 0 {
                    break;
                }
                body.extend_from_slice(&buffer[..count]);
            }
            assert_eq!(body, b"body!");
            response.close();
            response.close();
            assert!(matches!(
                response.read(&mut buffer),
                Err(Error::HttpTransport { .. })
            ));
            assert!(recv_head(&fixture).starts_with("GET /vectors HTTP/1.1"));
            fixture._thread.join().unwrap();
        }
        let fixture = serve_once(
            "HTTP/1.1 503 Original checked refusal",
            b"unused",
            Duration::ZERO,
            "Connection: close\r\n",
        );
        assert!(matches!(
            open_stream(
                &HttpRequest::get(&fixture.url).unwrap(),
                "checked",
                allow(),
                &TransportConfig::default()
            ),
            Err(Error::HttpStatus { status: 503, .. })
        ));
        fixture._thread.join().unwrap();
    }

    #[test]
    fn retained_header_owner_follows_redirect_and_final_interim_buffer_consumption() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let server = std::thread::spawn(move || {
            for path in ["/initial", "/final"] {
                let (mut stream, _) = listener.accept().unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(10)))
                    .unwrap();
                let mut request = Vec::new();
                let mut byte = [0; 1];
                while !request.ends_with(b"\r\n\r\n") {
                    assert_eq!(stream.read(&mut byte).unwrap(), 1);
                    request.push(byte[0]);
                }
                assert!(
                    String::from_utf8(request)
                        .unwrap()
                        .starts_with(&format!("GET {path} HTTP/1.1"))
                );
                if path == "/initial" {
                    stream.write_all(b"HTTP/1.1 302 Actual redirect\r\nLocation: /final\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
                } else {
                    stream.write_all(b"HTTP/1.1 100 Continue\r\n\r\nHTTP/1.1 410 Final custom reason\r\nX-Final: yes\r\nContent-Length: 4\r\nConnection: close\r\n\r\nlast").unwrap();
                }
            }
        });
        let url = format!("http://{address}/initial");
        let mut response = open_stream_with_status(
            &HttpRequest::get(&url).unwrap(),
            "redirect owner",
            allow(),
            &TransportConfig::default(),
        )
        .unwrap();
        assert_eq!(response.status(), 410);
        let head = response.received_head().unwrap();
        assert_eq!(head.url(), format!("http://{address}/final"));
        assert_eq!(head.status_line(), b"HTTP/1.1 410 Final custom reason\r\n");
        assert!(head.header_bytes().starts_with(b"X-Final: yes\r\n"));
        let mut body = [0; 4];
        assert_eq!(response.read(&mut body).unwrap(), 4);
        assert_eq!(&body, b"last");
        server.join().unwrap();
    }

    #[test]
    fn status_owner_is_returned_before_the_server_provides_any_body_byte() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let (release, wait) = std::sync::mpsc::channel();
        let server = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .unwrap();
            let mut request = Vec::new();
            let mut byte = [0; 1];
            while !request.ends_with(b"\r\n\r\n") {
                assert_eq!(stream.read(&mut byte).unwrap(), 1);
                request.push(byte[0]);
            }
            stream.write_all(b"HTTP/1.1 429 Wait for caller\r\nContent-Length: 4\r\nConnection: close\r\n\r\n").unwrap();
            wait.recv_timeout(Duration::from_secs(10)).unwrap();
            stream.write_all(b"late").unwrap();
        });
        let mut response = open_stream_with_status(
            &HttpRequest::get(&format!("http://{address}/late")).unwrap(),
            "unread body",
            allow(),
            &TransportConfig::default(),
        )
        .unwrap();
        assert_eq!(response.status(), 429);
        assert_eq!(
            response.received_head().unwrap().status_line(),
            b"HTTP/1.1 429 Wait for caller\r\n"
        );
        release.send(()).unwrap();
        let mut bytes = [0; 4];
        let mut count = 0;
        while count < bytes.len() {
            let read = response.read(&mut bytes[count..]).unwrap();
            assert_ne!(read, 0);
            count += read;
        }
        assert_eq!(&bytes, b"late");
        server.join().unwrap();
    }

    #[test]
    fn acquired_wire_sends_only_caller_bytes_and_keeps_all_response_bytes() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let request = b"GET /caller-only HTTP/1.1\r\nHost: owned\r\nConnection: close\r\n\r\n";
        let reply=b"HTTP/1.1 103 Actual first\r\nX-Fold: one\r\n two\r\n\r\nHTTP/1.1 200 Later\r\nContent-Encoding: gzip\r\nContent-Length: 4\r\n\r\nwire";
        let server = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .unwrap();
            let mut actual = vec![0; request.len()];
            stream.read_exact(&mut actual).unwrap();
            assert_eq!(&actual, request);
            stream.write_all(reply).unwrap();
        });
        let mut wire = connect_response_wire(
            &format!("http://{address}/unused-agent-selector"),
            "actual wire",
            allow(),
            &TransportConfig::default(),
        )
        .unwrap();
        wire.write_all(request).unwrap();
        let mut bytes = Vec::new();
        let mut chunk = [0; 17];
        loop {
            let count = wire.read(&mut chunk).unwrap();
            if count == 0 {
                break;
            }
            bytes.extend_from_slice(&chunk[..count]);
        }
        assert_eq!(bytes, reply);
        wire.close();
        wire.close();
        assert!(matches!(
            wire.read(&mut chunk),
            Err(Error::HttpTransport { .. })
        ));
        server.join().unwrap();
    }

    #[test]
    fn url_gate_rejects_before_any_socket_work() {
        for bad in [
            "",
            "ftp://example.com/x",
            "file:///etc/hosts",
            "https://",
            "http://",
            "ssd.jpl.nasa.gov/api",
        ] {
            match HttpRequest::get(bad) {
                Err(Error::HttpInvalidUrl { url }) => assert_eq!(url, bad),
                other => panic!("{bad:?} must be HttpInvalidUrl, got {other:?}"),
            }
            match HttpRequest::post(bad, vec![1]) {
                Err(Error::HttpInvalidUrl { .. }) => {}
                other => panic!("{bad:?} POST must be HttpInvalidUrl, got {other:?}"),
            }
        }
        // Ordinary http/https shapes (host, port, IPv4/IPv6 loopback,
        // query) pass the gate.
        for good in [
            "https://ssd.jpl.nasa.gov/api/horizons.api?format=json",
            "http://127.0.0.1:9/vectors",
            "http://[::1]:8080/x",
            "https://example.com",
        ] {
            assert!(HttpRequest::get(good).is_ok(), "{good:?} must pass");
        }
    }

    #[test]
    fn timeout_and_limit_configs_are_checked() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            match TransportConfig::with_timeout_secs(bad) {
                Err(Error::InvalidTimeout { seconds }) => {
                    if bad.is_nan() {
                        assert!(seconds.is_nan(), "NaN timeout must round-trip");
                    } else {
                        assert_eq!(seconds, bad);
                    }
                }
                other => panic!("{bad} must be InvalidTimeout, got {other:?}"),
            }
        }
        match TransportConfig::new().unwrap().with_max_body_bytes(0) {
            Err(Error::InvalidBodyLimit { limit_bytes }) => assert_eq!(limit_bytes, 0),
            other => panic!("zero bound must be InvalidBodyLimit, got {other:?}"),
        }
        let cfg = TransportConfig::with_timeout_secs(120.0).unwrap();
        assert_eq!(cfg.connect_timeout(), Duration::from_secs(120));
        assert_eq!(cfg.overall_timeout(), Duration::from_secs(120));
        assert_eq!(cfg.max_body_bytes(), DEFAULT_MAX_BODY_BYTES);
        assert!(!cfg.user_agent().is_empty());
    }

    #[test]
    fn sealed_policy_refuses_before_socket_work() {
        // An unroutable port proves no connection was attempted: the sealed
        // gate fires first.
        let req = HttpRequest::get("http://127.0.0.1:9/vectors").unwrap();
        let cfg = TransportConfig::new().unwrap();
        match fetch(
            &req,
            "Horizons state vector for 499",
            EffectiveNetworkPolicy::Sealed,
            &cfg,
        ) {
            Err(Error::NetworkSealed { purpose }) => {
                assert_eq!(purpose, "Horizons state vector for 499");
            }
            other => panic!("sealed fetch must refuse, got {other:?}"),
        }
        // Same through the context wrapper: an explicit sealed-LEB context
        // refuses, an explicit-allow context proceeds to the (failing) socket.
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert_eq!(
            sealed.effective_network_policy(),
            EffectiveNetworkPolicy::Sealed
        );
        assert!(matches!(
            fetch_with_context(&sealed, &req, "probe", &cfg),
            Err(Error::NetworkSealed { .. })
        ));
        let allowed = EphemerisContext::new().with_network_policy(Some(NetworkPolicy::Allow));
        match fetch_with_context(&allowed, &req, "probe", &cfg) {
            Err(Error::NetworkSealed { .. }) => panic!("allow context must reach the socket"),
            Err(_) => {}
            Ok(_) => panic!("unroutable port must not succeed"),
        }
    }

    #[test]
    fn get_returns_status_headers_and_body() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"{\"x\":1}", Duration::ZERO, "");
        let cfg = TransportConfig::new().unwrap();
        let req = HttpRequest::get(&fixture.url)
            .unwrap()
            .header("X-Test", "yes");
        let resp = fetch(&req, "probe", allow(), &cfg).expect("loopback GET");
        assert_eq!(resp.status(), 200);
        assert_eq!(resp.body(), b"{\"x\":1}");
        assert_eq!(resp.body_text().unwrap(), "{\"x\":1}");
        let head = recv_head(&fixture);
        assert!(head.starts_with("GET /vectors HTTP/1.1\r\n"), "{head}");
        assert!(head.contains("x-test: yes\r\n"), "{head}");
        assert!(
            head.contains(&format!("user-agent: {}\r\n", DEFAULT_USER_AGENT)),
            "{head}"
        );
    }

    #[test]
    fn caller_user_agent_wins_over_default() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"ok", Duration::ZERO, "");
        let cfg = TransportConfig::new().unwrap();
        let req = HttpRequest::get(&fixture.url)
            .unwrap()
            .header("User-Agent", "libephemeris-download/1.0");
        fetch(&req, "probe", allow(), &cfg).expect("loopback GET");
        let head = recv_head(&fixture);
        assert!(
            head.contains("user-agent: libephemeris-download/1.0\r\n"),
            "{head}"
        );
        assert!(!head.contains(DEFAULT_USER_AGENT), "{head}");
    }

    #[test]
    fn post_sends_body_bytes() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"stored", Duration::ZERO, "");
        let cfg = TransportConfig::new().unwrap();
        let req = HttpRequest::post(&fixture.url, vec![1, 2, 3, 4]).unwrap();
        let resp = fetch(&req, "probe", allow(), &cfg).expect("loopback POST");
        assert_eq!(resp.status(), 200);
        assert_eq!(resp.body(), b"stored");
        let raw = fixture
            .received
            .recv_timeout(Duration::from_secs(10))
            .expect("server saw the request");
        let head = String::from_utf8_lossy(&raw).into_owned();
        assert!(head.starts_with("POST /vectors HTTP/1.1\r\n"), "{head}");
        assert!(raw.ends_with(&[1, 2, 3, 4]), "body bytes must arrive");
    }

    #[test]
    fn error_status_maps_to_typed_variant() {
        for (line, code) in [
            ("HTTP/1.1 400 Bad Request", 400),
            ("HTTP/1.1 404 Not Found", 404),
            ("HTTP/1.1 500 Internal Server Error", 500),
        ] {
            let fixture = serve_once(line, b"nope", Duration::ZERO, "");
            let cfg = TransportConfig::new().unwrap();
            let req = HttpRequest::get(&fixture.url).unwrap();
            match fetch(&req, "SBDB lookup", allow(), &cfg) {
                Err(Error::HttpStatus { status, purpose }) => {
                    assert_eq!(status, code);
                    assert_eq!(purpose, "SBDB lookup");
                }
                other => panic!("{code} must be HttpStatus, got {other:?}"),
            }
        }
    }

    #[test]
    fn over_bound_body_is_discarded_not_truncated() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"0123456789abcdef", Duration::ZERO, "");
        let cfg = TransportConfig::new()
            .unwrap()
            .with_max_body_bytes(8)
            .unwrap();
        let req = HttpRequest::get(&fixture.url).unwrap();
        match fetch(&req, "probe", allow(), &cfg) {
            Err(Error::HttpResponseTooLarge {
                purpose,
                limit_bytes,
            }) => {
                assert_eq!(purpose, "probe");
                assert_eq!(limit_bytes, 8);
            }
            other => panic!("over-bound body must refuse, got {other:?}"),
        }
        // A body of exactly the bound passes through intact.
        let exact = serve_once("HTTP/1.1 200 OK", b"01234567", Duration::ZERO, "");
        let exact_req = HttpRequest::get(&exact.url).unwrap();
        let resp = fetch(&exact_req, "probe", allow(), &cfg).expect("exact bound passes");
        assert_eq!(resp.body(), b"01234567");
    }

    #[test]
    fn streamed_response_owns_headers_chunks_and_idempotent_close() {
        let body: Vec<_> = (0..131073).map(|i| (i % 251) as u8).collect();
        let fixture = serve_once(
            "HTTP/1.1 200 OK",
            &body,
            Duration::ZERO,
            "X-Owner: stream\r\n",
        );
        let cfg = TransportConfig::new()
            .unwrap()
            .with_max_body_bytes(1)
            .unwrap();
        let request = HttpRequest::get(&fixture.url).unwrap();
        let response = open_stream(&request, "large asset", allow(), &cfg).unwrap();
        assert_eq!(response.status(), 200);
        assert!(
            response
                .headers()
                .iter()
                .any(|(name, value)| name == "x-owner" && value == "stream")
        );
        // Moving the owned reader keeps its TLS/timeout/body resource valid.
        let actual = std::thread::spawn(move || {
            let mut response = response;
            let mut chunk = [0u8; 65536];
            let mut actual = Vec::new();
            loop {
                let size = response.read(&mut chunk).unwrap();
                if size == 0 {
                    break;
                }
                actual.extend_from_slice(&chunk[..size]);
            }
            response.close();
            response.close();
            assert!(matches!(
                response.read(&mut chunk),
                Err(Error::HttpTransport { .. })
            ));
            actual
        })
        .join()
        .unwrap();
        assert_eq!(actual, body);
        assert!(recv_head(&fixture).starts_with("GET /vectors HTTP/1.1\r\n"));
        fixture._thread.join().unwrap();
    }

    #[test]
    fn stream_sealed_policy_precedes_connection() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let request =
            HttpRequest::get(&format!("http://{}/sealed", listener.local_addr().unwrap())).unwrap();
        assert!(matches!(
            open_stream(
                &request,
                "asset",
                EffectiveNetworkPolicy::Sealed,
                &TransportConfig::default()
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert_eq!(
            listener.accept().unwrap_err().kind(),
            std::io::ErrorKind::WouldBlock
        );
    }

    #[test]
    fn stream_timeout_is_per_io_while_fetch_retains_total_budget() {
        fn slow_chunks() -> (String, std::thread::JoinHandle<()>) {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let url = format!("http://{}/chunks", listener.local_addr().unwrap());
            let thread = std::thread::spawn(move || {
                let (mut stream, _) = listener.accept().unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(2)))
                    .unwrap();
                let mut header = Vec::new();
                let mut byte = [0];
                while !header.ends_with(b"\r\n\r\n") {
                    stream.read_exact(&mut byte).unwrap();
                    header.push(byte[0]);
                }
                stream
                    .write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 12\r\n\r\n")
                    .unwrap();
                for byte in b"abcdefghijkl" {
                    std::thread::sleep(Duration::from_millis(50));
                    if stream.write_all(&[*byte]).is_err() {
                        break;
                    }
                }
            });
            (url, thread)
        }
        let config = TransportConfig::with_timeout_secs(0.2).unwrap();
        let (url, worker) = slow_chunks();
        let mut response = open_stream(
            &HttpRequest::get(&url).unwrap(),
            "chunked asset",
            allow(),
            &config,
        )
        .unwrap();
        let mut bytes = Vec::new();
        let mut buffer = [0; 64];
        loop {
            let count = response.read(&mut buffer).unwrap();
            if count == 0 {
                break;
            }
            bytes.extend_from_slice(&buffer[..count]);
        }
        assert_eq!(bytes, b"abcdefghijkl");
        worker.join().unwrap();
        let (url, worker) = slow_chunks();
        assert!(matches!(
            fetch(
                &HttpRequest::get(&url).unwrap(),
                "bounded payload",
                allow(),
                &config
            ),
            Err(Error::HttpTimeout { .. })
        ));
        worker.join().unwrap();
    }

    #[test]
    fn slow_server_hits_the_overall_budget() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"late", Duration::from_millis(1500), "");
        let cfg = TransportConfig::with_timeout_secs(0.2).unwrap();
        let req = HttpRequest::get(&fixture.url).unwrap();
        match fetch(&req, "probe", allow(), &cfg) {
            Err(Error::HttpTimeout { purpose, seconds }) => {
                assert_eq!(purpose, "probe");
                assert!((seconds - 0.2).abs() < 1e-9, "{seconds}");
            }
            other => panic!("slow server must time out, got {other:?}"),
        }
    }

    #[test]
    fn refused_connection_is_transport_not_sealed() {
        // Nothing listens on port 9 (discard): an allowed policy reaches
        // the socket and reports the refusal as a transport failure.
        let req = HttpRequest::get("http://127.0.0.1:9/vectors").unwrap();
        let cfg = TransportConfig::new().unwrap();
        match fetch(&req, "probe", allow(), &cfg) {
            Err(Error::HttpTransport { purpose, detail }) => {
                assert_eq!(purpose, "probe");
                assert!(!detail.is_empty());
            }
            other => panic!("refused connection must be HttpTransport, got {other:?}"),
        }
    }
}
