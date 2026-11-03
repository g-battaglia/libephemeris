//! Closed primary selector and native-array transport enums.

/// Closed Calendar contract. Adding a variant requires a coordinated major version.
#[derive(Clone, Debug, PartialEq, Copy, Eq)]
pub enum Calendar {
    /// Julian selector or actual native transport.
    Julian,
    /// Gregorian selector or actual native transport.
    Gregorian,
}

/// Closed Target contract. Adding a variant requires a coordinated major version.
#[derive(Clone, Debug, PartialEq)]
pub enum Target {
    /// Body selector or actual native transport.
    Body(i32),
    /// Star selector or actual native transport.
    Star(String),
}

/// Closed LocalSolarTimes contract. Adding a variant requires a coordinated major version.
#[derive(Clone, Debug, PartialEq)]
pub enum LocalSolarTimes {
    /// Seven selector or actual native transport.
    Seven([f64; 7]),
    /// Ten selector or actual native transport.
    Ten([f64; 10]),
}

/// Closed RiseTransitTimes contract. Adding a variant requires a coordinated major version.
#[derive(Clone, Debug, PartialEq)]
pub enum RiseTransitTimes {
    /// One selector or actual native transport.
    One([f64; 1]),
    /// Ten selector or actual native transport.
    Ten([f64; 10]),
    /// Absent selector or actual native transport.
    Absent,
}

/// Closed ErrorKind contract. Adding a variant requires a coordinated major version.
#[derive(Clone, Debug, PartialEq, Copy, Eq)]
pub enum ErrorKind {
    /// InvalidInput selector or actual native transport.
    InvalidInput,
    /// MissingData selector or actual native transport.
    MissingData,
    /// OutOfRange selector or actual native transport.
    OutOfRange,
    /// Calculation selector or actual native transport.
    Calculation,
    /// Configuration selector or actual native transport.
    Configuration,
    /// Unsupported selector or actual native transport.
    Unsupported,
    /// Unavailable selector or actual native transport.
    Unavailable,
}

impl Calendar {
    /// Decode only the Julian/Gregorian zero/one selectors; all other values fail.
    pub fn from_flag(flag: i32) -> Result<Self, super::Error> {
        match flag {
            0 => Ok(Self::Julian),
            1 => Ok(Self::Gregorian),
            value => Err(super::Error::from_engine(
                "Calendar::from_flag",
                crate::errors::Error::InvalidCalendar { value },
            )),
        }
    }
    /// Return the original zero/one calendar selector.
    pub fn flag(self) -> i32 {
        match self {
            Self::Julian => 0,
            Self::Gregorian => 1,
        }
    }
    pub(crate) fn native(self) -> crate::time::Calendar {
        match self {
            Self::Julian => crate::time::Calendar::Julian,
            Self::Gregorian => crate::time::Calendar::Gregorian,
        }
    }
}
