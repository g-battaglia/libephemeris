// SPDX-License-Identifier: AGPL-3.0
//! Explicit process boundary for the optional standalone offline fitter.
//!
//! Runtime astronomy stays native. This command starts the separately installed
//! fitting application with literal arguments; it never embeds Python, installs
//! dependencies or falls back to a source checkout.

use std::path::Path;
use std::process::Command;

pub(super) const EXECUTABLE: &str = "libephemeris-generate";

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Request {
    arguments: Vec<String>,
}

pub(super) fn parse<'a>(args: impl Iterator<Item = &'a str>) -> Result<Request, String> {
    let arguments: Vec<String> = args.map(str::to_owned).collect();
    let tier = arguments
        .first()
        .ok_or("leb generate needs a tier (base|medium|extended|all)")?;
    let options_start = if tier == "all" {
        1
    } else {
        if !matches!(tier.as_str(), "base" | "medium" | "extended") {
            return Err(format!(
                "unknown leb tier '{tier}' (base|medium|extended|all)"
            ));
        }
        let mode = arguments.get(1).ok_or("leb generate <tier> needs a mode (groups|full|single|body <name>|planets|asteroids|exotics|analytical|merge)")?;
        match mode.as_str() {
            "body" => {
                if arguments.get(2).is_none_or(|value| value.starts_with('-')) {
                    return Err("leb generate <tier> body needs a body name or id".into());
                }
                3
            }
            "groups" | "full" | "single" | "planets" | "asteroids" | "exotics" | "analytical"
            | "merge" | "core" | "apogee" | "--help" | "-h" => 2,
            _ => return Err(format!("unknown leb generate mode '{mode}'")),
        }
    };
    if let Some(extra) = arguments
        .get(options_start)
        .filter(|value| !value.starts_with('-'))
    {
        return Err(format!("unexpected argument '{extra}'"));
    }
    Ok(Request { arguments })
}

fn command(request: Request, program: &Path) -> Command {
    let mut command = Command::new(program);
    command.args(request.arguments);
    command
}

pub(super) fn run(request: Request, program: Option<&Path>) -> i32 {
    let program = program.unwrap_or_else(|| Path::new(EXECUTABLE));
    match command(request, program).status() {
        Ok(status) => status.code().unwrap_or(super::EXIT_FAILED),
        Err(error) => {
            eprintln!(
                "libephemeris: cannot start offline LEB generator '{}': {error}. Install libephemeris-leb-tools 4.0.0b1 in the same environment and put its executable on PATH for the native binary",
                program.display()
            );
            super::EXIT_FAILED
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arguments_cross_the_process_boundary_literally() {
        let words = [
            "base",
            "body",
            "sun,MOON,14",
            "--output",
            "a space;$(false).leb",
        ];
        let request = parse(words.into_iter()).unwrap();
        let child = command(request, Path::new("a directory/libephemeris-generate"));
        assert_eq!(child.get_program(), "a directory/libephemeris-generate");
        assert_eq!(child.get_args().collect::<Vec<_>>(), words);
    }

    #[cfg(unix)]
    #[test]
    fn normal_child_failure_status_is_preserved() {
        let request = Request {
            arguments: vec!["-c".into(), "exit 23".into()],
        };
        assert_eq!(run(request, Some(Path::new("/bin/sh"))), 23);
    }

    #[test]
    fn missing_explicit_generator_is_an_operational_failure() {
        let request = parse(["base", "full"].into_iter()).unwrap();
        assert_eq!(
            run(
                request,
                Some(Path::new("/no-such-libephemeris-directory/generate"))
            ),
            1
        );
    }
}
