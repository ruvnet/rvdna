# ADR-001: metaharness-generated agent-harness CLI for rvdna

- Status: accepted
- Date: 2026-06-17

## Context

rvdna was extracted from `ruvnet/ruvector` (ADR-257) and now builds as a
standalone Rust workspace: a quantum/AI genomics engine plus the `.rvdna`
cognitive-container format. The Rust crates expose the engine, but there was no
agent-facing entry point — no CLI an LLM host (Claude Code) could drive to boot
the engine, verify its environment, or coordinate agents over it.

rUv's `npx metaharness` generator scaffolds exactly this: a non-interactive,
signed agent-harness CLI built on `@metaharness/kernel` (a cross-platform kernel
with a `native | wasm | js` backend) and a host adapter
(`@metaharness/host-claude-code`). We want that surface for rvdna without
disturbing the Rust workspace.

## Decision

Ship a metaharness-generated agent-harness CLI for rvdna, targeting the kernel's
WASM backend where available.

- The CLI lives in `cli/` (a subdirectory, not the repo root) so it stays out of
  the Rust workspace and cannot clash with `Cargo.toml` / crate builds.
- Generated with: `npx metaharness rvdna --target cli
  --template vertical:health --host claude-code --force`.
- Commands: `init` (boot the kernel + host adapter, report status) and `doctor`
  (verify the install end-to-end and report which backend answered).
- Runtime dependencies: `@metaharness/kernel` and
  `@metaharness/host-claude-code`.
- `bin/cli.js` is plain ESM and runs with no build step; `npm run build` (tsc)
  is only needed when extending the TypeScript under `src/`.
- The kernel selects its backend at load time in the order `native > wasm > js`.
  We prefer wasm; the pure-JS backend is the guaranteed floor.

## Consequences

- The kernel and host adapter are beta `0.1.x` packages and are an **optional**
  surface — additive to, and independent of, the Rust crates. Nothing in the
  existing workspace depends on them; removing `cli/` leaves rvdna unchanged.
- The beta `0.1.0` kernel published to npm ships only the pure-JS backend: the
  wasm `pkg/` artifact and the per-platform NAPI-RS native binaries are built by
  separate CI jobs and are absent from a plain `npm install`. As a result
  `doctor` currently reports the `js` backend on a clean install. When the kernel
  publishes its wasm/native artifacts, the same CLI will transparently upgrade to
  `wasm` (or `native`) with no code change, because backend selection is dynamic.
- `cli/node_modules` and `cli/dist` are build/install artifacts and are
  git-ignored.
- The CLI's version line is pinned to a beta range (`^0.1.0`); expect churn until
  the kernel reaches a stable release.
