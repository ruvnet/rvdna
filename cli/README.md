# rvdna

rvdna genomics CLI — agent harness over the rvdna engine

> **Health & Wellness** — Intake → triage → coordinate, with a knowledge MCP. Hard-codes "see a clinician" for anything clinical.
>
> Generated with [`create-agent-harness`](https://github.com/ruvnet/agent-harness-generator). WASM kernel, multi-host support, witness-signed releases.

## Install

```bash
npm install -g rvdna
rvdna init
rvdna doctor
```

## Agents

| Agent | Role |
|---|---|
| `intake` | Collects structured intake, flags red flags. |
| `triage` | Routes to the right resource, not a diagnosis. |
| `care-coordinator` | Organises logistics and reminders. |

This harness ships with the **claude-code** adapter.

## License

MIT
