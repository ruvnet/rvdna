# Embedded typefaces

Both faces are variable-weight woff2 (latin subset) and are inlined as data URIs
by `../build.mjs`, because the Artifact CSP blocks font CDNs and a silent
fallback would lose the typography the Cognitum design system is built on.

| File | Family | Licence |
|---|---|---|
| `Outfit.woff2` | Outfit — display and body | [SIL Open Font License 1.1](https://openfontlicense.org) |
| `JetBrainsMono.woff2` | JetBrains Mono — data, labels, code | [SIL Open Font License 1.1](https://openfontlicense.org) |

The OFL permits embedding. Neither file is modified.
