# Coercive Control Dashboard

Interactive reference on abuse patterns, digital threats, and recording consent law by state.

Built for Academy (Substack) — opener for the Digital Safety series.

---

## Local Preview

The dashboard fetches `data.json` at runtime. Must be served via HTTP, not opened as `file://`.

```bash
cd _Claude/Outputs/coercive-control-dashboard/
python3 -m http.server 8080
```

Then open: `http://localhost:8080`

---

## GitHub Pages Deploy

Two options — choose based on repo structure.

### Option A: Dedicated repo (recommended)

```bash
# Create a new repo: coercive-control-dashboard
git init
git add index.html data.json README.md
git commit -m "initial dashboard"
git remote add origin https://github.com/[username]/coercive-control-dashboard.git
git push -u origin main
```

Enable GitHub Pages: repo Settings → Pages → Source: `main` branch, `/ (root)` folder.

Live URL: `https://[username].github.io/coercive-control-dashboard/`

### Option B: Subfolder of existing repo

Place files in `/docs/coercive-control-dashboard/`. Enable GitHub Pages from `/docs` folder.

Live URL: `https://[username].github.io/[repo]/coercive-control-dashboard/`

---

## Substack Embed

Paste into a Substack post or page (HTML block):

```html
<iframe
  src="https://[username].github.io/coercive-control-dashboard/"
  width="100%"
  height="820"
  frameborder="0"
  style="border-radius:4px;border:none;"
  title="Coercive Control Resource Dashboard"
></iframe>
```

**Notes:**
- Substack's iframe support varies by plan. Test in a draft post before publishing.
- If embedding is blocked, link directly to the GitHub Pages URL instead.
- Height 820px fits the dashboard header + first section without a scroll; adjust down if the footer disclaimer clips.

---

## Updating State Law Data

`index.html` is data-agnostic — all state and resource data lives in `data.json`. When laws change:

1. Edit `data.json` — update the relevant `states[]` entry
2. Push to GitHub — GitHub Pages redeploys automatically (usually within 60 seconds)
3. Update `meta.updated` field in `data.json`

No changes to `index.html` required for data updates.

---

## Critical: Verify Before Publishing

**Hawaii HRS § 709-906** — Criminal coercive control statute is a 5-year pilot enacted 2021 (expires ~2026). Check renewal status before publishing. See `scope_note` in `data.json` for the flagged entry.

Search: `Hawaii HRS 709-906 pilot renewal 2026`

If expired and not renewed, update `data.json`:
```json
"coercive_control": {
  "status": "none",
  "statute": null,
  "year": null,
  "scope_note": "Pilot statute HRS § 709-906 expired 2026 without renewal as of [date]."
}
```

---

## File Inventory

| File | Purpose |
|------|---------|
| `index.html` | Single-file dashboard — all CSS and JS inline |
| `data.json` | All state law, resource, and checklist data |
| `README.md` | This file |

Vault source: `_References/coercive-control-resource.md`
Design spec: `_Claude/Outputs/2026-05-22_COERCIVE-CONTROL-DASHBOARD-SPEC.md`
