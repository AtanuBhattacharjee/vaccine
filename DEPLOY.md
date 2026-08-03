# Deploying to Posit Connect Cloud

The app is one file, `app.R`, and needs only **shiny** and **shinydashboard**.
All the numerical work is base R — no Rcpp, no compiled dependencies, nothing
that will fail to build on the server.

---

## The short version

1. Put `app.R` and `manifest.json` in a **public GitHub repository**.
2. Go to <https://connect.posit.cloud>, sign in with GitHub.
3. **Publish → Shiny → R**, pick the repository, branch and `app.R`.
4. Publish. First build takes a few minutes while packages install.

---

## Regenerate the manifest first

`manifest.json` is included, but it was written on a machine whose R is 4.3.3
and whose `shiny` is 1.8.0 — versions that came from an Ubuntu archive rather
than CRAN. Connect Cloud pins to exactly what the manifest says, so it is worth
one line on your own machine before you publish:

```r
install.packages("rsconnect")
rsconnect::writeManifest(appDir = "path/to/app-folder")
```

Commit the regenerated `manifest.json` alongside `app.R`. If a build fails with
a package-resolution error, a stale manifest is almost always the cause.

Note that `writeManifest()` refuses to run if the folder contains both `app.R`
and a separate single-file Shiny app — keep the folder to `app.R`,
`manifest.json` and `DEPLOY.md` only.

---

## Alternative: publish straight from RStudio

If you would rather not go through GitHub:

```r
install.packages("rsconnect")
rsconnect::deployApp(appDir = "path/to/app-folder", server = "shinyapps.io")
```

This targets shinyapps.io rather than Connect Cloud. The free tier gives 25
active hours a month, which is ample for a paper companion but will suspend the
app if it is heavily used.

---

## Performance, and one thing to watch

The φ slider is the expensive control. With competition off, a full refresh
takes about 0.2 s. With competition on it takes roughly 1.3 s, because the
objective then has to be integrated over the panel-size distribution as well as
over immune capacity and epitope quality. A progress bar covers it.

On the free Connect Cloud tier the app runs on a single modest instance, so
several people moving the φ slider at once will queue. If that becomes a
problem, the cheapest fix is to drop the competition quadrature from 24 × 10
nodes to 16 × 8 — search for `GHq <- gh_nodes(10)` in `app.R`. The separation
changes in the second decimal place only.

The sensitivity table is behind a button rather than reactive on purpose: it
solves the model 25 times, which would otherwise fire on every slider nudge.

---

## Checks run before shipping

- `app.R` parses, the UI renders (23 kB of HTML), and every input and output id
  resolves.
- The server was exercised headlessly with `shiny::testServer`: all four value
  boxes and all three plots render, and the reset button and every slider
  update correctly.
- The model reproduces the manuscript figures: d\*_B = 396 µg and d\*_A = 444 µg
  at the reference configuration (12.1 %); exactly 0.0 % at τ = 0 with κ = 0;
  29.2 % at τ = 1; 20.3 % with panel means 24 and 8; and 396 / 451 µg (13.8 %)
  with competition at φ = 0.5. These match `epiOBD_sensitivity.R` to the second
  decimal place.

---

## Files

| file | purpose |
|---|---|
| `app.R` | the whole application — model, UI, server |
| `manifest.json` | dependency lock for Connect Cloud; regenerate before publishing |
| `DEPLOY.md` | this file |
