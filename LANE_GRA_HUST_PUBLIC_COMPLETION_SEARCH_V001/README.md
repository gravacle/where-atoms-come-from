# HUST-2018 public completion search V001

Date: 2026-08-27  
Status: **PASS — bounded search completed; no qualifying public completion root located**

## Purpose

This lane asks one narrow question left by the calibrated-source identifiability theorem: can a primary or authoritative public object supply the HUST-2018 physical-harmonic remainder, the physical maps needed to reconstruct it, or the campaign raw/correction/covariance packet needed for an independently owned (G) cross-check?

It does not repeat broad gravity-data searches and does not alter the calibrated source model. It searches the official article/release surface, DOI registries, major research-data repositories, Chinese scientific-data surfaces, HUST institutional surfaces, code hosts and campaign-specific thesis catalogues.

## Result

No object meeting the predeclared completion-root criteria was located on the bounded surfaces searched on 2026-08-27. The official Nature record still exposes seven associated public objects—the supplement, one supplementary-data workbook, and five source-data workbooks—and states that additional supporting data are available from the corresponding authors on reasonable request. Crossref exposes no related dataset; exact DOI searches in DataCite, Figshare and Zenodo produced either zero results or later works that merely cite the HUST paper.

Two campaign-specific acquisition leads remain useful, but they have different evidentiary types:

- Jun-Fei Wu, *Improved experiment for measuring the Newtonian gravitational constant G with the angular-acceleration-feedback method* (confirmed HUST doctoral dissertation lead, 2021; Chinese title in the search ledger);
- *Measurement of the Newtonian gravitational constant G by the time-of-swing method using high-Q silica fibres* (unverified title-only bibliographic lead). Its discoverable metadata do not authenticate the author, institution, year or document type.

No openly inspectable payload was located, so neither lead is promoted to evidence that the missing numeric maps or raw packet are present.

## Scientific consequence

The public-data (G) cross-check remains exactly where the parent theorem placed it:

\[
G_i^{\rm AAF}(r_i)=\frac{\alpha_i f_{m,i}}{K_i^{\rm partial}+r_i},
\qquad
G_i^{\rm ToS}(r_i)=\frac{\Delta\omega_i^2(1+c_i^{\rm anel}+c_i^{\rm mag})}{K_i^{\rm partial}+r_i}.
\]

One independently owned normalized remainder (r_i) per released row would close row-wise point evaluation. The underlying density/CMM/orientation/attachment/deformation maps would close a fuller physical-source reconstruction. Campaign samples plus the correction/event/operator and design/covariance packet would be needed for an independent raw-response refit.

## Exact ceiling

This is a dated, bounded search result—not a theorem that the data do not exist, are irrecoverable, or are absent from every public surface. Its API-backed endpoint checks are reproducible; its aggregate index-mediated entries are documented curator searches, not executable completeness certificates. The next high-yield step is targeted acquisition from the corresponding authors, lawful access to the confirmed dissertation, and authoritative resolution of the title-only lead, followed by payload inspection against the field list in `RESULT.json`.

Run:

```bash
python3 verify_hust_public_completion_search.py
```
