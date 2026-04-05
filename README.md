- [Wearable & Omics Biomarkers of Aging (Consortium-Aligned)](#wearable--omics-biomarkers-of-aging-consortium-aligned)
  - [Executive Summary](#executive-summary)
  - [Purpose](#purpose)
  - [Scope & Sources](#scope--sources)
  - [Quick Start](#quick-start)
  - [Core Definitions](#core-definitions)
  - [Consortium Framework: Roles & Validation](#consortium-framework-roles--validation)
  - [Curated Table A - Wearable & Functional Digital Biomarkers](#curated-table-a---wearable--functional-digital-biomarkers)
  - [Curated Table B - Clinical Chemistry & Omics Biomarkers](#curated-table-b---clinical-chemistry--omics-biomarkers)
  - [Algorithms Appendix](#algorithms-appendix)
  - [License](#license)
  - [References](#references)

# Wearable & Omics Biomarkers of Aging (Consortium-Aligned)

---


## Purpose

- Provide a **concise, reference-grounded** overview of biomarkers of aging most relevant to **wearables/functional measures** and **omics**.  
- Organize each biomarker by **domain, modality, validation stage, evidence type, and use cases**.  

## Scope & Sources

This README includes *only* claims supported by the three Biomarkers of Aging Consortium manuscripts:

- [1](#ref-1) **Cell (2023)** - *Biomarkers of aging for the identification and evaluation of longevity interventions.*  
- [2](#ref-2) **Nature Medicine (2024)** - *Validation of biomarkers of aging.*  
- [3](#ref-3) **Nature Aging (2024)** - *Challenges and recommendations for the translation of biomarkers of aging.*  

---

## Quick Start

- Read **Core Definitions** and **Framework** for context.  
- Use **Table A** (wearables/functional) and **Table B** (clinical chemistry + omics).  
- See **Algorithms Appendix** for metric definitions.  

---

## Core Definitions

- **Chronological age:** Time since birth.  
- **Biological age:** State of age-dependent changes; can diverge from chronological age.  
- **Biomarker of aging (BoA):** Quantitative measure predicting biological age or its trajectory.  
- **Healthspan:** Period of life in good health without chronic disease/disability.  
- **Age deviation (AgeDev):** **Difference between biological and chronological age**.  

---

## Consortium Framework: Roles & Validation

**Roles (Consortium definitions):**  
- **Predictive validity:** Forecasts outcomes beyond age/sex (mortality, morbidity, function).  
- **Response:** Tracks change when interventions influence aging.  
- **Prognostic:** Predicts course/progression in diseased populations.  
- **Mechanistic biomarkers:** Provide insight into aging biology, though not always outcome-predictive.  

**Validation stages:**  
- **Analytical validity:** Reliability/reproducibility; standardized computation.  
- **Predictive validity:** Association with clinically meaningful outcomes.  
- **Clinical utility:** Evidence that biomarker-guided actions improve outcomes.  

*These role categories are adapted from FDA-BEST and explicitly mapped to aging by the Consortium (Cell 2023).*  
*Guardrail:* No BoA is a qualified **surrogate endpoint**.

---

## Curated Table A - Wearable & Functional Digital Biomarkers



| Biomarker | Domain | Modality / Setting | Examples | Algorithms / Methods | Role | Validation Status | Evidence Type | Notes / Use Cases | Refs |
|---|---|---|---|---|---|---|---|---|---|
| **Sleep patterns** (duration, efficiency, fragmentation, circadian regularity) | Aging & Mental | Wearable actigraphy ± HR | Research actigraphs; consumer watches | Epoch scoring, staging, longitudinal trends | Predictive (exploratory); candidate Response | Exploratory: emerging | Longitudinal | Relevant to resilience, stress, mood outcomes | [1](#ref-1), [3](#ref-3) |
| **Sleep fragmentation index** | Aging | Actigraphy | Research actigraphs | WASO, immobility bouts | Predictive (exploratory) | Exploratory | Cross-sectional; Longitudinal | Associates with frailty risk in older adults | [1](#ref-1) |
| **Daily activity volume (steps)** | Aging & Mental | Accelerometry | Watches, bands, smartphones | Step detection, compositional metrics | Predictive; candidate Response | Exploratory: emerging | Longitudinal | Declines signal frailty; complements energy/mood tracking | [1](#ref-1), [2](#ref-2) |
| **Activity intensity & sedentary time** | Aging & Mental | Accelerometry ± HR | HR-enabled watches | Intensity zones, dose-response | Predictive; candidate Response | Emerging | Cross-sectional; Longitudinal | Links to cardiorespiratory fitness, energy | [1](#ref-1), [2](#ref-2) |
| **Circadian rest–activity regularity** | Mental (and Aging) | 24-h actigraphy | Actigraphs, watches | Relative amplitude, interdaily stability | Predictive (exploratory) | Exploratory | Longitudinal | Disruption linked to mood/cognition | [1](#ref-1), [3](#ref-3) |
| **Resting heart rate (RHR)** | Aging & Mental | PPG/ECG | Watches, chest straps | Basal extraction, trend analysis | Predictive (exploratory); candidate Response | Exploratory: emerging | Longitudinal; Interventional (limited) | Stress/recovery load indicator | [1](#ref-1) |
| **Heart rate variability (HRV)** | Aging & Mental | ECG/PPG | ECG watches, straps | RMSSD, SDNN, LF/HF | Predictive (exploratory); candidate Response | Exploratory: emerging | Longitudinal; Interventional (limited) | Autonomic stress load; mental health linkage | [1](#ref-1), [2](#ref-2) |
| **Gait speed (4-6 m walk)** | Aging & Mental | Clinical ± IMU | Stopwatch; phone/watch IMU | Distance/time; stride variability | Predictive; Clinical outcome | Clinically meaningful; not surrogate | Cross-sectional; Longitudinal | Strong morbidity/mortality predictor | [1](#ref-1), [2](#ref-2) |
| **6-minute walk distance** | Aging & Mental | Clinical protocol | Timed 6-min walk | Total distance; optional HR/SpO₂ | Predictive; Clinical outcome | Clinically meaningful outcome | Cross-sectional; Interventional | Functional aging and rehab measure | [1](#ref-1) |
| **Postural stability / balance** | Aging | Digital IMU | Standing tasks | IMU sway metrics | Predictive (exploratory) | Exploratory | Cross-sectional | Decline: falls risk; screening use | [1](#ref-1) |

---

## Curated Table B - Clinical Chemistry & Omics Biomarkers


| Biomarker (family) | Type | Outcomes | Role | Validation Status | Evidence Type | Strengths | Limitations | Refs |
|---|---|---|---|---|---|---|---|---|
| **Clinical blood chemistry composites** | Multi-analyte (e.g., phenotypic age indices) | Survival, morbidity, function | Predictive; some Response | Established in research | Longitudinal; Interventional | Uses routine clinical tests; actionable | Overlaps with known risk factors | [1](#ref-1), [2](#ref-2) |
| **Telomere length** | DNA | Modest associations with aging outcomes | Predictive (mixed) | Variable, weaker than DNAm | Cross-sectional; Longitudinal | Senescence link; easy assay | High variability; poor sensitivity | [1](#ref-1), [2](#ref-2) |
| **DNAmAge (Horvath, Hannum)** | Epigenetic | Mortality, disease, health | Predictive | Established; standardization ongoing | Cross-sectional; Longitudinal | Reproducible; simple | Trained on age, less responsive | [1](#ref-1), [2](#ref-2) |
| **DNAm PhenoAge** | Epigenetic | Mortality, multimorbidity, frailty | Predictive | Established | Longitudinal | Captures physiology | Still correlative | [1](#ref-1), [2](#ref-2) |
| **DNAm GrimAge** | Epigenetic | Mortality, function, exposures | Predictive | Established | Longitudinal | Strong outcome prediction | Embeds risk exposures | [1](#ref-1), [2](#ref-2) |
| **Pace-of-aging DNAm (DunedinPACE/PoAm)** | Epigenetic | Decline, mortality | Predictive; some Response | Emerging; replicated | Longitudinal; Interventional (limited) | Sensitive to change rate | Needs broader validation | [2](#ref-2) |
| **IgG glycans** | Glycomics | Inflammation, risk | Predictive (emerging) | Emerging | Cross-sectional | Immune/inflammation insight | Influenced by immune state | [1](#ref-1) |
| **Proteomic clocks** | Proteomics | Function, mortality | Predictive (emerging) | Emerging | Cross-sectional; Longitudinal | Rich biological info | Cost, batch effects | [2](#ref-2) |
| **Transcriptomic clocks** | Transcriptomics | Condition-specific signals | Predictive (research-stage) | Early research | Cross-sectional | Dynamic regulation | Tissue-sensitive, variable | [2](#ref-2), [3](#ref-3) |

---


## Algorithms Appendix

- **RMSSD / SDNN:** Time-domain HRV metrics.  
- **LF/HF ratio:** Frequency-domain HRV.  
- **Relative Amplitude (RA):** Day vs night activity ratio.  
- **Interdaily Stability (IS):** 24-h rhythm consistency.  
- **4-6 m gait speed:** Distance ÷ time.  
- **6-minute walk distance (6MWD):** Total meters walked.  
- **Sleep efficiency:** Sleep time ÷ time in bed.  
- **WASO:** Wake minutes after sleep onset.  

---



## License

MIT

---

## References

<a id="ref-1"></a>1. **Moqri, M., Herzog, C., Poganik, J.R. et al.** *Biomarkers of aging for the identification and evaluation of longevity interventions.* **Cell** (2023). DOI: [10.1016/j.cell.2023.08.003](https://doi.org/10.1016/j.cell.2023.08.003)  
<a id="ref-2"></a>2. **Moqri, M., Herzog, C., Poganik, J.R. et al.** *Validation of biomarkers of aging.* **Nature Medicine** 30, 360–372 (2024). DOI: [10.1038/s41591-023-02784-9](https://doi.org/10.1038/s41591-023-02784-9)  
<a id="ref-3"></a>3. **Biomarkers of Aging Consortium, Herzog C.M.S., Goeminne L.J.E. et al.** *Challenges and recommendations for the translation of biomarkers of aging.* **Nature Aging** 4, 1372–1383 (2024). DOI: [10.1038/s43587-024-00683-3](https://doi.org/10.1038/s43587-024-00683-3)

