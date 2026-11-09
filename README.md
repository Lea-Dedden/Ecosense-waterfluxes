# Ecosense-waterfluxes

## Soil moisture quality control

Automated quality flagging for near-surface soil moisture time series at
15-minute resolution, following the framework of Dorigo et al. (2013).

Flags are cumulative — one observation may carry several codes, e.g. `"D04,D09.2"`.
Only observations with no flag receive `G`. Thresholds are site-specific and in accordance to sensor installation depth.

| Flag | Meaning | Test |
|------|---------|------|
| `G` | Good quality | no other flag |
| `C01` | SM below moisture threshold | `SM < 1` [% VWC] |
| `C02` | SM below second threshold | *optional* |
| `C03` | SM above saturation point | `SM > 60` [% VWC] |
| `D01` | Soil temperature too low | `SoilTemp < 1` [°C] |
| `D04` | SM rise without precipitation | SM rise in 15-min interval + previous 24-h rise + previous 24-h rain sum = 0 |
| `D06` | SM spike (positive/negative) | Dorigo et al. (2013) Eq. 4–6 |
| `D07` | SM drop / jump | Dorigo et al. (2013) Eq. 7–9 |
| `D09.1` | SM invariant values ≥ 12 h | rolling variance `== 0` |
| `D09.2` | SM identical value ≥ 48 h | rolling variance `== 0` |

## Throughfall & stemflow quality control
tbd

## Sapflow quality control & baseline corrections
tbd

## Requirements

R ≥ 4.2 and:
install.packages(c("dplyr", "purrr", "stringr", "zoo", "yaml", "readr"))

## Citation
Please cite Dorigo et al. (2013) alongside the code for soil moisture quality control.
Please cite Seeger et al. (2026) alongside the code for sapflow density calculation.
#tbd own publications

## License
This code can be distributed freely under the GPL v3 license. Please read the LICENSE for further information.

© 2026, Lea Dedden
