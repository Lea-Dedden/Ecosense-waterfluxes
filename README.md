# Ecosense-waterfluxes

## soil moisture quality control

Automated quality flagging for near-surface soil moisture time series at
15-minute resolution, following the framework of Dorigo et al. (2013).

Flags are cumulative — one observation may carry several codes, e.g. `"D04,D09.2"`.
Only observations with no flag receive `G`.

| Flag | Meaning | Test |
|------|---------|------|
| `G` | Good quality | no other flag |
| `C01` | Below moisture threshold | `SM < 1` |
| `C02` | Below second threshold | *optional* |
| `C03` | Above saturation point | `SM > 60` |
| `D01` | Soil temperature too low | `SoilTemp < 1` °C |
| `D04` | Rise without precipitation | 15-min rise + 24-h rise + 24-h rain sum |
| `D06` | Spike (positive/negative) | Dorigo et al. (2013) Eq. 4–6 |
| `D07` | Drop / jump | Dorigo et al. (2013) Eq. 7–9 |
| `D09.1` | Invariant values ≥ 12 h | rolling variance `== 0` |
| `D09.2` | Identical value ≥ 48 h | constant run ≥ X steps |

## throughfall & stemflow quality control
tbd

## sapflow quality control & baseline corrections
tbd

## Requirements

R ≥ 4.2 and:
install.packages(c("dplyr", "purrr", "stringr", "zoo", "yaml", "readr"))

## Citation
Please cite Dorigo et al. (2013) alongside the code for soil moisture quality control.
Please cite Seeger et al. (2026) alongside the code for sapflow density calculation.
#tbd own publications

## License
tbd
