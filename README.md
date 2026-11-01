# MIDSSAR: MId Dispersion Spectrograph for Stellar Activity Research
Reduction module for MIDSSAR on Seimei telescope:
[https://www.kusastro.kyoto-u.ac.jp/~iwamuro/CaHa/index.html](https://www.kusastro.kyoto-u.ac.jp/~iwamuro/CaHa/index.html)


## Installation

```bash
pip install git+https://github.com/KaiIkuta/MIDSSAR.git
```

## Import (To be developed soon)

```python
import mzr
mzr.mzr_reduction(obs_dates=["20260425"],observer=["Ikuta"]).run_all()
```
### Memorandum
If skyflat was not taken, it should be downloaded from our folder (or copied from your recent analysis).

Calibration should be analysed automatically (easily), and the target should be designated from observer name.




## Analysis procedure (Japanese): 

[https://www.kusastro.kyoto-u.ac.jp/~iwamuro/MIDSSAR/index.html](https://www.kusastro.kyoto-u.ac.jp/~iwamuro/MIDSSAR/index.html)

## Paper

[Iwamuro et al., PASJ in press](https://arxiv.org/abs/2610.22834)
