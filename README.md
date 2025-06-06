# FM + RDS Decoder using Python + RTL-SDR

This project decodes FM radio and extracts **RDS (Radio Data System)** information from IQ samples or live RTL-SDR input. Built with **Python**, **NumPy**, and **signal processing** techniques.

---

## Project Files

| File | Description |
|------|-------------|
| `main.py` | Main decoder pipeline |
| `fm-receiver-rds-packets.py` | RDS packet extraction |
| `clean_FM_RDS_Info_Basic.py` | Simplified version |
| `FM_RDS_Thonny_v1.py` | Older version (Thonny IDE) |
| `fm.wav` | Sample FM radio IQ data |
| `rds_log.txt` | Log of decoded RDS output |

---

## Requirements

- Python 3.x
- NumPy
- SciPy
- (Optional) `pyrtlsdr` for RTL-SDR support

### Install

```bash
pip install numpy scipy pyrtlsdr
```


## Features

- FM Demodulation
- Stereo support
- RDS extraction at 57 kHz
- Group decoding (PS, RT, CT)

## Author
Soham Samantaray
GitHub: CoalitionStorm

## License
This project is open-source and intended for educational use.
