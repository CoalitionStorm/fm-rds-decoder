# FM + RDS Decoder using Python + RTL-SDR

This project decodes FM radio and extracts **RDS (Radio Data System)** information from IQ samples or live SDR input. Built using **Python**, **NumPy**, and **signal processing** techniques.

---

## Project Files

| File | Description |
|------|-------------|
| `main.py` | Main decoder script |
| `fm-receiver-rds-packets.py` | RDS packet extraction |
| `clean_FM_RDS_Info_Basic.py` | Simplified decoder |
| `FM_RDS_Thonny_v1.py` | Thonny IDE version |
| `fm.wav` | Sample FM IQ data |
| `rds_log.txt` | Decoded output log |

---

## Requirements

- Python 3.x
- NumPy
- SciPy
- (Optional) pyrtlsdr for RTL-SDR support

Install:
```bash
pip install numpy scipy pyrtlsdr

## How to Use

To run with a sample `.wav` file:
```bash
python main.py

## How to Use

To run with a sample `.wav` file:
```bash
python main.py

## Author

Soham Samantaray  
GitHub: [CoalitionStorm](https://github.com/CoalitionStorm)

## License

This project is open-source and intended for educational use.