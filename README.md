# CSI-RFF
This repository contains the reference code for the article [''CSI-RFF: Leveraging Micro-Signals on CSI for RF Fingerprinting of Commodity WiFi''](https://ieeexplore.ieee.org/document/10517677). 

If you find the project useful and you use this code, please cite our articles:
```
@article{ruiqi2024csirff,
  author={Kong, Ruiqi and Chen, He},
  journal={IEEE Transactions on Information Forensics and Security}, 
  title={CSI-RFF: Leveraging Micro-Signals on CSI for RF Fingerprinting of Commodity WiFi}, 
  year={2024},
  volume={19},
  number={},
  pages={5301-5315},
  doi={10.1109/TIFS.2024.3396375}}
```
```
@INPROCEEDINGS{ruiqi2023phy,
  author={Kong, Ruiqi and Chen, He Henry},
  booktitle={2023 IEEE 24th International Workshop on Signal Processing Advances in Wireless Communications (SPAWC)}, 
  title={Physical-Layer Authentication of Commodity Wi-Fi Devices via Micro-Signals on CSI Curves}, 
  year={2023},
  volume={},
  number={},
  pages={486-490},
  doi={10.1109/SPAWC53906.2023.10304542}}
```

## Usage
* `main.py` for authentication using forest, lof, ocsvm, knn algorithms.
* `DBSCAN_main.m` for authentication using the DBSCAN algorithm.
* `Fingerprint.m` for micro-csi extraction.
* `adr_calculate.m` for attack detection rate calculation.
  
## Contact
kr020@ie.cuhk.edu.hk
