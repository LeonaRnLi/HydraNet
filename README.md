# From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs
HydraNet can improve spatial reuse and optimize spectrum efficiency in LPWANs.

This webpage contains source code to deploy HydraNet.


## Overview
Low Power Wide-Area Networks (LPWANs) suffer from sever interference due to their long range transmissions. Traditional spatial reuse solutions share a common philosophy, i.e., interference mitigation. 

HydraNet is a new paradigm to facilitate simultaneous LPWAN transmissions on the same frequency,fundamentally evolving the current philosophy of interference management—from mitigating interference to tolerating its presence.

### Reference:
Ruonan Li, Ziyue Zhang, Xianjin Xia, et al. From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs, ACM MobiCom 2025. 

The PDF file in the main folder.
## Hardware Requirements and Dependencies
- **Hardware**
LoRa nodes: SX1276 + Arduino Uno (version: 1.8.X)

Open source gateway: External + USRP N210 with more than 1 antenna (or other Software Defined Radios with equivalent functionality) + GNU Radio Platform

Data trace processing: A workstation with Matlab (version: R2022a)
- This code is tested on Ubuntu 18.04.
## HydraNet Spatial Reuse Protocol
Below is the protocol of HydraNet, illustrating how it achieve spatial reuse in LPWANs.

## To Run the Artifact 

## To Run the HydraNet on Your Testbed

## Video Demo

## Contact & Citation
For questions, open an issue or contact:

- **Ruonan Li** - ruo-nan.li@connect.polyu.hk  
- **Ziyue Zhang** - ziyue.zhang@connect.polyu.hk  

If you use HydraNet in your research, please cite:

```bibtex
@inproceedings{hydranet2025,
  author    = {Ruonan Li and Ziyue Zhang and Xianjin Xia and Ningning Hou and 
               Wenchang Chai and Shiming Yu and Yuanqing Zheng and Tao Gu},
  title     = {From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs},
  booktitle = {Proceedings of the 31st Annual International Conference on Mobile Computing and Networking (ACM MobiCom ’25)},
  year      = {2025},
  publisher = {ACM},
  location  = {Hong Kong, China},
  doi       = {10.1145/3680207.3723483}
}
